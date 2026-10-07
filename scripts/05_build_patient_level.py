#!/usr/bin/env python3
"""构建患者级纳排表，作为后续所有分析的唯一输入。

设计原则（来自 PROJECT.md 执行纪律）：
- 分析单位是**患者**，不是样本；重复活检进入配对/混合结构，不当独立样本计数。
- 每行一个患者，携带纳排标记与排除理由，排除必须写明依据。
- 所有"重复/冲突"的判定都来自 2026-10-06 的回原文核对结论，不在此处重新推断。

输入（均在本仓库 meta/ 下）：
  GSE91061_samples.tsv, GSE91061_Riaz2017_sample_mapping.tsv
  GSE78220_samples.tsv, GSE244982_samples.tsv, GSE294272_samples.tsv, GSE308434_samples.tsv

输出：
  meta/patient_level.tsv        患者级纳排表（主表）
  meta/patient_level_summary.md 各队列纳入/排除统计与口径说明

用法：
  python3 scripts/05_build_patient_level.py --outdir /path/to/project
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

FIELDS = [
    "cohort", "patient_id", "role", "n_samples",
    "gsm_list", "sample_titles", "visits",
    "label_primary", "label_source", "key_covariates",
    "dup_flag", "dup_rule",
    "include_primary", "include_sensitivity", "exclude_reason", "notes",
]

# 作者临床表 IBOR 的原始写法 → 与样本表 BOR 同一套三分类口径。
# 临床表该字段大量缺失，且与 BOR 有冲突，只作回退与冲突标记，不作主口径。
IBOR_MAP = {
    "COMPLETE RESPONSE": "PRCR", "CR": "PRCR",
    "PARTIAL RESPONSE": "PRCR", "PR": "PRCR",
    "STABLE DISEASE": "SD", "SD": "SD",
    "PROGRESSION": "PD", "PD": "PD",
    "DEATH PRIOR TO DISEASE ASSESSMENT": "NE",
    "DEATH": "NE",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def char_val(row: dict[str, str], key: str, default: str = "") -> str:
    """从 characteristics 列里取值。列名为 'char::<key>'，两侧都做归一化比对。"""
    want = key.strip().lower()
    for k, v in row.items():
        name = (k or "").strip().lower()
        if name.startswith("char::"):
            name = name[len("char::"):]
        if name == want:
            return (v or "").strip()
    return default


def join(items) -> str:
    return "|".join(str(x) for x in items)


def build_gse91061(meta: Path) -> list[dict]:
    """主开发队列。109 样本 / 65 患者。
    核对结论：出处 Riaz 2017 Cell；nivolumab 单药；标签 = RECIST v1.1 BOR；
    PRCR = CR+PR；Pre/On 为同一病灶配对活检。
    """
    samples = read_tsv(meta / "GSE91061_samples.tsv")
    mapping = {
        r["geo_title"]: r
        for r in read_tsv(meta / "GSE91061_Riaz2017_sample_mapping.tsv")
    }

    by_pt: dict[str, list[dict]] = defaultdict(list)
    for s in samples:
        m = re.match(r"^(Pt\d+)_(Pre|On)", s["title"])
        if not m:
            raise ValueError(f"GSE91061 标题无法解析患者号：{s['title']}")
        by_pt[m.group(1)].append(s)

    out = []
    for pt in sorted(by_pt, key=lambda x: int(x[2:])):
        rows = by_pt[pt]
        titles = [r["title"] for r in rows]
        maps = [mapping[t] for t in titles if t in mapping]
        dup, rule, notes = "", "", ""

        # BOR 是患者级属性，同一患者的 Pre/On 应一致。
        # 作者表未收录的分装样本（NOT_IN_TABLE）不参与一致性判定，改用同患者的已知值。
        known = [m for m in maps if m["riaz_bor"] != "NOT_IN_TABLE"]
        bors = {m["riaz_bor"] for m in known}
        if len(bors) > 1:
            raise ValueError(f"{pt} 的 BOR 不一致：{bors}")
        bor = bors.pop() if bors else ""
        label = {"CR": "PRCR", "PR": "PRCR", "SD": "SD", "PD": "PD", "NE": "NE"}.get(bor, bor)

        # 标签来源与冲突记录：
        # BOR（作者样本表）是与 GEO response 字段 100% 对账的主口径；
        # IBOR（作者临床表）大量缺失且与 BOR 有 4 例冲突，只作回退与冲突标记，不作主口径。
        ibor_raw = {m["clinical_ibor"] for m in maps if m["clinical_ibor"]}
        ibor = IBOR_MAP.get(sorted(ibor_raw)[0].upper(), "") if ibor_raw else ""
        label_src = f"Riaz 2017 RECIST v1.1 BOR={bor}" if bor else ""
        if not bor and ibor:
            label = ibor
            label_src = f"⚠ 作者样本表未收录该样本，标签由临床表 IBOR={sorted(ibor_raw)[0]} 回填"
            dup = dup or "仅有 1 个样本且未收录于作者样本表"
            rule = rule or "标签为单一来源回填，须在敏感性分析中剔除复核"
            notes = (notes + "；" if notes else "") + "管号冲突样本（作者表把同一管号登记给 Pt7_Pre），标签来源单一"
        elif bor and ibor and ibor != label:
            notes = (notes + "；" if notes else "") + (
                f"⚠ 标签来源冲突：样本表 BOR={bor} vs 临床表 IBOR={sorted(ibor_raw)[0]}；"
                f"主分析采用 BOR（与 GEO response 对账），敏感性分析改用 IBOR 复核"
            )
        label_source = f"{label_src}（PMID 29033130；PRCR=CR+PR 合并）" if label_src else "标签不可得"

        visits = [re.match(r"^Pt\d+_(Pre|On)", t).group(1) for t in titles]
        prior_ici = {m["prior_ici"] for m in maps}
        pop = {m["pop_categ"] for m in maps}
        os_days = {m["clinical_os_days"] for m in maps}
        pfs = {m["clinical_ibor"] for m in maps}

        dup2, rule2, notes2 = dup, rule, notes
        if len(rows) > 2 or (len(rows) == 2 and visits.count("On") == 2):
            dup2 = "同一活检两个分装样本"
            rule2 = "合并或择一（默认保留 -5）；不得当两个独立样本"
            notes2 = "作者样本表未收录 -6，以 GEO 标题为准"
        elif len(rows) == 2:
            dup2 = "Pre/On 同一病灶配对"
            rule2 = "按患者内配对结构处理，不当独立样本"
        unmatched = [t for t in titles if t not in mapping]
        if unmatched and dup2:
            notes2 = (notes2 + "；" if notes2 else "") + f"作者表未收录：{join(unmatched)}"
        dup, rule, notes = dup2, rule2, notes2

        if label == "NE":
            # 疗效不可评估，既无应答也无进展，任何以疗效为终点的分析都用不上
            inc_pri, inc_sen, why = "否", "否", "RECIST NE（疗效不可评估），主分析与敏感性分析均无法归类"
        elif label == "SD":
            inc_pri, inc_sen, why = "否", "是", "SD 不进主分析（主终点为 PRCR vs PD），单列做敏感性分析"
        elif label in ("PRCR", "PD"):
            inc_pri, inc_sen, why = "是", "是", ""
        else:
            # 兜底：任何无法归入既定口径的标签一律不进主分析
            inc_pri, inc_sen, why = "否", "否", f"标签无法归类（label={label!r}），主分析与敏感性分析均排除"

        out.append({
            "cohort": "GSE91061", "patient_id": pt, "role": "主开发队列",
            "n_samples": len(rows),
            "gsm_list": join(r["gsm"] for r in rows),
            "sample_titles": join(titles),
            "visits": join(visits),
            "label_primary": label,
            "label_source": label_source,
            "key_covariates": join([
                f"prior_ici={join(sorted(prior_ici))}",
                f"pop_categ={join(sorted(pop))}",
                f"OS_days={join(sorted(os_days))}",
                f"IBOR={join(sorted(pfs))}",
            ]),
            "dup_flag": dup, "dup_rule": rule,
            "include_primary": inc_pri, "include_sensitivity": inc_sen,
            "exclude_reason": why, "notes": notes,
        })
    return out


def build_gse78220(meta: Path) -> list[dict]:
    """外部验证队列。28 样本 / 27 患者。
    核对结论：Pt27A/B 为同患者治疗前两个不同解剖部位活检；Pt16 仅 on-treatment。
    """
    samples = read_tsv(meta / "GSE78220_samples.tsv")
    by_pt: dict[str, list[dict]] = defaultdict(list)
    for s in samples:
        by_pt[char_val(s, "patient id")].append(s)

    out = []
    for pt in sorted(by_pt, key=lambda x: int(x[2:])):
        rows = by_pt[pt]
        titles = [r["title"] for r in rows]
        resp = {char_val(r, "anti-pd-1 response") for r in rows}
        if len(resp) > 1:
            raise ValueError(f"GSE78220 {pt} 应答标签不一致：{resp}")
        label = resp.pop()
        times = {char_val(r, "biopsy time") for r in rows}
        locs = [char_val(r, "anatomical location") for r in rows]

        dup = rule = notes = ""
        if len(rows) == 2:
            dup = "同患者治疗前两个不同解剖部位活检"
            rule = "患者内重复：group-aware 划分 + 仅取 Pt27A 的敏感性分析；不得平均或去重"
            notes = "论文原文：两灶 may not share the same transcriptomic profile"

        if times == {"on-treatment"}:
            inc_pri, inc_sen = "否", "是"
            why = "仅 on-treatment 样本，非治疗前基线，与开发队列时点不一致"
        else:
            inc_pri, inc_sen, why = "是", "是", ""

        out.append({
            "cohort": "GSE78220", "patient_id": pt, "role": "独立外部验证队列",
            "n_samples": len(rows),
            "gsm_list": join(r["gsm"] for r in rows),
            "sample_titles": join(titles),
            "visits": join(sorted(times)),
            "label_primary": label,
            "label_source": "Hugo 2016 Cell（PMID 26997480）RECIST best overall response；GEO anti-pd-1 response 字段",
            "key_covariates": join([
                f"site={char_val(rows[0], 'study site')}",
                f"biopsy_time={join(sorted(times))}",
                f"location={join(locs)}",
                f"OS_days={char_val(rows[0], 'overall survival (days)')}",
                f"vital={char_val(rows[0], 'vital status')}",
                f"braf={char_val(rows[0], 'braf') or 'WT/NA'}",
                f"lib={char_val(rows[0], 'stranded/unstranded rnaseq')}",
            ]),
            "dup_flag": dup, "dup_rule": rule,
            "include_primary": inc_pri, "include_sensitivity": inc_sen,
            "exclude_reason": why, "notes": notes,
        })
    return out


def build_gse244982(meta: Path) -> list[dict]:
    """耐药队列。41 患者，全部为 ICB 进展后取样，无应答者。
    核对结论：prior.ctla4i 语义为并集，标记不可用；11/41 处于 BRAFi 第 7 天暴露。
    """
    samples = read_tsv(meta / "GSE244982_samples.tsv")
    by_pt: dict[str, list[dict]] = defaultdict(list)
    for s in samples:
        by_pt[s["title"]].append(s)

    out = []
    for pt in sorted(by_pt, key=lambda x: int(x[3:])):
        r = by_pt[pt][0]
        groups = char_val(r, "groups")
        tx = char_val(r, "treatment")
        primary = char_val(r, "type.of.primary")
        brafi = char_val(r, "previous. brafi")

        inc_pri, inc_sen, why = "是", "是", ""
        if groups == "mucosal":
            inc_pri = "否"
            why = "原发黏膜黑色素瘤，原文已将其排除出统计"
        elif groups == "NA":
            inc_pri, inc_sen = "否", "是"
            why = "groups 字段漏填（Pat31）；该例为 CTLA4res + BRAFi day7"

        out.append({
            "cohort": "GSE244982", "patient_id": pt, "role": "耐药态刻画与方向一致性检验",
            "n_samples": len(by_pt[pt]),
            "gsm_list": join(x["gsm"] for x in by_pt[pt]),
            "sample_titles": pt,
            "visits": "post-progression",
            "label_primary": tx or groups,
            "label_source": "Nat Commun 2024（PMID 38594286）临床进展后取样；非 RECIST 判定，影像学标准不可得",
            "key_covariates": join([
                f"groups={groups}",
                f"primary={primary}",
                f"stage={char_val(r, 'stage')}",
                f"gender={char_val(r, 'gender')}",
                f"brafi={brafi}" + ("（BRAFi 联合暴露，强混杂）" if "During" in brafi else ""),
                "prior_ctla4i=不可用（并集语义，与 treatment 共线）",
            ]),
            "dup_flag": "", "dup_rule": "",
            "include_primary": inc_pri, "include_sensitivity": inc_sen,
            "exclude_reason": why,
            "notes": "本队列无应答样本，只能检验模型在进展样本上是否给出预期低分，不能验证区分能力",
        })
    return out


def build_gse294272(meta: Path) -> list[dict]:
    """29 样本 / 26 例独立患者。核对结论：HM003 Large/Small 同患者同天两灶；
    Mel-TIL-019 与 HM019、Mel-TIL-026 与 HM026 各为同一患者。
    """
    samples = read_tsv(meta / "GSE294272_samples.tsv")

    # 同一患者的不同写法归并
    alias = {
        "Mel-TIL-019": "HM019",
        "Mel-TIL-026": "HM026",
        "HM003-Large": "HM003",
        "HM003-Small": "HM003",
        "Mel-TIL-002": "MT002", "Mel-TIL-011": "MT011",
        "Mel-TIL-019_": "HM019", "Mel-TIL-022": "MT022",
        "Mel-TIL-026_": "HM026", "Mel-TIL-030": "MT030",
    }

    by_pt: dict[str, list[dict]] = defaultdict(list)
    for s in samples:
        base = s["title"].replace("_RNA-seq", "")
        pt = alias.get(base, base)
        by_pt[pt].append(s)

    out = []
    for pt in sorted(by_pt):
        rows = by_pt[pt]
        titles = [r["title"].replace("_RNA-seq", "") for r in rows]
        groups = [char_val(r, "treatment") for r in rows]

        dup = rule = notes = ""
        inc_pri, inc_sen, why = "是", "是", ""
        if pt == "HM003":
            dup = "同患者同一天的两个不同瘤灶"
            rule = "合并或择一；论文按样本统计存在伪重复"
            notes = "补充数据原文：different tumors collected on the same day from the same patient"
        elif pt in ("HM019", "HM026"):
            dup = "跨时点重复患者（Mel-TIL-### 与 HM### 同一人）"
            rule = "同一患者不得同时进入 Untreated 与 ICI 组对比；保留一个时点或整例剔除"
            notes = f"对应样本 {join(titles)}，分组 {join(groups)}"
            if pt == "HM026":
                notes += "；HM026 即 resistant/responder 泄漏"
            if len(set(groups)) > 1:
                inc_pri, inc_sen = "否", "是"
                why = "同一患者跨组重复，直接污染 Responder vs Resistant 对比"

        out.append({
            "cohort": "GSE294272", "patient_id": pt, "role": "去卷积参考谱的独立验证集",
            "n_samples": len(rows),
            "gsm_list": join(r["gsm"] for r in rows),
            "sample_titles": join(titles),
            "visits": join(groups),
            "label_primary": join(sorted(set(groups))),
            "label_source": "Nat Commun 2026（PMID 42277002）补充数据 1：treatment 分组",
            "key_covariates": "tissue=lymph node（含瘤整块，SPIAT 分 Tumour/Stroma）；batch 字段 GEO 未提供",
            "dup_flag": dup, "dup_rule": rule,
            "include_primary": inc_pri if dup else "是",
            "include_sensitivity": inc_sen if dup else "是",
            "exclude_reason": why if dup else "",
            "notes": notes,
        })
    return out


def build_gse308434(meta: Path) -> list[dict]:
    """42 样本 / 34 患者。核对结论：无应答标签；F(25 例) 与 R(9 例) 为两套独立编号。
    """
    samples = read_tsv(meta / "GSE308434_samples.tsv")

    def series(title: str) -> str:
        return "F" if re.match(r"^Melanoma_F\d", title) else "R"

    def base_pt(title: str) -> str:
        body = title.replace("Melanoma_", "")
        m = re.match(r"^(F\d+|R\d+)_", body)
        return m.group(1) if m else body

    by_pt: dict[str, list[dict]] = defaultdict(list)
    for s in samples:
        by_pt[base_pt(s["title"])].append(s)

    out = []
    for pt in sorted(by_pt, key=lambda x: (series(x + "_"), int(re.sub(r"\D", "", x)))):
        rows = by_pt[pt]
        titles = [r["title"] for r in rows]
        tx = sorted({char_val(r, "treatment") for r in rows})
        geno = sorted({char_val(r, "genotype") for r in rows})

        multi = len(rows) > 1
        dup = "同患者多活检" if multi else ""
        rule = "配对结构处理；不得当独立样本计数" if multi else ""
        notes = ""
        if "post1_pre2" in join(titles):
            notes = "post1_pre2 = 一线结束后、二线开始前的洗脱窗口活检（推断级，无作者字典）"
        if "naive" in join(tx):
            notes = (notes + "；" if notes else "") + "本队列唯一未治疗锚点"

        out.append({
            "cohort": "GSE308434", "patient_id": pt, "role": "模块 2 肿瘤内在状态×TCR 配对分析",
            "n_samples": len(rows),
            "gsm_list": join(r["gsm"] for r in rows),
            "sample_titles": join(titles),
            "visits": join(tx),
            "label_primary": "无（该系列无 response 字段）",
            "label_source": "SOFT 全字段确认无 response/RECIST/结局字段；Columbia Izar 组 GSE308433/434/435 三姊妹系列",
            "key_covariates": join([f"treatment={join(tx)}", f"genotype={join(geno)}",
                                    f"series={'F' if pt.startswith('F') else 'R'}"]),
            "dup_flag": dup, "dup_rule": rule,
            "include_primary": "是", "include_sensitivity": "是",
            "exclude_reason": "",
            "notes": (notes + "；" if notes else "") + "只做描述性/关联分析，不做应答监督预测",
        })
    return out


def build_gse115821() -> list[dict]:
    return [{
        "cohort": "GSE115821", "patient_id": "(队列级排除)", "role": "IMPRES 开发队列",
        "n_samples": 37, "gsm_list": "", "sample_titles": "", "visits": "PRE/On",
        "label_primary": "R 3 / NR 34",
        "label_source": "作者 GitHub MGH.mat 逐样本字段（Supp. Table 9 无法获取）",
        "key_covariates": "去重后仅 8 例独立患者；3 个 R 样本中 2 个属同一患者 pre/on 配对",
        "dup_flag": "", "dup_rule": "",
        "include_primary": "否", "include_sensitivity": "否",
        "exclude_reason": "独立应答患者约 2 例，统计上不足以建模或做外部验证",
        "notes": "详见 docs/02_回原文核对结论.md 第七之二节",
    }]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    meta = Path(args.outdir) / "meta"

    rows = []
    rows += build_gse91061(meta)
    rows += build_gse78220(meta)
    rows += build_gse244982(meta)
    rows += build_gse294272(meta)
    rows += build_gse308434(meta)
    rows += build_gse115821()

    out_tsv = meta / "patient_level.tsv"
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"患者级纳排表已写入 {out_tsv}（{len(rows)} 行）")

    # 汇总
    lines = ["# 患者级纳排统计", "",
             "由 `scripts/05_build_patient_level.py` 自动生成，与 `meta/patient_level.tsv` 同源。", "",
             "| 队列 | 患者数 | 样本数 | 主分析纳入 | 敏感性分析纳入 | 排除 | 角色 |",
             "|---|---:|---:|---:|---:|---:|---|"]
    for coh in ["GSE91061", "GSE78220", "GSE244982", "GSE294272", "GSE308434", "GSE115821"]:
        rs = [r for r in rows if r["cohort"] == coh]
        if not rs:
            continue
        ns = sum(int(r["n_samples"]) for r in rs)
        inc = sum(1 for r in rs if r["include_primary"] == "是")
        sen = sum(1 for r in rs if r["include_sensitivity"] == "是")
        exc = len(rs) - inc
        lines.append(f"| {coh} | {len(rs)} | {ns} | {inc} | {sen} | {exc} | {rs[0]['role']} |")

    lines += ["", "## 排除明细", ""]
    for coh in ["GSE91061", "GSE78220", "GSE244982", "GSE294272", "GSE308434", "GSE115821"]:
        exc = [r for r in rows if r["cohort"] == coh and r["include_primary"] != "是"]
        if not exc:
            continue
        lines.append(f"### {coh}（{len(exc)} 例）")
        for r in exc:
            lines.append(f"- `{r['patient_id']}` — {r['exclude_reason'] or '（见 notes）'}")
        lines.append("")

    lines += ["## 重复/冲突标记汇总", ""]
    for coh in ["GSE91061", "GSE78220", "GSE294272", "GSE308434"]:
        dups = [r for r in rows if r["cohort"] == coh and r["dup_flag"]]
        lines.append(f"- **{coh}**：{len(dups)} 例带重复标记")
        for r in dups:
            lines.append(f"  - `{r['patient_id']}`（{r['n_samples']} 样本）— {r['dup_flag']}；处理：{r['dup_rule']}")
    lines.append("")

    (meta / "patient_level_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"汇总已写入 {meta / 'patient_level_summary.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())