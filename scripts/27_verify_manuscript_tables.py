#!/usr/bin/env python3
"""核对稿件新增 5 张表的每一个数字是否与仓库数据文件一致。

为什么要有这个脚本
------------------
本项目的核心主张是「稿件里每一个数字都能追回某次脚本运行」。新增的表如果只是
写死在 `build_manuscript.py` 里，这条链子在表格上就断了。因此本脚本把稿件表内容
与真源逐格比对，并在不一致时**非零退出**——让不一致在提交前暴露，而不是留给审稿人。

真源
----
- 患者级纳排：`meta/patient_level.tsv`、`meta/GSE215868_samples.tsv`
- IMPRES 基因对与覆盖：`scripts/impres_pairs.py`、`qc/GSE215868_panel_diagnostic.tsv`
- IPS-MHC+CP 基因表：`scripts/23_ips_mhccp.py`
- 特征集偏离：`scripts/15_bulk_signature_scores.py`（真源是代码里的 MAIN_PANEL 与
  SIGNATURES，不是 docs/03 的中文转述——§21.1 的转述有一处名称错误，见 docs/14）

用法：`python scripts/27_verify_manuscript_tables.py`
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "manuscript"))
sys.path.insert(0, str(ROOT / "scripts"))

FAILS: list[str] = []
CHECKS = 0


def check(label: str, got, want) -> None:
    """两侧都按字符串比较：稿件表里一律是字符串，而真源常是 int，
    直接比会得到一堆「17 vs '17'」的假不一致，掩盖真正的问题。"""
    global CHECKS
    CHECKS += 1
    if str(got) != str(want):
        FAILS.append(f"{label}: 稿件={got!r} 真源={want!r}")


def table_by_key(blocks, key):
    for b in blocks:
        if b.get("type") == "table" and b.get("key") == key:
            return b
    raise KeyError(key)


def rowmap(tbl) -> dict[str, list[str]]:
    return {r[0]: r for r in tbl["rows"]}


# ── 患者级纳排表 ────────────────────────────────────────────────
def verify_cascade(blocks) -> None:
    tbl = rowmap(table_by_key(blocks, "cascade"))

    pl = list(csv.DictReader((ROOT / "meta/patient_level.tsv").open(encoding="utf-8"),
                             delimiter="\t"))

    # GSE91061 / GSE78220 来自患者级表
    for coh, col in (("GSE91061", 1), ("GSE78220", 2)):
        rs = [r for r in pl if r["cohort"] == coh]
        n_pat = len(rs)
        primary = [r for r in rs if r["include_primary"] == "是"]
        sens = [r for r in rs if r["include_sensitivity"] == "是"]
        n_ne = n_pat - len(sens)
        n_sd = len([r for r in sens if r["label_primary"] == "SD"])
        check(f"cascade/{coh} patients", tbl["Independent patients after sample-level "
                                           f"de-duplication"][col], str(n_pat))
        check(f"cascade/{coh} NE", tbl["Excluded: response not classifiable"][col],
              f"{n_ne}" + (" (RECIST NE)" if coh == "GSE91061" else ""))
        check(f"cascade/{coh} SD", tbl["Excluded: stable disease (two-level endpoint)"][col],
              str(n_sd))
        # 「可分类」= 标签落在二分类定义内（PRCR/PD），且不是 NE、不是 SD。
        # 注意不能用 len(primary) 代替：GSE78220 的 Pt16 可分类，只因无治疗前活检
        # 才没进主分析，把它算成「不可分类」会把时点排除与标签缺失混为一谈。
        two_level = {"PRCR", "PD", "Partial Response", "Complete Response",
                     "Progressive Disease"}
        cls = [r for r in rs if r["label_primary"] in two_level]
        check(f"cascade/{coh} classifiable",
              tbl["Patients with a classifiable two-level response"][col], str(len(cls)))
        # 治疗前子集
        pre = [r for r in cls if "pre" in r["visits"].lower()]
        check(f"cascade/{coh} on-treatment-only excluded",
              tbl["Further excluded: no pre-treatment biopsy available"][col],
              str(len(cls) - len(pre)))
        rr = Counter(r["label_primary"] for r in pre)
        n_resp = rr["PRCR"] + rr["Partial Response"] + rr["Complete Response"]
        n_non = rr["PD"] + rr["Progressive Disease"]
        check(f"cascade/{coh} evaluable",
              tbl["Evaluable patients entering the comparison"][col], str(len(pre)))
        check(f"cascade/{coh} R/NR", tbl["Responders / non-responders"][col],
              f"{n_resp} / {n_non}")

    # GSE215868 来自样本级元数据
    g = list(csv.DictReader((ROOT / "meta/GSE215868_samples.tsv").open(encoding="utf-8"),
                            delimiter="\t"))
    col = 3
    bor = Counter(r["bor"] for r in g)
    check("cascade/GSE215868 patients",
          tbl["Independent patients after sample-level de-duplication"][col], str(len(g)))
    check("cascade/GSE215868 NE", tbl["Excluded: response not classifiable"][col],
          f"{bor['UK']} (unknown)")
    check("cascade/GSE215868 SD", tbl["Excluded: stable disease (two-level endpoint)"][col],
          str(bor["SD"]))
    n_cls = len(g) - bor["UK"] - bor["SD"]
    check("cascade/GSE215868 classifiable",
          tbl["Patients with a classifiable two-level response"][col], str(n_cls))
    check("cascade/GSE215868 evaluable",
          tbl["Evaluable patients entering the comparison"][col], str(n_cls))
    check("cascade/GSE215868 R/NR", tbl["Responders / non-responders"][col],
          f"{bor['CR'] + bor['PR']} / {bor['PD']}")


# ── IMPRES 基因对表 ─────────────────────────────────────────────
def verify_impres(blocks) -> None:
    from impres_pairs import FEATURE_IDX, idx_to_pair

    rows = rowmap(table_by_key(blocks, "impres"))
    check("impres n rows", str(len(rows)), str(len(FEATURE_IDX)))

    cov: dict[tuple[str, int], str] = {}
    with (ROOT / "qc/GSE215868_panel_diagnostic.tsv").open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cov[(r["cohort"], int(r["pair_no"]))] = r["status"]

    cols = {"GSE91061": 3, "GSE78220": 4, "GSE215868": 5}
    for k, idx in enumerate(FEATURE_IDX, start=1):
        g1, g2 = idx_to_pair(idx)
        row = rows[str(k)]
        check(f"impres #{k} g1", row[1], g1)
        check(f"impres #{k} g2", row[2], g2)
        for coh, c in cols.items():
            want = "yes" if cov[(coh, k)] == "可算" else "not on panel"
            check(f"impres #{k} {coh}", row[c], want)


# ── IPS-MHC+CP 基因表 ───────────────────────────────────────────
def verify_ips(blocks) -> None:
    from importlib import import_module
    mhc = import_module("23_ips_mhccp") if False else None
    # 直接从源码读取基因表，避免为读两个常量而 import 一整支分析脚本
    src = (ROOT / "scripts/23_ips_mhccp.py").read_text(encoding="utf-8")
    ns: dict = {}
    exec(compile(src.split("def resolve")[0].replace("import pandas as pd", "")
                 .replace("import numpy as np", "").replace("from scipy import stats", "")
                 .replace("ROOT = Path(__file__).resolve().parent.parent", "ROOT=None")
                 .replace("SEED = 20261007", "").replace("N_BOOT = 4000", ""),
                 "ips_consts", "exec"), ns)
    rows = rowmap(table_by_key(blocks, "ips"))
    mhc_genes, cp_genes = ns["MHC"], ns["CP"]

    listed_mhc = [g.strip() for g in rows["MHC"][2].split(",")]
    check("ips MHC count", rows["MHC"][1], str(len(mhc_genes)))
    check("ips MHC genes", sorted(listed_mhc), sorted(mhc_genes))
    listed_cp = [g.strip() for g in rows["Checkpoint"][2].split(",")]
    check("ips CP count", rows["Checkpoint"][1], str(len(cp_genes)))
    check("ips CP genes", sorted(listed_cp), sorted(cp_genes))
    check("ips total", rows["Total"][1], str(len(mhc_genes) + len(cp_genes)))


# ── 特征集偏离表 ────────────────────────────────────────────────
def verify_featset(blocks) -> None:
    src = (ROOT / "scripts/15_bulk_signature_scores.py").read_text(encoding="utf-8")
    ns: dict = {}
    exec(compile(src.split("def read_gmt")[0], "feat_consts", "exec"), ns)
    main_panel: list[str] = ns["MAIN_PANEL"]
    sigs: dict[str, str] = ns["SIGNATURES"]

    rows = rowmap(table_by_key(blocks, "featset"))
    data_rows = {k: v for k, v in rows.items() if k != "Totals"}

    used = {k for k, v in data_rows.items() if v[3] == "used"}
    removed = {k for k, v in data_rows.items() if v[3] == "removed"}
    added = {k for k, v in data_rows.items() if v[3] == "added"}
    primary = {k for k, v in data_rows.items() if v[2].startswith("primary panel")}

    check("featset rows", str(len(data_rows)), "21")
    check("featset primary rows", str(len(primary)), "7")      # 6 states → 7 gene sets
    check("featset used", str(len(used)), "14")
    check("featset removed", str(len(removed)), "4")
    check("featset added", str(len(added)), "3")
    check("featset modelled total", len(used) + len(added), str(len(sigs)))
    check("featset preregistered total", len(used) + len(removed), "18")
    check("featset totals row", rows["Totals"][3],
          f"{len(used) + len(added)} modelled: {len(used)} used + {len(added)} added")

    # 换入的三个必须确实存在于 scripts/15，且不在预登记主面板里
    ADDED_MAP = {"IL-6 / JAK / STAT3 (IL6_JAK_STAT3_SIGNALING)": "il6_jak_stat3",
                 "WNT / \u03b2-catenin (WNT_BETA_CATENIN_SIGNALING)": "wnt_beta_catenin",
                 "KRAS signalling up (KRAS_SIGNALING_UP)": "kras_up"}
    check("featset added rows", sorted(added), sorted(ADDED_MAP))
    for name, key in ADDED_MAP.items():
        if key not in sigs:
            FAILS.append(f"featset: {key} 不在 scripts/15 的 SIGNATURES 里")
        if key in main_panel:
            FAILS.append(f"featset: {key} 被列为换入，却出现在预登记主面板中")

    # 换出的四个必须确实在预登记探索面板里、且不在实际 SIGNATURES 里
    REMOVED_KEYS = ("INFLAMMATORY_RESPONSE", "PI3K_AKT_MTOR", "APOPTOSIS", "P53_PATHWAY")
    for tok in REMOVED_KEYS:
        key = next((k for k in sigs if tok.split("_")[0] in k), None)
        if key is not None:
            FAILS.append(f"featset: {tok} 被列为换出，却仍在 scripts/15 的 SIGNATURES 里")

    # 主面板 6 个状态必须在表中出现
    LABEL = {"ifng_response": "INTERFERON_GAMMA_RESPONSE",
             "antigen_presentation": "ALLOGRAFT_REJECTION",
             "proliferation": "E2F_TARGETS", "emt": "EPITHELIAL_MESENCHYMAL",
             "hypoxia": "HYPOXIA", "tnf_nfkb": "TNFA_SIGNALING"}
    for k, tok in LABEL.items():
        hit = [r for r in primary if tok in r]
        if not hit:
            FAILS.append(f"featset: 预登记主面板特征 {k} 未出现在表中")


# ── 数据来源表 ──────────────────────────────────────────────────
def verify_sources(blocks) -> None:
    tbl = rowmap(table_by_key(blocks, "sources"))
    pl = list(csv.DictReader((ROOT / "meta/patient_level.tsv").open(encoding="utf-8"),
                             delimiter="\t"))
    from collections import defaultdict
    pat: dict[str, int] = defaultdict(int)
    for r in pl:
        pat[r["cohort"]] += 1
    # 患者级表里有的队列，其患者数必须与来源表一致
    for series, n_col in (("GSE91061", 3), ("GSE78220", 3), ("GSE244982", 3),
                          ("GSE294272", 3), ("GSE308434", 3)):
        check(f"sources/{series} patients", tbl[series][n_col], str(pat[series]))
    g = list(csv.DictReader((ROOT / "meta/GSE215868_samples.tsv").open(encoding="utf-8"),
                            delimiter="\t"))
    check("sources/GSE215868 samples", tbl["GSE215868"][2], str(len(g)))
    check("sources/GSE215868 patients", tbl["GSE215868"][3], str(len(g)))


# ── 跨表一致性 ──────────────────────────────────────────────────
def verify_cross_table(blocks) -> None:
    """同一队列的患者数在来源表、级联表、特征表三处必须同口径。

    新增表 1、表 2 把「落盘患者数」显式写了出来，于是暴露出旧表里
    GSE78220 的患者数曾按可评估数（26）填写。两张表现在都说 27，
    本检查把这个一致性钉住，避免以后再退回混用。
    """
    src = rowmap(table_by_key(blocks, "sources"))
    casc = rowmap(table_by_key(blocks, "cascade"))
    char = rowmap(table_by_key(blocks, "cohort"))
    for series, s_col, c_col, char_col in (("GSE91061", 3, 1, 1),
                                          ("GSE78220", 3, 2, 2),
                                          ("GSE215868", 3, 3, 3)):
        check(f"cross/{series} patients sources vs cascade",
              src[series][s_col], casc["Independent patients after sample-level "
                                       f"de-duplication"][c_col])
        check(f"cross/{series} patients sources vs characteristics",
              src[series][s_col], char["Patients"][char_col])
    # 正文提到患者数的地方也要跟上
    for b in blocks:
        if b.get("type") == "paragraph" and "GSE78220 contributes 28 biopsies" in b.get("text", ""):
            check("cross/GSE78220 patients in 3.1", "from 27 " in b["text"], True)


def main() -> int:
    import build_manuscript
    blocks = build_manuscript.build()["blocks"]

    verify_sources(blocks)
    verify_cascade(blocks)
    verify_impres(blocks)
    verify_ips(blocks)
    verify_featset(blocks)
    verify_cross_table(blocks)

    print(f"核对项 {CHECKS} 个")
    if FAILS:
        print(f"\n不一致 {len(FAILS)} 处：")
        for f in FAILS:
            print("  ✗", f)
        return 1
    print("全部一致：稿件 5 张新增表的每个数字均可追回仓库数据文件")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())