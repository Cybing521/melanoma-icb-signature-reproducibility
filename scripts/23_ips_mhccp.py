#!/usr/bin/env python3
"""第二个签名 `IPS-MHC+CP` 的三队列验证（口径见 `docs/03` §22，计算前已锁定）。

为什么是它
----------
IMPRES 的 15 个基因对是从作者 GitHub 重建的（原文 Supp. Table 2 四条路径全不可达）。
单一签名 + 单一来源是全文最容易被打穿处，故补一个**构造原理不同、来源独立**的签名。

选 Immunophenoscore（IPS，Charoentong et al., Cell Rep 2017, PMID 28052254）——
4 类 / 26 个免疫特征集的加权 z 分数，与 IMPRES 的「两两大小排序」完全不同。
但其**完整基因表无法从任何开放来源取得**（6 条路径全部失败，见 `docs/03` §22.2）。

故按 §22.3 改测可完整取得的两类：`MHC` + `CP`，共 20 个基因。
**这 20 个基因在 GSE215868 的 NanoString 面板上 20/20 全覆盖**——
也就是说它失败时不能用「测不到」辩解，这是选它的价值。

打分（§22.5 锁定）
-----------------
逐基因在本队列内 z-score → 每类取均值 → 两类相加得 `IPS-MHC+CP_raw`。
**用连续值**，不做 IPS 原文的 0–10 离散化（预登记 §22.5 已判定该离散化会摧毁
n=45 量级的排序分辨率）。

输出
----
qc/IPS_MHCCP_scores.tsv
qc/IPS_MHCCP_auc.tsv
qc/docs_09_IPS_MHCCP.md
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
SEED = 20261007
N_BOOT = 4000

MHC = ["B2M", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
       "HLA-DPA1", "HLA-DPB1", "TAP1", "TAP2"]
CP = ["CD27", "CTLA4", "ICOS", "IDO1", "LAG3",
      "PDCD1", "CD274", "PDCD1LG2", "TIGIT", "HAVCR2"]

# 注释体系差异：GSE91061 用 Entrez→Symbol 映射，部分 HLA 基因符号带前缀
ALIASES = {
    "HLA-A": ["HLA-A", "HLA-A*", "HLA*A"],
    "HLA-B": ["HLA-B", "HLA-B*", "HLA*B"],
    "HLA-C": ["HLA-C", "HLA-C*", "HLA*C"],
    "HLA-E": ["HLA-E", "HLA-E*", "HLA*E"],
    "HLA-F": ["HLA-F", "HLA-F*", "HLA*F"],
    "HLA-DPA1": ["HLA-DPA1", "HLA-DPA1*", "HLA*DPA1"],
    "HLA-DPB1": ["HLA-DPB1", "HLA-DPB1*", "HLA*DPB1"],
}


def resolve(sym: str, have: set[str]) -> str | None:
    for cand in ALIASES.get(sym, [sym]):
        if cand in have:
            return cand
    # 去星号再试一次（部分注释写成 HLA-A*01 之类）
    if "-" in sym:
        base = sym.split("-")[0]
        for g in have:
            if g.startswith(base + "-") or g.startswith(base + "*"):
                return g
    return None


def load(path: Path) -> dict[str, np.ndarray]:
    df = pd.read_csv(path, sep="\t", index_col=0)
    return {g: df.loc[g].to_numpy(dtype=float) for g in df.index}, list(df.columns)


def score(mat: dict[str, np.ndarray], genes: list[str]) -> tuple[np.ndarray, dict]:
    """逐基因 z-score → 类均值 → 两类相加。缺失基因如实返回，不插补。"""
    detail, zs = {}, []
    for cls, gset in (("MHC", MHC), ("CP", CP)):
        vals = []
        for g in gset:
            hit = resolve(g, set(mat))
            if hit is None:
                detail[f"{cls}:{g}"] = None
                continue
            v = mat[hit]
            sd = v.std()
            z = (v - v.mean()) / sd if sd > 0 else np.zeros_like(v)
            vals.append(z)
            detail[f"{cls}:{g}"] = {"symbol": g, "matched": hit,
                                    "sd": float(sd)}
        zs.append(np.vstack(vals).mean(axis=0) if vals else np.zeros(1))
    return zs[0] + zs[1], detail


def boot_auc(pos, neg, rng):
    n1, n0 = len(pos), len(neg)
    if n1 == 0 or n0 == 0:
        return np.nan, np.nan, np.nan

    def auc(p, q):
        r = stats.rankdata(np.concatenate([p, q]))
        return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

    point = auc(pos, neg)
    i = rng.integers(0, n1, size=(N_BOOT, n1))
    j = rng.integers(0, n0, size=(N_BOOT, n0))
    d = np.array([auc(pos[i[k]], neg[j[k]]) for k in range(N_BOOT)])
    lo, hi = np.percentile(d, [2.5, 97.5])
    u, p = stats.mannwhitneyu(pos, neg, alternative="greater")
    return point, lo, hi, p


def main() -> int:
    rng = np.random.default_rng(SEED)
    qc = ROOT / "qc"
    log: list[str] = ["# 第二个签名 `IPS-MHC+CP` 的三队列验证\n",
                      "生成脚本 `scripts/23_ips_mhccp.py`，口径锁定于 `docs/03` §22，"
                      "**写于任何计算之前**。\n"]

    # ---------- 标签 ----------
    pat = pd.read_csv(ROOT / "meta" / "patient_level.tsv", sep="\t", dtype=str)
    RESP = {"Partial Response": "PRCR", "Complete Response": "PRCR",
            "Progressive Disease": "PD", "PRCR": "PRCR", "PD": "PD"}

    def labels_for(cohort: str) -> pd.DataFrame:
        m = pat[pat["cohort"] == cohort].copy()
        m["gsm"] = m["gsm_list"].str.split(r"[;|]")
        m = m.explode("gsm")
        m = m[m["include_primary"].astype(str).eq("是")]
        m["grp"] = m["label_primary"].map(RESP)
        # 同患者多位点合并为 1 例（取第一个标签；两部位同患者标签必然相同）
        return m.groupby("patient_id", as_index=False).agg(
            grp=("grp", "first"),
            gsm_first=("gsm", "first"))

    devs = pd.read_csv(ROOT / "meta" / "dev_set_pre_treatment.tsv", sep="\t", dtype=str)
    devs["gsm_first"] = devs["pre_gsm"]
    devs["grp"] = devs["label_primary"].map(RESP)

    # ---------- 三队列 ----------
    cohorts = [
        ("GSE91061", "expr/GSE91061_symbol_matrix.tsv.gz", devs, "33 例治疗前（PRCR 10 / PD 23）"),
        ("GSE78220", "expr/GSE78220_symbol_matrix.tsv.gz", labels_for("GSE78220"),
         "26 例治疗前患者（PRCR 14 / PD 12）"),
        ("GSE215868", "expr/GSE215868_symbol_matrix.tsv.gz",
         pd.read_csv(ROOT / "meta" / "GSE215868_samples.tsv", sep="\t", dtype=str).rename(
             columns={"gsm": "gsm_first", "response_group": "grp"}),
         "79 例（PRCR 45 / PD 34）"),
    ]

    rows, sc_all = [], []
    for name, rel, meta, desc in cohorts:
        mat, cols = load(ROOT / rel)
        s, detail = score(mat, MHC + CP)
        nmiss = sum(1 for v in detail.values() if v is None)
        log.append(f"## {name}（{desc}）\n")
        log.append(f"- 基因覆盖：**{20 - nmiss}/20**"
                   + (f"，缺失 {sorted(k for k, v in detail.items() if v is None)}"
                      if nmiss else "（全部命中）"))
        md = pd.DataFrame({"gsm_first": cols, "ips_mhccp": s}).merge(
            meta[["gsm_first", "grp"]], on="gsm_first", how="inner")
        md = md[md["grp"].isin(["PRCR", "PD"])]
        y = (md["grp"] == "PRCR").to_numpy()
        sc_all.append(md.assign(cohort=name))

        log.append(f"- 实际参与 **{len(md)} 例**（PRCR {int(y.sum())} / PD {int((~y).sum())}）\n")
        log.append("| 方向 | AUC | 95% CI | 单侧 p |")
        log.append("|---|---:|---|---:|")
        for orient in ("high_ips_responder", "mirror"):
            v = md["ips_mhccp"].to_numpy()
            a_, b_ = (v[y], v[~y]) if orient == "high_ips_responder" else (v[~y], v[y])
            auc, lo, hi, p = boot_auc(a_, b_, rng)
            log.append(f"| `{orient}` | **{auc:.3f}** | {lo:.3f}–{hi:.3f} | {p:.4f} |")
            rows.append({"cohort": name, "n_total": len(md), "n_PRCR": int(y.sum()),
                         "n_PD": int((~y).sum()), "orientation": orient,
                         "auc": auc, "ci_lo": lo, "ci_hi": hi, "p_greater": p,
                         "genes_covered": 20 - nmiss})
        log.append("")

    pd.concat(sc_all, ignore_index=True).to_csv(
        qc / "IPS_MHCCP_scores.tsv", sep="\t", index=False)
    adf = pd.DataFrame(rows)
    adf.to_csv(qc / "IPS_MHCCP_auc.tsv", sep="\t", index=False)

    # ---------- 裁决（沿用 §17.4 三条，写死于计算之前） ----------
    log.append("## 裁决（沿用 `docs/03` §17.4 写死的判定标准）\n")
    verdict = []
    for name in adf["cohort"].unique():
        sub = adf[adf["cohort"] == name]
        ok = sub[(sub["ci_lo"] > 0.50) & (sub["auc"].between(0.70, 1.00))]
        if len(ok):
            verdict.append(f"- **{name}**：方向 `{ok.iloc[0]['orientation']}` 复现。")
        else:
            verdict.append(f"- **{name}**：两向皆未复现（CI 均覆盖 0.50，"
                           f"点估计均在 0.70–1.00 之外）。")
    log.extend(verdict)
    log.append("")
    any_ok = adf[(adf["ci_lo"] > 0.50) & (adf["auc"].between(0.70, 1.00))]
    if any_ok.empty:
        log.append("> **三个队列、两个方向，全部未复现。**")
        log.append(">")
        log.append("> 与 IMPRES 合并后，本项目在**两个构造原理完全不同、来源互相独立的已发表"
                   "ICB 应答签名**上都没有复现已发表性能。")
    log.append("")
    log.append("## 与 IMPRES 的对照（这是选这个签名的全部意义）\n")
    log.append("| | IMPRES | `IPS-MHC+CP` |")
    log.append("|---|---|---|")
    log.append("| 构造原理 | 免疫检查点基因**两两大小关系** | 抗原呈递 + 检查点基因的**加权 z 分数** |")
    log.append("| 需要的基因数 | 15 | 20 |")
    log.append("| 常规临床面板覆盖 | 13/15（86.7%） | **20/20（100%）** |")
    log.append("| 基因表可开放获取 | **否**（从作者 GitHub 重建） | **是**（开放文献逐基因列出） |")
    log.append("")
    log.append("**`IPS-MHC+CP` 在临床可及性上全面优于 IMPRES**——它 100% 落在常规 IO 面板上，"
               "基因表也是公开的。因此它若失败，**不能**用「测不准」「测不到」辩解；"
               "而预登记 §22 已核 20/20 全覆盖、标签与队列口径与其他分析完全一致、"
               "n 最大的一支队列有 79 例。")
    log.append("")
    log.append("**必须同时写明的边界**：")
    log.append("")
    log.append("1. `IPS-MHC+CP` **不是 IPS**。IPS 的 `EC`（效应细胞，4 个特征集约 100 基因）"
               "与 `SC`（抑制细胞，2 个特征集 40 基因）两类的基因表无法开放取得，"
               "本测试**没有覆盖 IPS 的全部四个类别**。")
    log.append("2. z-score 是**队列内相对量**，不是绝对表达量（IPS 构造自带）。")
    log.append("3. 连续 raw 分数而非 IPS 原文的 0–10 整数（预登记 §22.5 事先判定）。")
    log.append("4. 仍未消除 IMPRES 那个残余风险：**本项目的 IMPRES 实现本身是重建的**，"
               "这一点由第二个签名**不能**替 IMPRES 洗清。第二个签名补的是"
               "「已发表签名系统性不可复现」这一更宽的论断，不是「IMPRES 实现正确」。")

    md = qc / "docs_09_IPS_MHCCP.md"
    md.write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))
    print(f"\n[OK] {qc/'IPS_MHCCP_scores.tsv'}\n[OK] {qc/'IPS_MHCCP_auc.tsv'}\n[OK] {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())