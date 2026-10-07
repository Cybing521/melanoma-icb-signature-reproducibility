#!/usr/bin/env python3
"""在四个 bulk 队列上计算候选特征——模块 4 建模的输入。

特征集（17 个 MSigDB Hallmark 签名，全部可从常规 bulk 转录组算得）
----------------------------------------------------------------
主面板 6 个（`docs/03` 第十一节已锁定）：
    IFN-γ 应答 / 抗原呈递 / 增殖 / EMT / 缺氧 / TNFα-NF-κB
探索面板 11 个（只作探索性，不用于支持主要结论）

打分方法
--------
**样本内基因 z-score 后取签名内均值**（ssGSEA 之外最常见的一路）。
必须用样本内 z-score 而非样本间标准化：我们要的是"该签名在这个样本内部
有多突出"，而非跨样本的表达高低。

⚠️ 恒零量：样本内 z-score 后，签名的**样本均值在数学上恒等于 0**，没有信息。
因此只输出 `*_sd`（样本内状态分离散度）、`*_p90_p10`、`*_hi_pct`，
不输出 `*_mean`。

输入：expr/<cohort>_symbol_matrix.tsv.gz
输出：qc/signature_scores_<cohort>.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path

SIGNATURES: dict[str, str] = {
    # ---- 主面板 6 个（docs/03 第十一节已锁定）----
    "ifng_response": "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "antigen_presentation": "HALLMARK_ALLOGRAFT_REJECTION",   # MSigDB v7 起原集下架
    "proliferation": "HALLMARK_E2F_TARGETS",
    "emt": "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "hypoxia": "HALLMARK_HYPOXIA",
    "tnf_nfkb": "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    # ---- 探索面板 11 个 ----
    # 选取依据是**已发表的黑色素瘤 ICB 耐药机制**，不是在本队列上试出来的。
    # 集名按 MSigDB Hallmark 2024.1.0 实际清单核对过，该版本里
    # HALLMARK_NGF 与 HALLMARK_ESTROGEN_RESPONSE 都不存在
    # （后者已拆为 _EARLY / _LATE），不要凭印象写名字。
    "interferon_alpha": "HALLMARK_INTERFERON_ALPHA_RESPONSE",   # IFN-α 与 ICB 应答相关
    "il2_stat5": "HALLMARK_IL2_STAT5_SIGNALING",                # T/NK 活化
    "il6_jak_stat3": "HALLMARK_IL6_JAK_STAT3_SIGNALING",         # 炎症性旁路，常与 IFN-γ 竞争
    "tgf_beta": "HALLMARK_TGF_BETA_SIGNALING",                  # 经典 ICB 耐药介质
    "wnt_beta_catenin": "HALLMARK_WNT_BETA_CATENIN_SIGNALING",   # β-catenin 活化是已知的免疫逃逸机制
    "kras_up": "HALLMARK_KRAS_SIGNALING_UP",                    # 黑色素瘤主导驱动
    "angiogenesis": "HALLMARK_ANGIOGENESIS",
    "oxidative_phosphorylation": "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "glycolysis": "HALLMARK_GLYCOLYSIS",
    "mtorc1": "HALLMARK_MTORC1_SIGNALING",
    "g2m_checkpoint": "HALLMARK_G2M_CHECKPOINT",                # 与 E2F 增殖互补
}
MAIN_PANEL = ["ifng_response", "antigen_presentation", "proliferation", "emt", "hypoxia", "tnf_nfkb"]


def read_gmt(p: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    with p.open(encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                out[parts[0]] = {g.strip().upper() for g in parts[2:] if g.strip()}
    return out


def read_matrix(path: Path) -> tuple[list[str], list[str], dict[str, list[float]]]:
    op = gzip.open if path.suffix == ".gz" else open
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        cols = next(csv.reader(fh, delimiter="\t"))
        samples = cols[1:]
        genes, expr = [], {}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[0].upper())
            expr[parts[0].upper()] = [float(v) for v in parts[1:]]
    return samples, genes, expr


def zscore_by_sample(values: list[float]) -> list[float]:
    n = len(values)
    mu = sum(values) / n
    var = sum((v - mu) ** 2 for v in values) / (n - 1) if n > 1 else 0.0
    sd = math.sqrt(var)
    if sd == 0:
        return [0.0] * n
    return [(v - mu) / sd for v in values]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--gmt", default="resources/hallmark_2024.1.Hs.gmt")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    gmt = read_gmt(outdir / args.gmt)
    mtx = outdir / "expr" / f"{args.cohort}_symbol_matrix.tsv.gz"
    samples, genes, expr = read_matrix(mtx)

    # 全基因的样本内 z-score 只算一次
    print(f"{args.cohort}: {len(samples)} 样本 × {len(genes)} 基因；先算全基因样本内 z-score…")
    Z: dict[str, list[float]] = {}
    for g in genes:
        Z[g] = zscore_by_sample(expr[g])

    rows = []
    for idx, s in enumerate(samples):
        row: dict[str, object] = {"sample": s}
        for short, full in SIGNATURES.items():
            gs = sorted(g for g in gmt[full] if g in Z)
            if not gs:
                row[f"{short}_genes"] = 0
                continue
            vals = [Z[g][idx] for g in gs]
            n = len(vals)
            m = sum(vals) / n
            sd = math.sqrt(sum((v - m) ** 2 for v in vals) / (n - 1)) if n > 1 else 0.0
            sv = sorted(vals)
            row[f"{short}_sd"] = round(sd, 5)
            row[f"{short}_p90_p10"] = round(sv[int(0.9 * (n - 1))] - sv[int(0.1 * (n - 1))], 5)
            row[f"{short}_hi_pct"] = round(sum(1 for v in vals if v > 0) / n, 5)
            row[f"{short}_genes"] = n
        rows.append(row)

    cols = ["sample"]
    for short in SIGNATURES:
        cols += [f"{short}_sd", f"{short}_p90_p10", f"{short}_hi_pct", f"{short}_genes"]
    out = outdir / "qc" / f"signature_scores_{args.cohort}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"特征 {len(SIGNATURES)} 个（主面板 {len(MAIN_PANEL)} + 探索 {len(SIGNATURES) - len(MAIN_PANEL)}）")
    for short in SIGNATURES:
        k = f"{short}_genes"
        n = rows[0].get(k)
        print(f"  {short:<26} 命中基因 {n}")
    print(f"已写入 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
