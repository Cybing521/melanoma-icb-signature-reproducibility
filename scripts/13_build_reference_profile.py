#!/usr/bin/env python3
"""模块 3 前置：从 GSE294273 的 12 例淋巴结单细胞数据构建去卷积参考谱。

为什么需要这一步
----------------
方案模块 3 要求把"细胞比例型特征"映射到常规 bulk 转录组。去卷积需要一张
**细胞类型 × 基因**的参考谱。GSE294273 是方案指定的参考谱来源（12 例供体，
未治疗 4 / ICI 耐药 4 / ICI 应答 4，全部为含瘤整块淋巴结）。

数据事实（已核）
----------------
* GEO deposit 为标准 10x 三件套：12 个样本 × (barcodes/features/matrix)，共 36 个文件
* CellRanger 6.0.2，GRCh38，features 文件同时给 Ensembl ID 与 symbol
* **GEO 未提供任何细胞类型注释**，元数据只有 tissue 与 treatment 两个字段

因此细胞类型必须从头标注。本脚本的做法是：聚类后用**公开的经典标志基因集**
对每个簇打分，按最高分指派细胞类型，并把每簇 top markers 一并导出供人工复核。
标志基因集写死在脚本里，便于第三方复核与替换。

质控阈值（依本队列自身分布，出结果前设定）
------------------------------------------
不套用其他队列的经验阈值；先看单细胞常规下限的实际剔除量，阈值记在输出里。
本脚本先按 `n_genes ≥ 200 且 n_umi ≥ 500 且 mt_pct ≤ 20` 做一次，
该组合对 scRNA 属常规宽松档，预期主要剔除线粒体比例高的尾部细胞。

输出
----
qc/scrna_qc_per_sample.tsv       每样本质控前后细胞数
qc/scrna_cluster_markers.tsv     每簇 top 30 markers
qc/scrna_cluster_annotation.tsv  每簇的指派细胞类型与打分
reference/reference_profile.tsv  细胞类型 × 基因 的参考谱（CPM）
reference/README.md              方法、版本与来源
"""

from __future__ import annotations

import argparse
import csv
import gzip
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# 公开的经典标志基因集（不含本项目数据的任何结果推导，不做主观增删）
MARKERS: dict[str, list[str]] = {
    "T_CD8": ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "GZMK", "CCL5", "GZMA"],
    "T_CD4": ["CD3D", "CD3E", "CD3G", "CD4", "IL7R", "CCR7", "SELL", "LEF1"],
    "T_reg": ["FOXP3", "IL2RA", "CTLA4", "TNFRSF4", "IKZF2", "TNFRSF18"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "PRF1", "TYROBP"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD19", "IGHM", "TCL1A", "BANK1"],
    "Plasma": ["MZB1", "JCHAIN", "DERL3", "IGHG1", "SDC1", "XBP1", "PRDM1"],
    "Macrophage": ["CD68", "CD163", "MRC1", "MSR1", "LST1", "C1QA", "C1QB", "AIF1"],
    "DC": ["CD1C", "FCER1A", "CLEC9A", "CST3", "HLA-DQA1", "HLA-DPA1", "ITGAX"],
    "pDC": ["LILRA4", "IRF7", "CLEC4C", "GZMB", "JCHAIN", "TCF4"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "EGFL7", "ERG", "RAMP2"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "COL3A1", "FAP"],
    "Melanocyte": ["MLANA", "PMEL", "TYR", "TYRP1", "DCT", "MITF", "S100B"],
}

MIN_GENES = 200
MIN_UMI = 500
MAX_MT = 20.0
N_PCS = 30
RESOLUTION = 1.0
SEED = 20261007


def extract_izar10x(tar_path: Path, out_root: Path, expect: int = 36) -> Path:
    """把 RAW.tar 里的 10x 三件套按样本解到 out_root/<GSM>_<patient>/ 下。"""
    out_root.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:*") as tf:
        for m in tf:
            if not m.isfile():
                continue
            gsm, rest = m.name.split("_", 1)
            patient = rest.split("_")[0]
            d = out_root / f"{gsm}_{patient}"
            d.mkdir(exist_ok=True)
            # CellRanger 需要目录内命名为 barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz
            kind = "barcodes" if rest.endswith("barcodes.tsv.gz") else \
                   "features" if rest.endswith("features.tsv.gz") else \
                   "matrix" if rest.endswith("matrix.mtx.gz") else None
            if kind is None:
                continue
            with open(d / f"{kind}.tsv.gz" if kind != "matrix" else d / "matrix.mtx.gz",
                      "wb") as out:
                out.write(tf.extractfile(m).read())
    n = sum(1 for _ in out_root.glob("*/*.mtx.gz"))
    print(f"解出 {n} 个样本的 matrix.mtx.gz（预期 {expect}）")
    return out_root


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tar", default="data/GSE294273/GSE294273_RAW.tar")
    ap.add_argument("--work-dir", default="work/scrna")
    ap.add_argument("--samples-tsv", default="meta/GSE294273_samples.tsv")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    work = outdir / args.work_dir
    qc = outdir / "qc"
    qc.mkdir(parents=True, exist_ok=True)

    extract_izar10x(Path(args.tar), work)

    samples = {r["gsm"]: r for r in csv.DictReader(
        Path(args.samples_tsv).open(encoding="utf-8"), delimiter="\t")}

    print("\n读入 10x 矩阵…")
    adatas = {}
    for d in sorted(work.glob("*")):
        gsm = d.name.split("_")[0]
        a = sc.read_10x_mtx(d, var_names="gene_symbols", make_unique=True, cache=False)
        a.obs["gsm"] = gsm
        a.obs["patient"] = d.name.split("_", 1)[1]
        a.obs["treatment"] = samples.get(gsm, {}).get("char::treatment", "")
        adatas[gsm] = a
        print(f"  {gsm} {a.n_obs:>6} 细胞 × {a.n_vars} 基因  [{a.obs['treatment'].iloc[0]}]")

    adata = sc.concat(adatas, label="cohort", index_unique=None, join="outer")
    adata.var_names_make_unique()
    # 各样本的 barcode 在各自文件里都是 `AAACCTGAGAAACCAT-1` 起的，跨样本必然重复。
    # 不去重的话 obs_names 撞车，聚类赋值与后续按细胞名索引都会错位。
    adata.obs_names_make_unique()
    n_dup_fixed = adata.obs_names.duplicated().sum()
    print(f"\n合并后 {adata.n_obs} 细胞 × {adata.n_vars} 基因，"
          f"{adata.obs['patient'].nunique()} 位供体；obs_names 去重后重复 {n_dup_fixed} 个")

    # ---- 质控 ----
    # calculate_qc_metrics(qc_vars=["mt"]) 要求 var 里先有 mt 布尔列，
    # 没有这一步会直接 KeyError: 'mt'
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    print(f"线粒体基因 {int(adata.var['mt'].sum())} 个")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    n0 = adata.n_obs
    adata = adata[
        (adata.obs.n_genes_by_counts >= MIN_GENES)
        & (adata.obs.total_counts >= MIN_UMI)
        & (adata.obs.pct_counts_mt <= MAX_MT)
    ].copy()
    print(f"质控（n_genes≥{MIN_GENES}, n_umi≥{MIN_UMI}, mt≤{MAX_MT}%）："
          f"{n0} → {adata.n_obs}，剔除 {n0 - adata.n_obs}（{(n0 - adata.n_obs) / n0 * 100:.1f}%）")

    # AnnData 没有 groupby，逐样本统计走 obs
    rows = []
    for gsm, sub in adata.obs.groupby("gsm", observed=True):
        rows.append({
            "gsm": gsm, "patient": sub["patient"].iloc[0],
            "treatment": sub["treatment"].iloc[0],
            "n_cells_after_qc": int(len(sub)),
            "n_genes_median": round(float(sub["n_genes_by_counts"].median()), 1),
            "n_umi_median": round(float(sub["total_counts"].median()), 1),
            "mt_pct_median": round(float(sub["pct_counts_mt"].median()), 2),
            "mt_pct_p95": round(float(sub["pct_counts_mt"].quantile(0.95)), 2),
        })
    qc_df = pd.DataFrame(rows).sort_values("gsm")
    qc_df.to_csv(qc / "scrna_qc_per_sample.tsv", sep="\t", index=False)
    print(qc_df.to_string(index=False))

    # ---- 标准化 / HVG / PCA / 聚类 ----
    print("\n标准化与降维…")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="gsm")
    # 必须先切到 HVG 再 scale：全矩阵 85k × 33.5k 稠密化后是 22 GB，
    # scale 会把稀疏矩阵稠密化，既慢又没必要。
    adata_hv = adata[:, adata.var.highly_variable].copy()
    print(f"  HVG {int(adata.var.highly_variable.sum())} 个 → scale/PCA 子集 "
          f"{adata_hv.n_obs} × {adata_hv.n_vars}")
    sc.pp.scale(adata_hv, max_value=10)
    sc.tl.pca(adata_hv, n_comps=N_PCS, svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(adata_hv, n_neighbors=15, n_pcs=N_PCS, random_state=SEED)
    sc.tl.leiden(adata_hv, resolution=RESOLUTION, key_added="cluster",
                 random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
    print(f"聚类完成：{adata_hv.obs['cluster'].nunique()} 个簇（resolution={RESOLUTION}, seed={SEED}）")

    adata.obs = adata_hv.obs.copy()
    adata.obsm["X_pca"] = adata_hv.obsm["X_pca"]
    adata.obsp["connectivities"] = adata_hv.obsp["connectivities"]
    adata.obsp["distances"] = adata_hv.obsp["distances"]
    adata.uns["neighbors"] = adata_hv.uns["neighbors"]

    adata.obs.to_csv(qc / "scrna_obs.tsv", sep="\t")
    adata.write_h5ad(outdir / "work/scrna_annotated.h5ad", compression="gzip")

    # ---- top markers ----
    # 只在 HVG 上做差异检验：全矩阵 33538 个基因 × 85k 细胞 × ~20 簇，
    # wilcoxon 实测跑 10 分钟仍未出结果。marker 只用于人工复核簇身份，
    # 用 3000 个 HVG 完全够，且不改变任何聚类结果与细胞类型指派。
    print(f"\n计算各簇 top markers（仅 HVG {int(adata.var.highly_variable.sum())} 基因）…")
    adata_mk = adata[:, adata.var.highly_variable].copy()
    sc.tl.rank_genes_groups(adata_mk, "cluster", method="wilcoxon", n_genes=30)
    mk = adata_mk.uns["rank_genes_groups"]
    out = []
    for cl in mk["names"].dtype.names:
        names = list(mk["names"][cl])
        scores = list(mk["scores"][cl])
        lfc = list(mk["logfoldchanges"][cl])
        for r, (g, s, l) in enumerate(zip(names, scores, lfc), 1):
            out.append({"cluster": cl, "rank": r, "gene": g,
                        "wilcoxon_z": round(float(s), 3), "log2FC": round(float(l), 3)})
    pd.DataFrame(out).to_csv(qc / "scrna_cluster_markers.tsv", sep="\t", index=False)
    print(f"各簇 top markers 已写入 qc/scrna_cluster_markers.tsv")

    # ---- 标志基因打分指派细胞类型 ----
    print("\n用标志基因集给每个簇打分…")
    present = {g: g in adata.raw.var_names for gs in MARKERS.values() for g in gs}
    missing = sorted(g for g, ok in present.items() if not ok)
    print(f"  标志基因 {sum(present.values())}/{len(present)} 个存在于矩阵"
          + (f"；缺失 {missing}" if missing else ""))

    # 只取标志基因那几十个做 scale。全矩阵 scale 会稠密化成 22 GB 纯属浪费，
    # score_genes 本来也只用这些基因。
    marker_genes = sorted({g for gs in MARKERS.values() for g in gs
                           if g in adata.raw.var_names})
    e = adata.raw[:, marker_genes].to_adata()
    e.obs = adata.obs[["cluster"]].copy()
    sc.pp.scale(e, max_value=10)
    scores = {}
    for ct, genes in MARKERS.items():
        use = [g for g in genes if g in e.var_names]
        if not use:
            continue
        sc.tl.score_genes(e, use, score_name=f"_s_{ct}", random_state=SEED)
        scores[ct] = e.obs[f"_s_{ct}"].values
    S = pd.DataFrame(scores, index=e.obs_names)
    S["cluster"] = e.obs["cluster"].values
    per_cluster = S.groupby("cluster", observed=True).mean()
    assign = per_cluster.idxmax(axis=1)
    margin = per_cluster.max(axis=1) - per_cluster.apply(
        lambda r: sorted(r.values)[-2], axis=1)

    ann = []
    for cl in per_cluster.index:
        row = {"cluster": cl, "assigned": assign[cl], "margin": round(float(margin[cl]), 3)}
        for ct in per_cluster.columns:
            row[f"score_{ct}"] = round(float(per_cluster.loc[cl, ct]), 3)
        ann.append(row)
    ann_df = pd.DataFrame(ann).sort_values("cluster", key=lambda s: s.astype(int))
    ann_df.to_csv(qc / "scrna_cluster_annotation.tsv", sep="\t", index=False)

    print(f"\n{'簇':>4} {'指派类型':<14} {'与次高差距':>10}")
    for _, r in ann_df.iterrows():
        print(f"{r['cluster']:>4} {r['assigned']:<14} {r['margin']:>10.3f}")

    print(f"\n质控与聚类产物已写入 {qc}/，AnnData 已保存 {outdir}/work/scrna_annotated.h5ad")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
