#!/usr/bin/env python3
"""任务 10：单细胞簇指派的客观 tie-break（`docs/03` §19，规则先于结果锁定）。

为什么需要
----------
`scripts/13` 用 `sc.tl.score_genes` 的簇均值取 argmax 指派细胞类型。该法用
**表达量纲的绝对值**，被少数高表达基因拉偏，41 个簇里有 6 个 margin < 0.3，
其中簇 9 实为 NK、簇 34 实为成纤维细胞，**Fibroblast 在参考谱中计数为 0**。

本脚本按 §19.2 锁定的主判据重做指派：**top-30 markers × marker 集的单侧
超几何富集**，并对无显著富集的簇给出明确的弃权标签 `Unassigned`。

**不重跑聚类**：41 簇的 Leiden 结果与 top-30 markers 来自 `scripts/13`，
分辨率与随机种子早已写定（§19.4 规定不因指派结果回头改聚类）。

输出
----
qc/scrna_cluster_tiebreak.tsv        每簇 × 每类型的富集结果与最终指派
qc/scrna_reference_profile_v2.tsv    参考谱组成（按细胞数加权）
qc/scrna_tiebreak_report.md
work/scrna_cluster_celltype_v2.tsv   每细胞一列的指派结果（供去卷积用）
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
QC = ROOT / "qc"
ALPHA = 0.05
TOP_N = 30

# 与 scripts/13 完全一致的 marker 集，不在此另写一份
pro13 = import_module("13_build_reference_profile")
MARKERS: dict[str, list[str]] = pro13.MARKERS


def detected_genes() -> list[str]:
    """被检测到的基因全集（`adata.n_vars`）。

    **必须用全集，不能用「top-30 里出现过的基因」**——后者已被差异表达
    筛选富集过，拿它当超几何零假设的总体会把富集人为压低（见 §19.2 第 3 条）。
    """
    import anndata as ad
    h5 = ROOT / "work" / "scrna_annotated.h5ad"
    a = ad.read_h5ad(h5, backed="r")
    genes = list(a.var_names)
    a.file.close()
    return genes


def load_cell_clusters():
    """从已落盘的 h5ad 里只读 obs（backed 模式，不载入 1.4 GB 表达矩阵）。

    若已存在 TSV 版本则优先用它。
    """
    tsv = ROOT / "work" / "scrna_cluster_celltype_v1.tsv"
    if tsv.exists():
        return pd.read_csv(tsv, sep="\t")
    h5 = ROOT / "work" / "scrna_annotated.h5ad"
    if not h5.exists():
        return None
    import anndata as ad
    a = ad.read_h5ad(h5, backed="r")
    obs = a.obs[["cluster"]].copy()
    celltype = "celltype" if "celltype" in a.obs.columns else None
    if celltype:
        obs["celltype"] = a.obs[celltype].astype(str).to_numpy()
    a.file.close()
    return obs.reset_index(names="cell")


def main() -> int:
    ann = pd.read_csv(QC / "scrna_cluster_annotation.tsv", sep="\t")
    mk = pd.read_csv(QC / "scrna_cluster_markers.tsv", sep="\t")

    # §19.2 第 3 条：超几何零假设的总体必须是**被检测到的基因全集**。
    # （初版实现误用了「top-30 出现过的基因」作总体，该总体已被差异表达
    #   富集过，会把富集人为压低，已修正——这是实现 bug，不是口径调整。）
    universe = set(detected_genes())
    N = len(universe)

    log: list[str] = []
    log.append("# 单细胞簇指派的客观 tie-break（`docs/03` §19）\n")
    log.append("生成脚本 `scripts/21_scrna_cluster_tiebreak.py`。")
    log.append("**规则写定于任何 tie-break 结果产生之前**；聚类与 top-30 markers "
               "沿用 `scripts/13` 已锁定的结果，不重跑、不改分辨率。\n")
    log.append(f"- 超几何检验总体：**{N} 个被检测基因全集**（§19.2 第 3 条）")
    log.append(f"- 每簇抽取 `top-{TOP_N}` 差异表达基因，检验其在 marker 集中的超几何富集")
    log.append(f"- 显著性阈值 α = {ALPHA}；未达阈值一律标 `Unassigned`\n")
    log.append("> **实现更正留痕**：本脚本初版把超几何总体错设为「top-30 里出现过的基因」"
               f"（{N} 的一个子集），该总体已被差异表达筛选富集过，会把富集人为压低，"
               "导致约四成细胞被弃权。已按 §19.2 第 3 条改回**被检测基因全集**。"
               "这是实现 bug 的修复，**不是**看到结果后改口径；规则本身一字未动。\n")

    rows: list[dict] = []
    assign: dict[int, str] = {}
    n_unassigned = 0

    for cl in sorted(ann["cluster"].astype(int)):
        top = set(mk.loc[mk["cluster"].astype(int) == cl, "gene"])
        old = ann.loc[ann["cluster"].astype(int) == cl].iloc[0]
        cands = []
        for ct, genes in MARKERS.items():
            present = [g for g in genes if g in universe]
            k = len(top & set(present))
            if k == 0:
                p = 1.0
            else:
                p = float(stats.hypergeom.sf(k - 1, N, len(present), len(top)))
            cands.append({"cluster": cl, "celltype": ct, "overlap": k,
                          "marker_n_present": len(present), "top_n": len(top),
                          "p_hyper": p,
                          "score_genes": float(old.get(f"score_{ct}", np.nan))})
        df = pd.DataFrame(cands)
        best = df.sort_values(["p_hyper", "overlap", "score_genes"],
                              ascending=[True, False, False]).iloc[0]
        sig = best["p_hyper"] < ALPHA
        final = best["celltype"] if sig else "Unassigned"
        if not sig:
            n_unassigned += 1
        assign[cl] = final
        for c in cands:
            c["final"] = final
            c["is_final_choice"] = (c["celltype"] == final)
            c["old_assigned"] = old["assigned"]
            c["old_margin"] = float(old["margin"])
        rows.extend(cands)

    detail = pd.DataFrame(rows)

    detail.to_csv(QC / "scrna_cluster_tiebreak.tsv", sep="\t", index=False)

    # ---- 与原指派对比 ----
    cmp_df = ann.copy()
    cmp_df["cluster"] = cmp_df["cluster"].astype(int)
    cmp_df["final"] = cmp_df["cluster"].map(assign)
    cmp_df["changed"] = cmp_df["assigned"] != cmp_df["final"]

    log.append("## 一、指派变化（只列发生变化的簇）\n")
    ch = cmp_df[cmp_df["changed"]].sort_values("cluster")
    if ch.empty:
        log.append("无变化。\n")
    else:
        log.append("| 簇 | 原指派 | 新指派 | 原 margin | 最优类型 | 命中基因 | 富集 p |")
        log.append("|---|---|---|---:|---|---:|---|")
        for _, r in ch.iterrows():
            sub = detail[detail["cluster"] == r["cluster"]]
            best = sub.sort_values(["p_hyper", "overlap", "score_genes"],
                                   ascending=[True, False, False]).iloc[0]
            note = "（未达 α，弃权）" if r["final"] == "Unassigned" else ""
            log.append(f"| {r['cluster']} | {r['assigned']} | **{r['final']}**{note} | "
                       f"{r['margin']:.3f} | {best['celltype']} | "
                       f"{int(best['overlap'])}/{int(best['top_n'])} | "
                       f"{best['p_hyper']:.2e} |")
    log.append("")

    # 逐细胞归属要在诊断段落里用到，先加载
    obs = load_cell_clusters()
    if obs is not None:
        obs["cluster"] = obs["cluster"].astype(int)
        obs["final"] = obs["cluster"].map(assign)
        frac_un = float((obs["final"] == "Unassigned").mean())
        frac_fib = float((obs["final"] == "Fibroblast").mean())
    else:
        frac_un = frac_fib = float("nan")

    log.append("### 1.1 弃权簇的 top-5 markers（说明为什么弃权）\n")
    log.append("弃权不是因为「没有生物学意义」，而是**当前 marker 集覆盖不到**。"
               "这一列是留给读者判断 marker 集是否需要扩充的唯一依据。\n")
    log.append("| 簇 | 原指派 | top-5 markers | 12 类 marker 集命中数 |")
    log.append("|---|---|---|---:|")
    un = cmp_df[cmp_df["final"] == "Unassigned"].sort_values("cluster")
    for _, r in un.iterrows():
        top5 = mk[(mk["cluster"].astype(int) == r["cluster"])
                  & (mk["rank"] <= 5)]["gene"].tolist()
        best = detail[detail["cluster"] == r["cluster"]].sort_values(
            ["p_hyper", "overlap"], ascending=[True, False]).iloc[0]
        log.append(f"| {r['cluster']} | {r['assigned']} | "
                   f"{', '.join(f'`{g}`' for g in top5)} | "
                   f"{int(best['overlap'])}/{int(best['top_n'])} |")
    log.append("")
    log.append("### 1.2 弃权簇的真实身份（逐个查证，不靠印象）\n")
    log.append("逐簇核对 top-30 后，6 个弃权簇分属三种完全不同的情况，"
               "**不能用一句「marker 集不够用」概括**：\n")
    log.append("| 情况 | 簇 | 判据 |")
    log.append("|---|---|---|")
    log.append("| **低质量 / ambient 簇** | 4、5、38 | top markers 是 "
               "`RPL34`/`RPS3A`/`RPS27A`（核糖体）、`MALAT1`/`MT-CO1`/`MT-CO2`、"
               "`GAPDH`/`IFI27` —— 无细胞身份信息，**弃权正确，且这些细胞本就不该进参考谱** |")
    log.append("| **增殖态细胞** | 24 | `HMGB2`/`MKI67`/`STMN1`/`TUBA1B`/`TMPO` —— "
               "是增殖细胞的标记，不是某个细胞**类型**的标记 |")
    log.append("| **marker 集覆盖不到的基质/血管亚型** | 13、34 | 见下 |")
    log.append("")
    log.append("**簇 34 的真实身份（此前一度被误记为成纤维细胞，此处更正）**："
               "其 top-30 为 `IGFBP2`、`AKAP12`、`POSTN`、`AEBP1`、**`COL6A1`**、"
               "`AFAP1`、**`COL4A2`**、`CALD1`、`SPARCL1`、`PTPRG`、`CAV1`、"
               "`ITGB5`、`VEGFA`、`MYH10` —— 这是**血管/周细胞样基质程序**，"
               "既不是成纤维细胞（无 `COL1A1`/`COL1A2`/`DCN`/`LUM`），"
               "也不是本项目 `Endothelial` 集所列的经典内皮"
               "（无 `PECAM1`/`VWF`/`CDH5`/`CLDN5`）。"
               "**12 个候选集全都覆盖不到它，因此判 `Unassigned` 是正确行为**——"
               "若沿用 `scripts/13` 的打分法，它会被错分给 Melanocyte。")
    log.append("")
    log.append("**簇 37 是真正的成纤维细胞**：top-30 含 `COL1A2`、`COL11A1`、`TNC`、"
               "`SFRP1`、`IGFBP7`、`A2M``，命中 `Fibroblast` 集的 `COL1A2`，"
               f"富集 p = {detail[(detail['cluster']==37)&(detail['celltype']=='Fibroblast')].iloc[0]['p_hyper']:.1e}。"
               "**参考谱中 Fibroblast 从 0% 恢复到 "
               f"{frac_fib:.1%}，正是这次修复的目的。**")
    log.append("")
    log.append("**但不回避它的后果**：`Unassigned` 占 "
               f"{frac_un:.1%}，"
               "其中含一整簇真实的血管/周细胞样基质细胞。"
               "去卷积时这部分表达会被归到其余谱上，"
               "**在补齐 marker 集之前，本参考谱不能声称是 12 类精细分解**。\n")

    n_ch = int(cmp_df["changed"].sum())
    log.append(f"**{n_ch}/{len(cmp_df)} 个簇的指派发生变化，"
               f"{n_unassigned} 个被判为 `Unassigned`。**\n")

    # ---- 参考谱组成 ----
    log.append("## 二、参考谱组成变化\n")
    comp_old = comp_new = None
    if obs is not None:
        n_cells = len(obs)
        # 原指派：h5ad 的 obs 里若没有 celltype 列，用 scripts/13 的结果回填
        if "celltype" in obs.columns:
            comp_old = obs["celltype"].value_counts(normalize=True)
        else:
            comp_old = (obs["cluster"].map(
                ann.set_index(ann["cluster"].astype(int))["assigned"])
                .value_counts(normalize=True))
        comp_new = obs["final"].value_counts(normalize=True)
        cmp2 = pd.concat(
            [comp_old.rename("old_frac"), comp_new.rename("new_frac")],
            axis=1).fillna(0).sort_values("new_frac", ascending=False)
        log.append("| 细胞类型 | 原占比 | 新占比 | 变化 |")
        log.append("|---|---:|---:|---:|")
        for ct in cmp2.index:
            log.append(f"| {ct} | {cmp2.loc[ct,'old_frac']:.1%} | "
                       f"**{cmp2.loc[ct,'new_frac']:.1%}** | "
                       f"{cmp2.loc[ct,'new_frac']-cmp2.loc[ct,'old_frac']:+.1%} |")
        log.append("")
        pd.DataFrame({"celltype": cmp2.index,
                      "old_frac": cmp2["old_frac"].to_numpy(),
                      "new_frac": cmp2["new_frac"].to_numpy()}).to_csv(
            QC / "scrna_reference_profile_v2.tsv", sep="\t", index=False)
        obs[["cluster", "final"]].rename(columns={"final": "celltype"}).to_csv(
            ROOT / "work" / "scrna_cluster_celltype_v2.tsv", sep="\t", index=False)
        log.append(f"参考谱细胞总数 {n_cells}；`Unassigned` 的细胞在去卷积时"
                   "按**未知组分**处理，**不得就近指派**（§19.4）。\n")
    else:
        log.append("未找到逐细胞的簇归属文件，本节跳过"
                   "（`scripts/13` 的 obs 未落盘时属正常，不影响指派结果本身）。\n")

    log.append("## 三、这个改动影响什么、不影响什么\n")
    log.append("| 项 | 是否受影响 |")
    log.append("|---|---|")
    log.append("| 去卷积参考谱（GSE294273 → bulk 分解） | **受影响**，这就是本次修复的目的 |")
    log.append("| 模块 2 肿瘤内在状态打分（GSE308433/434） | **不受影响**——该模块只用肿瘤核，不做细胞类型指派 |")
    log.append("| 模块 2b 基因剂量（GSE308433） | **不受影响**——纯 CNV 统计量，与细胞组成无关 |")
    log.append("| IMPRES 与建模（`docs/05`–`docs/07`） | **不受影响**——纯 bulk 表达 |")
    log.append("")
    log.append("**已记录的残余限制**：富集检验只解决**指派**，不改变**聚类**。"
               "若某类细胞在转录组上与其邻近类型本就连续（如 CD8 T 与 NK），"
               "Leiden 会不会切出边界本身就有不确定性，"
               "这一点不因 tie-break 而消失（§19.4）。\n")

    rep = QC / "scrna_tiebreak_report.md"
    rep.write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))
    print(f"\n[OK] {QC/'scrna_cluster_tiebreak.tsv'}")
    print(f"[OK] {rep}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())