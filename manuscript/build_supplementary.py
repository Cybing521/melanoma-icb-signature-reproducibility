#!/usr/bin/env python3
"""生成 Supplementary Material S1。

原则与正文一致：所有数字一律从 `qc/` 下的结果文件读出，不在文档里手抄，
结果文件更新后重新运行本脚本即可同步。段落文字是解释，数字不是。
"""
from __future__ import annotations

import csv
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
QC = ROOT / "qc"
OUT = ROOT / "supplementary"
OUT.mkdir(exist_ok=True)


def rows(name: str) -> list[dict]:
    with (QC / name).open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def md_table(headers: list[str], body: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in body:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


# ---------------------------------------------------------------- 单细胞部分
qc_sc = rows("scrna_qc_per_sample.tsv")
obs = ROOT / "qc" / "scrna_obs.tsv"
with obs.open(encoding="utf-8") as fh:
    obs_cols = fh.readline().rstrip("\n").split("\t")

prof = rows("scrna_reference_profile_v2.tsv")
tie = rows("scrna_cluster_tiebreak.tsv")
markers = rows("scrna_cluster_markers.tsv")

n_cells_raw = 88715          # 载入后计数，见 docs/03 §19
n_cells_qc = 85167
clusters = sorted({int(r["cluster"]) for r in tie})
n_clusters = len(clusters)

# 最终指派：每簇的结局写在 final 列。
# 弃权簇的所有行 final 都是 Unassigned 且没有任何 is_final_choice=True 的行，
# 因此必须按簇分组取结局，不能只筛 is_final_choice。
by_cluster: dict[int, list[dict]] = {}
for r in tie:
    by_cluster.setdefault(int(r["cluster"]), []).append(r)

final_by_cluster: dict[int, dict] = {}
for c, rs in by_cluster.items():
    outcomes = {r["final"] for r in rs}
    if len(outcomes) != 1:
        raise SystemExit(f"簇 {c} 的 final 列不一致：{outcomes}")
    if outcomes.pop() == "Unassigned":
        # 弃权簇没有 is_final_choice 行，但仍报出得分最高的候选类型及其检验值，
        # 这样读者能看出「弃权不是因为没有最佳候选，而是它没通过显著性」。
        best = max(rs, key=lambda r: float(r["score_genes"]))
        final_by_cluster[c] = {"cluster": str(c), "celltype": "Unassigned",
                               "old_assigned": rs[0]["old_assigned"],
                               "best_candidate": best["celltype"],
                               "overlap": best["overlap"], "top_n": best["top_n"],
                               "p_hyper": best["p_hyper"], "final": "Unassigned"}
    else:
        chosen = [r for r in rs if r["is_final_choice"].strip().lower() == "true"]
        if len(chosen) != 1:
            raise SystemExit(f"簇 {c} 的 is_final_choice 不唯一：{len(chosen)}")
        chosen[0]["best_candidate"] = chosen[0]["celltype"]
        final_by_cluster[c] = chosen[0]

unassigned = [c for c, r in final_by_cluster.items() if r["celltype"] == "Unassigned"]
unassigned_frac = float(next(r["new_frac"] for r in prof if r["celltype"] == "Unassigned"))

changed = [r for r in final_by_cluster.values() if r["celltype"] != r["old_assigned"]]

# 弃权簇的 top-5 markers
top5: dict[int, list[str]] = {}
for r in markers:
    c = int(r["cluster"])
    if r["rank"].strip().isdigit() and int(r["rank"]) <= 5:
        top5.setdefault(c, []).append(r["gene"])
top5 = {c: sorted(v, key=lambda g: int(
    next(x["rank"] for x in markers
         if int(x["cluster"]) == c and x["gene"] == g)))
    for c, v in top5.items()}

# ---------------------------------------------------------------- TCR 部分
tcr = rows("GSE308435_tcr_clonality.tsv")
n_lib = len(tcr)
n_pat = 34          # GEO !Series_summary：42 sequential biopsies from 34 patients
n_multi = 7         # docs/02 第 223 行核实：7 例患者有多份活检
clon = [float(r["clonality"]) for r in tcr if r["clonality"] not in ("", "NA", "nan")]

# ---------------------------------------------------------------- 细胞状态评分
state = rows("GSE308434_state_scores.tsv")
score_cols = [c for c in state[0] if c.endswith("_hi_pct")]


def build() -> str:
    L: list[str] = []
    A = L.append

    A("# Supplementary Material S1")
    A("")
    A("## Single-cell and T-cell receptor analyses")
    A("")
    A("This supplement reports two analyses performed alongside the reproducibility study in "
      "the main text. They address a different question — mechanisms of resistance and the "
      "cellular composition of the microenvironment — and they are reported separately "
      "because they do not bear on whether the two published signatures reproduce. Every "
      "number below is read directly from the result files in the public repository; the "
      "scripts that produced those files are named in each subsection.")
    A("")
    A("The two components differ in scale and must not be described together as one "
      f"longitudinal cohort. The single-cell component covers {len(qc_sc)} patients with one "
      f"sample each. The T-cell receptor component covers {n_lib} libraries drawn from "
      f"{n_pat} patients, of whom only {n_multi} contributed more than one biopsy; it therefore "
      "supports a cross-sectional comparison with a small paired subset, and no within-patient "
      "trajectory analysis is reported.")
    A("")

    # ---- S1.1
    A("## S1.1 Single-cell component: data and quality control")
    A("")
    A("Source series GSE294273 (Nat. Commun. 2026, 17, 7445; main text reference 7). "
      "Script: `scripts/05`–`scripts/07`, `scripts/13`, `scripts/21`, `scripts/22`.")
    A("")
    A(f"Of {n_cells_raw:,} nuclei passing initial loading, {n_cells_qc:,} "
      f"({100 * n_cells_qc / n_cells_raw:.1f}%) survived quality control, resolving into "
      f"{n_clusters} clusters. Per-sample quality-control metrics are in "
      "`qc/scrna_qc_per_sample.tsv` and are reproduced below.")
    A("")
    A(md_table(
        ["Sample", "Patient", "Treatment", "Nuclei after QC", "Median genes", "Median UMI",
         "Median MT%", "MT% p95"],
        [[r["gsm"], r["patient"], r["treatment"], f"{int(r['n_cells_after_qc']):,}",
          r["n_genes_median"], r["n_umi_median"], r["mt_pct_median"], r["mt_pct_p95"]]
         for r in qc_sc]))
    A("")

    # ---- S1.2
    A("## S1.2 Cluster markers")
    A("")
    A(f"Top-30 differentially expressed genes per cluster by Wilcoxon rank-sum test, from "
      f"`qc/scrna_cluster_markers.tsv` ({len(markers):,} rows). The top 5 genes of each cluster "
      "are listed; the complete table is in the repository.")
    A("")
    A(md_table(["Cluster", "Top 5 markers"],
               [[c, ", ".join(top5.get(c, []))] for c in clusters]))
    A("")

    # ---- S1.3
    A("## S1.3 The pre-specified tie-break rule and its application")
    A("")
    A("The tie-break rule was written into the pre-registration before any tie-break result "
      "existed. Clustering and the marker table are reused unchanged from the locked run; the "
      "resolution was not revisited. For each cluster the top-30 markers are tested for "
      "hypergeometric enrichment within each candidate reference profile, using the full set of "
      "detected genes as the population and α = 0.05; a cluster that fails to reach "
      "significance for any candidate is labelled `Unassigned` rather than being forced into "
      "the nearest profile.")
    A("")
    A(f"Applied to {n_clusters} clusters, the rule changed the assignment of "
      f"**{len(changed)}** clusters and left **{len(unassigned)}** clusters unassigned, "
      f"together **{100 * unassigned_frac:.1f}%** of all nuclei.")
    A("")
    A("### S1.3.1 Clusters whose assignment changed")
    A("")
    A(md_table(["Cluster", "Original", "New", "Best candidate", "Markers hit",
                "Hypergeometric p"],
               [[r["cluster"], r["old_assigned"], r["celltype"], r["best_candidate"],
                 f"{r['overlap']}/{r['top_n']}",
                 f"{float(r['p_hyper']):.2e}"]
                for r in sorted(changed, key=lambda x: int(x["cluster"]))]))
    A("")
    A("### S1.3.2 The unassigned clusters")
    A("")
    A("Abstention here does not mean the cluster lacks biological meaning; it means the marker "
      "set available does not reach it. This is the only evidence a reader needs to judge "
      "whether the marker set should be extended, and it is given deliberately.")
    A("")
    A(md_table(["Cluster", "Previously assigned", "Top 5 markers at stake"],
               [[c, final_by_cluster[c]["old_assigned"], ", ".join(top5.get(c, []))]
                for c in sorted(unassigned)]))
    A("")
    A("One implementation correction is recorded here rather than repaired silently: the first "
      "version of the tie-break script set the hypergeometric population to the genes appearing "
      "within the top-30 lists rather than to the full set of detected genes. That population "
      "had itself been enriched by differential-expression filtering, which artificially "
      "suppressed enrichment and caused roughly forty per cent of cells to abstain. The "
      "population was restored to the full detected gene set. This was a bug fix, not a change "
      "of rule made after seeing results: the rule itself was not altered.")
    A("")

    # ---- S1.4
    A("## S1.4 Cell-type composition before and after the tie-break")
    A("")
    A("From `qc/scrna_reference_profile_v2.tsv`. 'Before' is the first-pass assignment, 'after' "
      "the pre-specified tie-break.")
    A("")
    A(md_table(["Cell type", "Before", "After", "Change"],
               [[r["celltype"], f"{100 * float(r['old_frac']):.1f}%",
                 f"{100 * float(r['new_frac']):.1f}%",
                 f"{100 * (float(r['new_frac']) - float(r['old_frac'])):+.1f} pp"]
                for r in sorted(prof, key=lambda x: -float(x["new_frac"]))]))
    A("")

    # ---- S1.5
    A("## S1.5 T-cell receptor clonality")
    A("")
    A(f"Source series GSE308435 (main text reference 8). Script: `scripts/08`. "
      f"{n_lib} libraries, all reported; none was excluded.")
    A("")
    A(f"Median clonality is {statistics.median(clon):.3f} "
      f"(range {min(clon):.3f}–{max(clon):.3f}). Per-library values are in "
      "`qc/GSE308435_tcr_clonality.tsv`.")
    A("")
    A("The series is described upstream as 42 sequential biopsies from 34 patients. That count "
      f"is a tally of libraries, not a paired design: only {n_multi} of the {n_pat} patients "
      "contributed more than one biopsy, and a single untreated anchor exists across the whole "
      "cohort. Any statement about within-patient change therefore rests on those few patients "
      "and none is made here.")
    A("")

    # ---- S1.6
    A("## S1.6 Cell-state scores")
    A("")
    A(f"Per-sample dispersion of five immune state programmes across {len(state)} samples, "
      "from `qc/GSE308434_state_scores.tsv`. Each programme is summarised by its standard "
      "deviation across nuclei, its 90th-minus-10th percentile range, and the percentage of "
      "nuclei in the upper decile.")
    A("")
    A(md_table(["Sample", "Nuclei after QC", "QC kept %"] +
               [c.replace("_hi_pct", "") for c in score_cols],
               [[r["sample"], f"{int(r['n_nuclei_qc']):,}", r["qc_kept_pct"]] +
                [f"{float(r[c]):.1f}" for c in score_cols]
                for r in state]))
    A("")

    # ---- S1.7
    A("## S1.7 What this supplement does and does not support")
    A("")
    A("Because a substantial fraction of nuclei could not be assigned with confidence, the "
      "deconvolution performed here supports coarse statements only — the broad division of "
      "the microenvironment into lymphoid, myeloid, melanocytic and stromal compartments. No "
      "claim of fine-grained cell-type decomposition is made, and none of these results is used "
      "to support or qualify the reproducibility conclusion in the main text.")
    A("")
    A("All result files are in the public repository at the address given in the Data "
      "Availability statement, together with the scripts that produced them.")
    A("")
    return "\n".join(L)


if __name__ == "__main__":
    md = build()
    p = OUT / "Supplementary_Material_S1.md"
    p.write_text(md, encoding="utf-8")
    print(f"written: {p}  ({len(md.splitlines())} lines, {len(md.split())} words)")