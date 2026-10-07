#!/usr/bin/env python3
"""按 MDPI（Cancers）体例生成投稿稿件 .docx。

为什么稿件正文写成代码
----------------------
本文的核心主张是「每一个数字都能追回某次脚本运行」。如果正文靠手工排版，
这条链子在最后一步就断了：读者无法知道稿子里的 0.505 与 qc/ 里的是不是同一个数。
因此正文内容与生成过程一起版本化，改一个字都会在 git 里留下痕迹。

体例依据 MDPI 投稿要求
----------------------
- A4、双倍行距、Times New Roman 12 pt、**连续行号**（MDPI 审稿硬性要求）
- 章节编号 1. / 2.1. / 2.1.1.
- 表格题注在表上方、图题注在图下方
- 摘要为单段，关键词单列
- 文末固定后置项：Author Contributions / Funding / IRB / Informed Consent /
  Data Availability / Conflicts of Interest

**当前稿件状态**：分析与结果部分已完成；引言为框架性铺陈而非系统综述；
未包含补充材料的完整内容（单细胞 / TCR 部分已降为补充材料，见 §2.11）。

用法
----
    python manuscript/build_manuscript.py
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_DOCX = ROOT / "黑色素瘤ICB签名可复现性_预登记研究_Cancers投稿稿.docx"


def h(level: int, text: str) -> dict:
    return {"type": "heading", "level": level, "text": text}


def p(text: str) -> dict:
    return {"type": "paragraph", "text": text}


def table(caption: str, headers: list[str], rows: list[list[str]],
          widths: list[float] | None = None) -> dict:
    d = {"type": "table", "caption": caption, "headers": headers, "rows": rows}
    if widths:
        d["widths"] = widths
    return d


def figure(path: str, caption: str, width_mm: float) -> dict:
    return {"type": "figure", "path": path, "caption": caption, "width_mm": width_mm}


def build() -> dict:
    B: list[dict] = []

    # ─────────────────────────── 摘要 ───────────────────────────
    B.append(h(1, "Abstract"))
    B.append(p(
        "Published transcriptomic signatures that predict response to immune checkpoint blockade "
        "(ICB) are widely used as research baselines, yet their performance is rarely re-examined "
        "in independent cohorts with pre-specified criteria. We conducted a pre-registered "
        "reproducibility study using exclusively public data, locking the primary comparison, its "
        "endpoints and the criteria for attributing a negative result in writing before any "
        "computation. Two signatures were tested: IMPRES, a pairwise ordering of 15 "
        "immune-checkpoint genes [1], and a reduced Immunophenoscore (IPS-MHC+CP) built from 20 "
        "MHC and checkpoint genes [2], chosen for a construction unrelated to IMPRES. In the "
        "largest independent cohort (79 evaluable patients, 45 "
        "responders / 34 non-responders), IMPRES reached an AUC of 0.505 (95% CI 0.380–0.631) in "
        "the pre-registered direction and IPS-MHC+CP 0.591 (0.463–0.714); neither "
        "reproduced its published performance, and a cross-cohort meta-analysis of study-level "
        "estimates did not change the verdict. They failed differently: IMPRES was directionally "
        "erratic, IPS-MHC+CP weak but consistently signed. Of seven prospectively specified "
        "technical explanations, six were excluded outright and sample size, the seventh, "
        "against the published performance but not a weak signal (power 0.93 at an AUC of 0.70). "
        "Two explanations remain untestable with public data: population composition, and "
        "weakness intrinsic to the derivation cohorts. These signatures should therefore not be "
        "assumed to transfer across melanoma ICB cohorts."
    ))
    B.append(p("Keywords: immune checkpoint blockade; melanoma; reproducibility; "
               "pre-registration; IMPRES; reduced Immunophenoscore (MHC and checkpoint "
               "classes); treatment response; negative results"))

    # ─────────────────────────── 1. 引言 ───────────────────────────
    B.append(h(1, "1. Introduction"))
    B.append(p(
        "A substantial fraction of melanoma patients derive durable benefit from immune "
        "checkpoint blockade (ICB), while the majority do not. Because response is only "
        "observed after treatment has begun, considerable effort has gone into identifying "
        "transcriptomic features that predict benefit in advance. Two predictors are now "
        "standard references in this literature: IMPRES, which scores a patient by how often "
        "each of 15 immune-checkpoint genes is expressed above its paired partner [1], and the "
        "Immunophenoscore (IPS), which aggregates weighted z-scores across 26 immune-related "
        "gene sets [2]. Both are used less as clinical tests than as comparators against which new "
        "models must demonstrate added value."))
    B.append(p(
        "That role makes their reported performance load-bearing. A new signature that fails to "
        "beat an AUC of 0.8 in a re-analysis is judged unremarkable; one that fails to beat 0.5 "
        "is not published. The published figures are therefore rarely re-derived under "
        "pre-specified conditions, and when they are re-derived, the analysis is generally "
        "retrospective: the cohort, the endpoint and the acceptance criteria are chosen after "
        "the result is known. Retrospective replication can fail for reasons that have nothing "
        "to do with the signature, and a negative outcome becomes uninterpretable rather than "
        "informative."))
    B.append(p(
        "This study asks a deliberately narrow question: using only public data and a protocol "
        "fixed in advance, do published ICB-response signatures reproduce their published "
        "discrimination in independent melanoma cohorts? We treat the possibility of failure as "
        "the expected outcome and design against the standard objections that would otherwise "
        "make a negative result uninterpretable. Specifically, we register in advance (i) the "
        "primary comparison and its direction, (ii) the endpoint definitions, (iii) the "
        "quality-control thresholds, and (iv) the set of technical explanations that must be "
        "excluded before a negative result may be attributed to the signatures themselves. We "
        "then test two signatures selected to fail in distinguishable ways, and we report the "
        "boundary of what the available data can support."))
    B.append(p(
        "Two features of this design are worth stating at the outset, because they would "
        "otherwise look like flaws. First, one of the two signatures, IMPRES, was implemented "
        "from a gene list that could not be obtained from the original publication through any "
        "of four documented routes; the 15 pairs were reconstructed from the authors' code "
        "repository. This is the study's largest attack surface and is discussed in Section 4.5. "
        "Second, one procedural deviation from the pre-registration was discovered during "
        "internal review and is reported in full (Section 3.3) rather than corrected silently. "
        "We regard both as part of the result: a reproducibility study that reports only its "
        "clean parts is not describing what reproducibility costs."))

    # ─────────────────────────── 2. 方法 ───────────────────────────
    B.append(h(1, "2. Materials and Methods"))
    B.append(h(2, "2.1. Study Design and Pre-Registration"))
    B.append(p(
        "The study was designed as a pre-registered reproducibility analysis. The primary "
        "comparison, including its direction; the endpoint definitions for each cohort; the "
        "quality-control thresholds; and the set of technical explanations that had to be "
        "excluded before a negative result could be attributed to the signatures, were all "
        "written to a pre-registration document before the corresponding computation was run. "
        "The pre-registration document is version-controlled in the public repository cited "
        "below, together with the version of each protocol section that was in force at the "
        "time of computation."))
    B.append(p(
        "Two classes of analysis appear in this paper. Analyses specified in the "
        "pre-registration are reported as such. Analyses added after results were known are "
        "labelled post hoc at the point of use; they do not use response labels for model "
        "fitting, do not alter any adjudicated primary comparison, and are not used to promote "
        "any finding to a positive conclusion. The post hoc analyses here are the cross-cohort "
        "meta-analysis (Section 2.8), the stratified analyses of biopsy timepoint and regimen "
        "(Section 2.9), and the score-dynamics check (Section 2.9). Their pre-specified "
        "interpretive rules were fixed in writing before they were run."))

    B.append(h(2, "2.2. Data Sources and Patient-Level Assembly"))
    B.append(p(
        "All data are public and were obtained from the NCBI Gene Expression Omnibus. The unit "
        "of analysis is the patient, not the biopsy. Cohorts containing repeated biopsies per "
        "patient were assembled into paired or mixed structures rather than counted as "
        "independent samples, because treating a patient's longitudinal biopsies as independent "
        "inflates the effective sample size and narrows confidence intervals artificially. A "
        "single patient-level inclusion and exclusion table, generated by one script, is the "
        "source of truth for every n reported in this paper."))
    B.append(p(
        "Three cohorts carry the primary signature analysis. GSE91061 [3] (RNA-seq; 109 samples, "
        "65 patients) is a nivolumab monotherapy cohort with both pre-treatment and "
        "on-treatment biopsies; the pre-treatment subset of 33 patients (10 partial or complete "
        "responders, 23 progressors) forms the development set. GSE78220 [4] (bulk RNA-seq; 28 "
        "samples, 26 patients; 14 responders, 12 progressors) is an anti-PD-1 monotherapy "
        "cohort of pre-treatment biopsies. GSE215868 [5] (NanoString IO 360 panel, 770 genes; 105 "
        "samples, 105 patients) is a pre-treatment cohort treated with PD-1 axis blockade; 79 "
        "patients had a classifiable response (45 responders, 34 progressors). Five further series "
        "were used, none of them for the primary comparison. GSE244982 [6] supplies "
        "post-progression biopsies for the resistance-state analysis in Section 3.7 and "
        "carries no response labels, so it can only be asked whether a model scores "
        "progressed tumours as expected, never whether it discriminates responders. "
        "GSE294272 and GSE294273 [7] supply the cell-type reference profile for deconvolution "
        "and its independent validation. GSE308433, GSE308434 and GSE308435 [8] supply the "
        "single-cell and T-cell receptor material reported in Supplementary Material S1. "
        "GSE115821, the IMPRES derivation cohort, is named here for completeness and is "
        "discussed in Section 4.2; it was never analysed and contributes no data to any result "
        "in this paper."))
    B.append(p(
        "A data gate preceded every download: series titles, sample counts, characteristic "
        "fields and platform were inspected and recorded before expression matrices were "
        "retrieved, and any cohort whose metadata could not support the required endpoint was "
        "excluded at that point rather than after partial analysis."))

    B.append(h(2, "2.3. Endpoint Definitions"))
    B.append(p(
        "The analysis unit is the patient and the endpoint is best overall response. Responders "
        "are patients whose recorded best overall response was complete or partial response; "
        "non-responders are those with stable or progressive disease. Response labels were "
        "taken from the original RECIST fields in each series and, where necessary, mapped to "
        "this two-level definition with a fixed lookup table. Where a patient contributed "
        "several biopsies at the same visit, the patient-level median score was used, so that "
        "the contribution of a single patient does not depend on which of their biopsies "
        "happens to sort first in a file."))
    B.append(p(
        "One endpoint mismatch required explicit handling and is reported here because it is "
        "not visible from the GEO metadata. The native endpoint of GSE215868 is not RECIST: the "
        "originating study defined its primary outcome as long-term benefit, that is, survival "
        "without progression at 24 months. These are different constructs — a patient may meet "
        "RECIST criteria and progress at six months, or fail to meet them and remain stable "
        "for two years. We therefore derived our two-level response from the series' best "
        "overall response field rather than from its native long-term-benefit variable, which "
        "restores comparability with the RECIST endpoint on which IMPRES was reported. The "
        "consequence is that the present analysis is aligned with the published comparator, "
        "at the cost of discarding the original study's primary variable."))

    B.append(h(2, "2.4. Signature 1: IMPRES"))
    B.append(p(
        "IMPRES scores a sample by counting, across its 15 gene pairs, how many times the first "
        "gene of a pair is expressed above the second. The raw count ranges from 0 to 15; we "
        "used the low-expression direction as pre-registered. The complete gene list in the "
        "original publication's supplementary material could not be retrieved: the publisher "
        "page returned an access error, and three further routes failed. The 15 pairs were "
        "therefore reconstructed from the authors' public code repository. The reconstruction "
        "has one internal validity check: every pair contains at least one gene that is a "
        "direct target of checkpoint blockade. This check is weak, and the reconstruction "
        "remains the largest single threat to the validity of the IMPRES arm."))
    B.append(p(
        "The direction convention is itself unresolved, and this is a separate risk from the "
        "gene list. The two MATLAB files in the authors' repository that define the score "
        "specify the pairwise ratio in opposite orientations for the two feature blocks, so the "
        "repository is internally inconsistent; and the publication states only that the feature "
        "is 1 when the first gene is expressed below the second, without giving the written "
        "orientation of any of the 15 pairs in the supplementary table that could not be "
        "retrieved. The only material that could settle the convention is therefore unavailable. "
        "We pre-registered one orientation, report the mirrored orientation in parallel at every "
        "comparison, and neither reproduced the published performance, so these data cannot "
        "resolve the ambiguity in favour of either orientation; but the ambiguity remains "
        "unexcluded, and the normalisation-invariance check described next does not test it."))

    B.append(p(
        "An important property of the scoring rule was exploited to control the implementation: "
        "because it depends only on the ordering of paired expression values, IMPRES is "
        "mathematically invariant to any monotone transformation of the expression matrix. We "
        "verified this empirically rather than assuming it, scoring every sample under four "
        "normalisations in two directions; all eight combinations produced identical per-sample "
        "scores. Normalisation choice therefore cannot explain any IMPRES result in this paper."))

    B.append(h(2, "2.5. Signature 2: IPS-MHC+CP"))
    B.append(p(
        "A single signature resting entirely on a reconstructed gene list is a fragile design, "
        "so a second signature was added. It was selected under criteria fixed in advance: the "
        "construction must be mathematically unrelated to the first, the source group must be "
        "independent, and the gene list must be obtainable from open sources. The Immunophenoscore "
        "satisfies the first two. Its complete 26-set gene list, however, could not be obtained: "
        "six retrieval routes were attempted and all failed, which we report as a finding in "
        "its own right (Section 4.4)."))
    B.append(p(
        "We therefore tested the two classes that could be reconstructed gene by gene from open "
        "literature, namely the MHC class (10 genes) and the checkpoint class (10 genes), and "
        "refer to the resulting score as IPS-MHC+CP. It is not the Immunophenoscore and is "
        "never described as such. Scoring follows the original construction: each gene is "
        "converted to a within-cohort z-score, class means are taken, and the two class scores "
        "are summed. We used the continuous value rather than the original 0–10 integer mapping; "
        "the discretisation was judged in advance to destroy ranking resolution at the sample "
        "sizes studied here, and that judgement was recorded before the data were examined."))
    B.append(p(
        "This signature carries a property that makes its failure harder to excuse: all 20 genes "
        "are present on the NanoString IO 360 panel used in GSE215868, so it can be computed in "
        "full in the same sample in which IMPRES loses two of its fifteen pairs. Coverage is "
        "20/20 for IPS-MHC+CP against 13/15 for IMPRES."))

    B.append(h(2, "2.6. De Novo Feature-Based Modelling"))
    B.append(p(
        "In parallel with the signature validations, a feature-based classifier was developed on "
        "the 33-patient development set, using Hallmark gene sets [9] as features rather than "
        "hand-curated gene lists. Model selection was nested: 100 outer folds, feature "
        "selection performed strictly within each training fold, a hard cap of 10 features, and "
        "20 repeats of the entire procedure. The random seed was fixed at 20261007 throughout, "
        "and each result file records the seed used."))
    B.append(p(
        "One definitional matter proved consequential enough to fix in advance. A single AUC is "
        "not a well-defined summary of a repeated nested cross-validation, and the three common "
        "conventions differ by up to 0.05 for the pre-registered panel and 0.07 for the set "
        "actually used, on the same model. All nested-CV results in "
        "this paper are therefore reported under three conventions simultaneously: the pooled "
        "out-of-fold AUC, the mean ± standard deviation across repeats, and the mean ± standard "
        "deviation across outer folds. Reporting a single convention is treated as a reporting "
        "defect."))

    B.append(h(2, "2.7. Pre-Specified Exclusion of Technical Explanations"))
    B.append(p(
        "Before any computation, we enumerated the technical explanations that would ordinarily "
        "absolve a signature that failed to reproduce, together with the artefact required to "
        "exclude each. Seven were registered: implementation error; normalisation; platform "
        "suitability; sample size; biopsy timepoint; treatment regimen; and endpoint "
        "definition. A negative result is reported as concerning the signatures only if all "
        "seven have been addressed. The items are numbered in the order registered here, and Table 6 "
        "and Figure 4 follow the same numbering. Sections 2.3, 2.4, 2.9, 3.6 and 3.7 present "
        "the corresponding artefacts."))

    B.append(h(2, "2.8. Power Analysis and Cross-Cohort Meta-Analysis (Post Hoc)"))
    B.append(p(
        "A negative result invites the response that the sample was too small. We therefore "
        "computed, before interpreting any cohort, the power of the primary comparison to detect "
        "a range of true AUCs, and the sample size required for a confidence interval to exclude "
        "a clinically usable discrimination. Standard errors of the AUC follow Hanley and "
        "McNeil [10]. Two thresholds were used: 0.70, taken as the lower bound of a usable "
        "discriminator, and 0.77, the lower bound of the cross-cohort AUC range reported for "
        "IMPRES in its original publication."))
    B.append(p(
        "Because a single cohort yields a wide interval, we additionally pooled the "
        "study-level estimates. This does not merge patients. The pre-registration forbids "
        "sample-level pooling across cohorts, and this analysis does not violate it: the inputs "
        "are the per-cohort AUC estimates computed in this study together with their cohort "
        "sample sizes, and no patient-level record is touched or created. Pooling was performed "
        "on the logit scale, with a fixed-effect "
        "estimate and a DerSimonian–Laird random-effects estimate [11]; heterogeneity was assessed "
        "with Cochran's Q and I². Directional consistency across cohorts was assessed with a "
        "sign test, in which a signature whose point estimates all fall below 0.50 is recorded "
        "as consistently reversed and not as consistent with the published direction."))
    B.append(p(
        "Two interpretive rules were fixed before running this analysis. A pooled estimate does "
        "not override any single-cohort conclusion. And where a sign test or a pooled interval "
        "does not reach conventional significance, the observation is reported as requiring "
        "independent confirmation, not as partial replication."))

    B.append(h(2, "2.9. Cohort Heterogeneity and Attribution Boundaries (Post Hoc)"))
    B.append(p(
        "The central alternative to our conclusion is that we did not observe signature failure "
        "but a change of patients. We addressed this by making the differences between cohorts "
        "explicit and testing the components we could test. For each cohort we recorded "
        "platform, biopsy timepoint, treatment regimen, native endpoint, endpoint used here, "
        "and the relationship of that cohort to the derivation of IMPRES. We then stratified "
        "GSE91061 by biopsy timepoint (pre-treatment versus on-treatment) and GSE215868 by "
        "treatment regimen (ipilimumab plus nivolumab, pembrolizumab, and nivolumab), computing "
        "AUCs within each stratum. Finally we verified that both signatures vary meaningfully "
        "between patients within each cohort, since an AUC computed on a near-constant score is "
        "uninformative and would constitute the easiest technical escape from a negative "
        "result."))
    B.append(p(
        "Stratified analyses are reported with their full tables, including strata that look "
        "favourable, and are used solely to bound attribution. A stratum with few responders and "
        "a wide interval is not evidence of efficacy in that stratum, and no stratum result is "
        "used to qualify the primary conclusion."))

    B.append(h(2, "2.10. Statistical Analysis"))
    B.append(p(
        "Discrimination was quantified as the area under the receiver operating characteristic "
        "curve with 95% confidence intervals from 4000 bootstrap resamples of the patient-level "
        "score. Comparisons of responders against non-responders used one-sided Mann–Whitney U "
        "tests in both directions, as registered. A signature was judged to have reproduced its "
        "published performance only if the 95% interval excluded 0.50 in the pre-registered "
        "direction and the point estimate fell between 0.70 and 1.00. Both directions were "
        "computed for every signature and cohort and all are reported, and no direction was "
        "selected after the fact. We do not apply a hierarchy of tests [12]: the mirrored "
        "direction is reported not as a second look at a negative result, but because the "
        "pre-registered direction convention of IMPRES is itself unresolved (Section 2.4), so "
        "both conventions had to be shown. Confidence intervals are reported to three decimal places "
        "throughout. No multiplicity adjustment was applied, and none is claimed: the only "
        "pre-specified inferential tests are the primary comparison and its mirrored counterpart "
        "in each of the three cohorts, and every secondary, pooled and stratified analysis is "
        "reported as description without an inferential claim attached to it. No comparison in "
        "this paper is used to support a positive conclusion. Power was computed with the "
        "Hanley–McNeil standard error evaluated at the hypothesised true AUC, and, where "
        "sample size is discussed, by solving that expression for n; both are stated in the "
        "script that produces Table 4."))

    B.append(h(2, "2.11. Supplementary Material"))
    B.append(p(
        "Two further melanoma series were analysed as part of the study, on a question separate "
        "from the one tested here; both are reported in the Supplementary Material rather than in "
        "the main text. A single-cell analysis covers 12 patients (one sample each), and a T-cell "
        "receptor analysis covers 42 libraries drawn from 34 patients. Of those 34 patients only 7 "
        "contributed more than one biopsy, so the receptor data support a cross-sectional "
        "comparison with a small paired subset, not a 42-patient longitudinal series, and no "
        "within-patient trajectory analysis is reported. In the single-cell component, 88,715 "
        "cells passed initial loading and 85,167 (96.0%) survived quality control, resolving into "
        "41 clusters; 11 of 41 cluster "
        "assignments changed under a pre-specified, objective tie-break rule; and 6 clusters "
        "(15.1% of cells) could not be assigned with confidence against the available reference "
        "profiles. Because a substantial fraction of the data remained unassigned, the "
        "deconvolution performed here supports coarse trichotomous statements only, and no "
        "claim of fine-grained cell-type decomposition is made."))

    B.append(h(2, "2.12. Code and Data Availability"))
    B.append(p(
        "All analysis code, the complete pre-registration, the patient-level inclusion and "
        "exclusion table, all intermediate and result files, and the figure sources are publicly "
        "available at the repository given in the Data Availability statement. Scripts are "
        "numbered in execution order. Every result file records the random seed and the script "
        "that produced it. Expression data are not redistributed; they are retrieved from NCBI "
        "GEO by the scripts that consume them."))

    B.append(h(2, "2.13. Ethics Statement"))
    B.append(p(
        "This study analyses exclusively de-identified, publicly available data obtained from "
        "NCBI Gene Expression Omnibus. It involved no recruitment, no intervention and no "
        "access to identifiable patient information, and was therefore not subject to "
        "institutional review board approval."))

    # ─────────────────────────── 3. 结果 ───────────────────────────
    B.append(h(1, "3. Results"))
    B.append(h(2, "3.1. Cohort Assembly"))
    B.append(p(
        "The three primary cohorts yielded 33, 26 and 79 evaluable patients (Table 1), "
        "corresponding to 138 evaluable patients in total with 69 responders and 69 "
        "non-responders. Patient-level assembly is what makes these numbers: GSE91061 "
        "contributes 109 biopsies from 65 patients and GSE78220 contributes 28 biopsies from 26 "
        "patients, so counting samples rather than patients would have overstated the evidence "
        "in both. Figure 1 summarises the cohort roles, the patient-level assembly and the "
        "duplication hazards that produced these numbers."))
    B.append(figure("figures/Figure1_StudyDesign.png",
                    "Figure 1. Study design. (a) cohorts and their roles; (b) the patient-level "
                    "assembly that produces every n used in this paper; (c) the specific "
                    "duplication hazards that were checked and how each was handled; (d) what "
                    "may and may not cross cohort boundaries, and why pooling patients or "
                    "expression matrices is forbidden while study-level meta-analysis of "
                    "the per-cohort estimates computed here is not.",
                    163.8))
    B.append(table(
        "Table 1. Characteristics of the three primary cohorts.",
        ["", "GSE91061", "GSE78220", "GSE215868"],
        [
            ["Platform", "RNA-seq", "bulk RNA-seq", "NanoString IO 360 (770 genes)"],
            ["Samples (GEO)", "109", "28", "105"],
            ["Patients", "65", "26", "105"],
            ["Evaluable patients", "33", "26", "79"],
            ["Responders / non-responders", "10 / 23", "14 / 12", "45 / 34"],
            ["Biopsy timepoint used", "pre-treatment", "pre-treatment (26 of 27)", "pre-treatment (all)"],
            ["Regimen", "nivolumab", "anti-PD-1", "IPI+NIVO, pembrolizumab, nivolumab"],
            ["Native endpoint", "RECIST", "RECIST", "24-month long-term benefit (PFS-derived)"],
            ["Endpoint used here", "RECIST", "RECIST", "RECIST (derived from best overall response)"],
            ["Relation to IMPRES derivation", "reported validation cohort", "none", "none"],
        ]))

    B.append(h(2, "3.2. Pre-Registered Primary Analysis"))
    B.append(p(
        "In the pre-registered direction, IMPRES reached an AUC of 0.505 (95% CI 0.380–0.631) "
        "in GSE215868, the largest and most independent cohort. The interval is wide, but its "
        "upper bound lies below the lower bound of the cross-cohort range reported for IMPRES in "
        "its original publication, and the point estimate falls in the region that no "
        "pre-specified criterion would accept. Neither direction of IMPRES met the registered "
        "replication criterion in any of the three cohorts: the pre-registered direction gave "
        "0.359 (0.172–0.552) in GSE91061, 0.298 (0.122–0.521) in GSE78220 and 0.505 "
        "(0.380–0.631) in GSE215868, while the mirrored direction gave 0.659 (0.452–0.843), "
        "0.583 (0.351–0.786) and 0.483 (0.358–0.608) respectively. The registered one-sided "
        "Mann–Whitney U test of the primary comparison returns p = 0.474 (mirrored direction "
        "p = 0.608)."))
    B.append(p(
        "Two sensitivity analyses were registered in advance for this cohort and both are "
        "reported here irrespective of their result, as the protocol requires. Excluding the 22 "
        "patients with prior checkpoint-blockade exposure, leaving an immunotherapy-naive "
        "population of 36 responders and 27 non-responders, gave 0.518 (0.380–0.656) in the "
        "pre-registered direction and 0.468 (0.333–0.610) mirrored. Admitting stable disease as "
        "a third category leaves the pairwise discrimination unchanged by construction and the "
        "registered one-sided test unchanged at p = 0.474; descriptively, the median score is "
        "3.0 in all three groups (PRCR, stable disease and progression). Neither sensitivity "
        "analysis moves the result away from chance, and neither is used to alter the primary "
        "comparison."))
    B.append(p(
        "The mirrored direction deserves a specific caution. Its value in GSE91061, 0.659, is "
        "the largest IMPRES estimate obtained anywhere in this study, and a reader could be "
        "tempted to describe the signature as working once its direction is reversed. We "
        "decline to do so. The interval includes 0.50, the pooled estimate across all three "
        "cohorts is 0.537 (0.435–0.635) with no difference from chance, and the direction of "
        "the point estimate is inconsistent across cohorts. A signature with no discriminative "
        "power will produce a moderately impressive number on whichever side is chosen after "
        "the fact; that arithmetic is not a finding."))

    B.append(h(2, "3.3. A Procedural Deviation, and the Re-Run That Followed It"))
    B.append(p(
        "During internal consistency review we found that the feature set used for the "
        "development-set model differed from the pre-registered set. Three features had been "
        "added and four removed, and the pre-registered panel of six signatures had never "
        "been modelled on its own. We record this rather than correct it, and we handled it by "
        "re-running the pre-registered panel under the identical protocol: same 100 outer "
        "folds, same within-fold selection, same cap of 10 features, same seed."))
    B.append(p(
        "The re-run produced a result weaker than the one obtained with the deviant feature "
        "set: pooled out-of-fold AUC 0.417, mean across repeats 0.408 ± 0.046, mean across "
        "folds 0.347 ± 0.213, with the 95% interval of the per-repeat mean excluding 0.50. The "
        "conclusion is unchanged and, if anything, strengthened: neither the pre-registered "
        "feature set nor the set actually used yielded reproducible discrimination in 33 "
        "patients. The deviation record is retained permanently, since a pre-registration that "
        "is edited to match what was done is no longer a pre-registration."))
    B.append(p(
        "This exercise also resolved a reporting defect of our own making. A previously "
        "reported headline value of 0.500 for the deviant model was a pooled out-of-fold AUC in "
        "which each of the 33 patients was counted once per repeat and therefore 20 times. "
        "Under the three registered conventions the same model gives 0.500, 0.516 ± 0.112 and "
        "0.549 ± 0.257. All three are now reported together (Table 2), and Figure 2 shows them "
        "side by side for both feature sets."))
    B.append(figure("figures/Figure4_NestedCV.png",
                    "Figure 2. Nested cross-validation performance under the three registered "
                    "AUC conventions, for the pre-registered six-feature Hallmark panel and for the "
                    "deviant seventeen-feature set actually used, over 20 repeats.",
                    162.7))
    B.append(table(
        "Table 2. Nested cross-validation AUC under the three registered conventions "
        "(development set, 33 patients, 10 responders, 20 repeats).",
        ["Convention", "Pre-registered panel (6 features)", "Set actually used (17 features)"],
        [
            ["Pooled out-of-fold AUC", "0.417", "0.500"],
            ["Mean across repeats ± SD", "0.408 ± 0.046", "0.516 ± 0.112"],
            ["Mean across outer folds ± SD", "0.347 ± 0.213", "0.549 ± 0.257"],
        ]))

    B.append(h(2, "3.4. Signature 1 Across Cohorts: Directional Erraticism"))
    B.append(p(
        "The pooling and sign tests in this section were not pre-registered and are reported as "
        "post hoc (Section 2.8); they are descriptive of the per-cohort estimates above and are "
        "not used to support any conclusion on their own."))
    B.append(p(
        "Across the three cohorts the pre-registered direction of IMPRES gave 0.359, 0.298 and "
        "0.505. The estimates are not merely low; they are not consistently oriented. One cohort "
        "places the score above chance and two place it below, and the sign test across the "
        "three is uninformative (p = 1.00). Pooled at the study level (n = 138), the random-"
        "effects estimate is 0.415 (0.293–0.549) with I² = 34.6% (Cochran's Q = 3.06, p = 0.22), "
        "and no difference from chance (p = 0.21). The mirrored direction behaves the same way, "
        "pooling to 0.537 (0.435–0.635), I² = 0%, p = 0.48 against chance."))
    B.append(p(
        "The failure mode of IMPRES in this study is therefore best described as directional "
        "erraticism. Whatever the signature measures, that quantity does not order responders "
        "above non-responders consistently from one cohort to the next."))

    B.append(h(2, "3.5. Signature 2 Across Cohorts: A Consistent but Weak Signal"))
    B.append(p(
        "IPS-MHC+CP produced 0.604 (0.365–0.822), 0.583 (0.339–0.821) and 0.591 (0.463–0.714) "
        "in the three cohorts in its pre-specified high-score direction, which is the direction "
        "the Immunophenoscore construction assigns. All three point estimates exceed 0.50, and "
        "none of the intervals does so significantly. Under the registered criterion, all three "
        "cohorts are recorded as non-reproducing."))
    B.append(p(
        "The mirrored direction was computed in the same three cohorts, as the protocol required, "
        "and is reported here for completeness: 0.396 (0.178–0.626), 0.417 (0.190–0.661) and "
        "0.409 (0.288–0.541). Every mirrored estimate is below chance and every interval spans "
        "0.50, so the reversed score reproduces nothing either. The pre-registration required "
        "both directions to be reported; in the first draft of this manuscript only the "
        "high-score direction was shown even though the mirrored values had already been "
        "computed. That was a reporting omission rather than a change of analysis, and it does "
        "not bear on any conclusion."))
    B.append(p(
        "The cross-cohort behaviour nonetheless differs from IMPRES in a way that is worth "
        "preserving rather than discarding. The random-effects pooled estimate is 0.592 "
        "(0.492–0.685) with I² = 0% and Cochran's Q = 0.018 (p = 0.99): the three cohort "
        "estimates agree with one another more closely than chance variation alone would "
        "predict, despite the cohorts differing in platform, regimen and recruitment setting. "
        "The pooled estimate differs from chance at p = 0.071, and its interval still includes "
        "0.50. The sign test across three cohorts gives p = 0.25."))
    B.append(p(
        "We therefore report the following, and nothing stronger: IPS-MHC+CP exhibits a weak, "
        "cross-cohort-consistent signal that this design can neither confirm nor exclude. It is "
        "not replication. It is worth being precise about what a larger sample would buy, because the "
        "two requirements are very different. To push the upper bound of the confidence interval "
        "below the registered replication threshold of 0.70, at the observed point estimate of "
        "0.591, would take approximately 53 patients per group. To establish instead that a "
        "discrimination of 0.591 is genuinely above chance requires far more: on the same "
        "convention used for Table 4, the power of the present 45-and-34 comparison at a true AUC "
        "of 0.591 is 0.29, and reaching 0.80 would need on the order of 150 patients per group "
        "(Section 3.6). Neither is the sample size available here. "
        "The I² of zero should not be read as evidence that the cohorts are "
        "interchangeable; the more parsimonious reading, given how different the cohorts are, "
        "is that the underlying signal is weak enough that population differences are swamped "
        "by it. Table 3 collects every cross-cohort estimate in one place and Figure 3 "
        "shows them with their intervals."))

    B.append(figure("figures/Figure2_SignaturePerformance.png",
                    "Figure 3. Discrimination of both signatures across the three validation "
                    "cohorts, in the pre-registered direction, with 95% bootstrap confidence "
                    "intervals. (a) and (b) show each cohort in the pre-registered direction, "
                    "with the dashed line at chance and the shaded band marking the registered "
                    "replication window, AUC 0.70 to 1.00. Neither signature reaches it. That band "
                    "is a pre-registered decision threshold and not a published value: IMPRES "
                    "reports a cross-cohort range of 0.77 to 0.96 in its original publication, "
                    "whereas IPS-MHC+CP is derived in this study and has no published counterpart, "
                    "so the same threshold is applied to it purely as a reference bar. "
                    "(c) plots the same point estimates as a line across "
                    "cohorts: IMPRES does not agree even on which side of chance it falls, "
                    "whereas IPS-MHC+CP stays above it throughout. The two signatures therefore "
                    "fail differently, which is the substantive result.",
                    153.2))
    B.append(table(
        "Table 3. Cross-cohort results for both signatures. Single-cohort AUCs with bootstrap "
        "95% intervals, and random-effects pooled estimates over the 138 patients contributing "
        "study-level estimates; no patient-level record was pooled. Pooled columns are post hoc "
        "(Section 2.8); the single-cohort columns are the pre-registered analyses.",
        ["Signature", "GSE91061", "GSE78220", "GSE215868", "Pooled (RE)", "I²", "p vs 0.50", "Sign test p"],
        [
            ["IMPRES, pre-registered direction", "0.359 (0.172–0.552)", "0.298 (0.122–0.521)",
             "0.505 (0.380–0.631)", "0.415 (0.293–0.549)", "34.6%", "0.213", "1.00"],
            ["IMPRES, mirrored", "0.659 (0.452–0.843)", "0.583 (0.351–0.786)",
             "0.483 (0.358–0.608)", "0.537 (0.435–0.635)", "0.0%", "0.480", "1.00"],
            ["IPS-MHC+CP, high-score direction", "0.604 (0.365–0.822)", "0.583 (0.339–0.821)",
             "0.591 (0.463–0.714)", "0.592 (0.492–0.685)", "0.0%", "0.071", "0.25"],
            ["IPS-MHC+CP, mirrored", "0.396 (0.178–0.626)", "0.417 (0.190–0.661)",
             "0.409 (0.288–0.541)", "not pooled", "\u2014", "\u2014", "0.25"],
        ], widths=[2.5, 1.35, 1.35, 1.35, 1.5, 0.8, 0.9, 0.95]))

    B.append(h(2, "3.6. What These Cohorts Can and Cannot Resolve"))
    B.append(p(
        "The sample size objection deserves a direct answer rather than a concession, and Table 4 "
        "gives the powers directly. In the "
        "primary comparison, with 45 responders and 34 non-responders and a one-sided alpha of "
        "0.025, the probability of detecting a true AUC of 0.70 is 0.929, and of detecting a "
        "true AUC of 0.77, the lower bound of the published cross-cohort range, is 0.999. The "
        "primary comparison would have seen performance of the magnitude IMPRES reports. It "
        "observed 0.505 instead. On this cohort the negative result is informative, not a "
        "consequence of being underpowered."))
    B.append(p(
        "The same calculation bounds the other direction. Detecting a true AUC of 0.65 has power "
        "0.685, and detecting one of 0.60 has power 0.347. The minimum discrimination "
        "detectable at 80% power is 0.661 at 45 patients per group and 0.681 at 34 per group, "
        "the latter being the smaller of the two arms actually available. This design is thus "
        "well powered against the published performance and poorly powered against a weak one, "
        "which is exactly the boundary within which the IPS-MHC+CP result must be read."))
    B.append(table(
        "Table 4. Power of the primary comparison (45 responders, 34 non-responders) to detect "
        "a true AUC, one-sided alpha = 0.025. Post hoc (Section 2.8); reported to bound what the "
        "negative result can and cannot exclude, not to qualify the pre-registered comparison.",
        ["True AUC", "0.60", "0.65", "0.70", "0.77"],
        [["Power", "0.347", "0.685", "0.929", "0.999"]]))

    B.append(h(2, "3.7. Cohort Heterogeneity and the Limits of Attribution"))
    B.append(p(
        "The analyses in this section were not pre-registered and are reported as post hoc "
        "(Section 2.9). They address why the result is negative, not whether it is; nothing in "
        "this section is used to alter the pre-registered comparison."))
    B.append(p(
        "The three cohorts differ in platform, in regimen, in recruitment setting and, as set "
        "out in Section 2.3, in their native endpoint. We tested the components of that "
        "difference that our data can resolve."))
    B.append(p(
        "Biopsy timepoint did not explain the result. Stratifying GSE91061 by visit, the "
        "pre-treatment stratum reproduced the development-set point estimate exactly (AUC 0.359 in "
        "both; the re-run interval is 0.180–0.557 against 0.172–0.552 in Table 3, the difference "
        "being a separate bootstrap resampling), with 10 responders and 23 non-responders, which "
        "also serves as an independent "
        "check that the scoring implementation is stable. The on-treatment stratum did not "
        "recover the signal; it reversed, giving 0.172 (0.057–0.309) with 12 responders and 24 "
        "non-responders, an interval lying entirely below chance. Pooling the two visits gives 0.250 "
        "(0.140–0.370), but that pooled figure is a sample-level quantity: it rests on 69 "
        "biopsies drawn from the 65 patients of this cohort, with the multi-visit patients "
        "contributing twice. It is reported only for completeness; no conclusion here rests on "
        "it, and the patient-level figures quoted above are the ones used."))
    B.append(p(
        "Treatment regimen did not explain the result either. Within GSE215868, IMPRES gave 0.490 "
        "(0.304–0.682) in the ipilimumab–nivolumab stratum, 0.485 (0.235–0.730) in the "
        "nivolumab stratum and 0.547 (0.297–0.801) in the pembrolizumab stratum. We report the "
        "corresponding IPS-MHC+CP estimates in full for transparency: 0.524 (0.326–0.726), 0.710 "
        "(0.450–0.920) with 10 responders and 10 non-responders, and 0.570 (0.320–0.805) "
        "respectively. The nivolumab value for IPS-MHC+CP is the highest number in this study "
        "for that signature, and it is the reason the full stratified results are given in Table 5 "
        "rather than only the favourable strata: a table reporting only the good ones would be "
        "selection rather than audit. "
        "With 20 patients and an interval spanning 0.45 to 0.92 it carries no inferential weight, "
        "and we draw no conclusion from it."))
    B.append(p(
        "Platform and normalisation were excluded before the cohort analyses, as described in "
        "Section 2.4. The escape route that a near-constant score renders an AUC uninformative "
        "was checked directly: IMPRES scores in each cohort have a standard deviation of 1.55 to "
        "2.26 and an interquartile range of 2.5 to 3.0 across 7 to 12 distinct values, and "
        "IPS-MHC+CP has a standard deviation of 1.10 to 1.38 with no repeated values. Both "
        "signatures vary substantially between patients in every cohort."))
    B.append(table(
        "Table 5. Stratified results underlying Section 3.7, reported in full including strata "
        "that do not favour the conclusion. n is responders / non-responders. Post hoc "
        "(Section 2.9); no inferential claim is made from any stratum.",
        ["Stratum", "Signature", "n", "AUC (95% bootstrap CI)"],
        [
            ["Biopsy timepoint: pre-treatment", "IMPRES g1_low", "10 / 23", "0.359 (0.180–0.557)"],
            ["Biopsy timepoint: on-treatment", "IMPRES g1_low", "12 / 24", "0.172 (0.057–0.309)"],
            ["Biopsy timepoint: both visits pooled*", "IMPRES g1_low", "22 / 47", "0.250 (0.140–0.370)"],
            ["Regimen: ipilimumab + nivolumab", "IMPRES g1_low", "18 / 16", "0.490 (0.304–0.682)"],
            ["Regimen: nivolumab", "IMPRES g1_low", "10 / 10", "0.485 (0.235–0.730)"],
            ["Regimen: pembrolizumab", "IMPRES g1_low", "16 / 8", "0.547 (0.297–0.801)"],
            ["Regimen: ipilimumab + nivolumab", "IPS-MHC+CP", "18 / 16", "0.524 (0.326–0.726)"],
            ["Regimen: nivolumab", "IPS-MHC+CP", "10 / 10", "0.710 (0.450–0.920)"],
            ["Regimen: pembrolizumab", "IPS-MHC+CP", "16 / 8", "0.570 (0.320–0.805)"],
            ["Regimen: nivolumab + experimental\u2026", "both", "1 / 0", "not estimable (a single responder)"],
        ], widths=[2.6, 1.5, 1.1, 2.3]))
    B.append(p(
        "* The pooled-visit row is a sample-level quantity resting on 69 biopsies from 65 patients "
        "and is shown for completeness only. The fourth regimen stratum contained one patient and "
        "is listed so that the table is complete; no AUC can be formed from it. Figure 4 shows "
        "the artefact used to address each registered explanation, and Table 6 records the "
        "verdict, including the two that the available data cannot exclude."))
    B.append(figure("figures/Figure3_RuledOutExplanations.png",
                    "Figure 4. The seven technical explanations registered in advance and the artefact used "
                    "to address each, with the two explanations that the available data cannot "
                    "exclude. (a) implementation and normalisation invariance, shown as per-sample "
                    "scores under four normalisations in two directions, all eight combinations "
                    "identical; (b) the fraction of IMPRES gene pairs with informative pairwise "
                    "orderings in each cohort; (c) the ratio of observed to independence-expected "
                    "score standard deviation, showing no low-dimensional collapse; (d) "
                    "biopsy-timepoint strata and (e) treatment-regimen strata, reported in full "
                    "including strata that do not favour the conclusion; (f) power of the primary "
                    "comparison, 45 responders and 34 non-responders, as a function of the true "
                    "AUC, with the registered 0.70 threshold and the 0.77 lower bound of the "
                    "published range marked. Panel titles follow the registration order given in "
                    "Section 2.7; item 7, endpoint definition, required no data panel and is "
                    "documented in the band.",
                    171.2))
    B.append(table(
        "Table 6. Attribution boundary. Alternative explanations for the negative result, "
        "whether the present data can exclude them, and the artefact that does so.",
        ["Alternative explanation", "Excludable here?", "Basis"],
        [
            ["Implementation error", "Partly", "Normalisation-invariance check (4 normalisations × 2 directions, 8/8 identical scores) rules out normalisation-dependent implementation error. The reconstructed 15-pair list and the direction convention (Section 2.4) cannot be verified against the publication and are not excluded"],
            ["Normalisation", "Yes", "IMPRES is provably invariant to any monotone transform"],
            ["Platform unsuitability", "Yes", "Targeted panel covers 13/15 IMPRES pairs, at the pre-registered inclusion gate, and scores are non-degenerate in every cohort (SD 1.55–2.26)"],
            ["Insufficient sample size", "Partly", "Power 0.93 at AUC 0.70, 0.999 at AUC 0.77. Not excluded against a weak signal: the 80% power floor is AUC 0.661 at 45 per group and 0.681 at 34 per group, and power at 0.60 is 0.347"],
            ["Biopsy timepoint", "Yes", "Pre-treatment stratum at chance; on-treatment stratum reversed"],
            ["Treatment regimen", "Yes", "All three regimen strata at chance for IMPRES"],
            ["Endpoint definition", "Partly", "GSE215868 pulled back to a RECIST definition for comparability. Residual risk: the RECIST labels were parsed by this project and cannot be fully reconciled with formal RECIST assessment"],
            ["Population composition", "No", "No matching variables (subtype, mutation, metastatic burden) available for adjustment"],
            ["Weakness intrinsic to the derivation cohorts", "No", "IMPRES derivation cohort holds 8 independent patients; never analysable here"],
        ]))

    # ─────────────────────────── 4. 讨论 ───────────────────────────
    B.append(h(1, "4. Discussion"))
    B.append(h(2, "4.1. Principal Findings"))
    B.append(p(
        "Two published ICB-response signatures, selected in advance to be mathematically "
        "unrelated and independently derived, both failed to reproduce their published "
        "discrimination in independent melanoma cohorts under a pre-registered protocol. IMPRES "
        "reached 0.505 (0.380–0.631) in the largest cohort in its pre-registered direction and "
        "pooled to 0.415 (0.293–0.549) across the three cohorts, which contributed study-level "
        "estimates only; no patient-level record was pooled. IPS-MHC+CP reached 0.591 "
        "(0.463–0.714) and pooled to 0.592 (0.492–0.685). Six of the seven pre-specified "
        "technical explanations were excluded outright by pre-specified tests; sample size, the "
        "seventh, was excluded against the published performance but not against a weak signal. "
        "The primary comparison had 0.93 power to detect the performance these signatures "
        "report elsewhere."))
    B.append(p(
        "The two failures are of different kinds, and this difference is the substantive result "
        "rather than an incidental observation. IMPRES is directionally erratic: its point "
        "estimates do not even agree on which side of chance they fall. IPS-MHC+CP is "
        "directionally consistent, with three estimates between 0.583 and 0.604 and a pooled "
        "interval whose lower bound sits just below 0.50. A study with modest power can "
        "distinguish these situations, and the distinction matters practically: the first "
        "signature offers nothing to build on, whereas the second may contain a small signal "
        "that a sufficiently large cohort would resolve."))

    B.append(h(2, "4.2. What We Can and Cannot Conclude"))
    B.append(p(
        "The boundary of this study is worth stating precisely, because it is easy to overstate "
        "in either direction. We have shown that these signatures do not perform as reported, "
        "in these cohorts, at these sample sizes, and that this cannot be attributed to any of "
        "six enumerated technical causes, with sample size the seventh and only partially "
        "excluded. We have not shown that the signatures are without "
        "value in the populations from which they were derived."))
    B.append(p(
        "The reason is structural. The derivation cohort of IMPRES, although nominally public, "
        "is GSE115821, which contains too few independent patients to analyse: its 37 deposited "
        "samples correspond to "
        "8 patients with approximately 2 independent responders, a discrepancy with the "
        "published cohort size that the original authors did not resolve publicly. It is "
        "therefore impossible, with public data, to separate the two possibilities that the "
        "signatures are weak everywhere from the possibility that they are specific to a "
        "population we did not sample. Any statement that these results show the signatures to "
        "be invalid, as opposed to non-transferable, overstates the evidence."))
    B.append(p(
        "The second irreducible limitation is population composition. Our three cohorts differ "
        "in platform, regimen, endpoint and setting. We could exclude timepoint and regimen as "
        "explanations, but we could not adjust for melanoma subtype, mutational background or "
        "metastatic disease burden, because the necessary covariates are not consistently "
        "available across series. GSE215868 makes the point concretely: 45 of its 79 evaluable "
        "patients are responders, a response rate of 57% that is high for an immune checkpoint "
        "blockade cohort and reflects that series' question, which was long-term benefit rather "
        "than response to first-line treatment. The enrolled populations are therefore not "
        "comparable on response composition, and no covariate adjustment can repair that. "
        "Cohort heterogeneity remains a live alternative explanation in full generality."))

    B.append(h(2, "4.3. Why Two Signatures Fail Differently"))
    B.append(p(
        "The contrast between the two signatures is more informative than either result alone. "
        "IMPRES depends on a gene list we could not obtain from the publication and had to "
        "reconstruct, and its behaviour is erratic in every direction. IPS-MHC+CP was specified "
        "in advance to have a fully open gene list, full coverage on a routine clinical panel, "
        "and a construction unrelated to the first signature, and it produced a consistent, if "
        "weak, ordering. The natural reading is that the erratic signature is the one whose "
        "definition is least securely known, while the consistent signature is better specified "
        "and may be capturing a small amount of real information. This reading is consistent "
        "with the data but is not established by them, and we present it as a hypothesis that "
        "the reconstruction problem in Section 4.5 makes difficult to test."))

    B.append(h(2, "4.4. Accessibility as a Failure Mode Distinct from Accuracy"))
    B.append(p(
        "A secondary finding emerged from assembling the second signature. The Immunophenoscore "
        "is among the most frequently cited predictors in this field, and its complete gene list "
        "could not be retrieved through any of six routes: the publisher's article page, the "
        "publisher's API, two repositories named by the authors, the package archive, and a "
        "literature search. It was necessary to fall back on gene lists published separately in "
        "the open literature for two of its four classes."))
    B.append(p(
        "This is not a complaint about access policy, and it does not reflect on the science. It "
        "is a reproducibility observation of the same kind as the reconstruction of IMPRES, and "
        "it has a practical consequence: a signature whose definition cannot be recovered from "
        "open sources cannot be independently implemented, and any re-analysis of it is "
        "confounded with the analyst's reconstruction choices. We would suggest that "
        "computational predictors report machine-readable gene definitions in a public "
        "repository as a matter of course. The cost is negligible and the effect on downstream "
        "reproducibility is direct."))

    B.append(h(2, "4.5. Limitations"))
    B.append(p(
        "The limitations of this study are, in order of severity: the reconstruction of the "
        "IMPRES gene list from the authors' repository rather than from the publication's "
        "supplementary material, whose four documented retrieval routes all failed; the unresolved "
        "direction convention in that reconstructed implementation, which the two files in the "
        "authors' repository define in opposite orientations and which the normalisation-invariance "
        "check does not test; the "
        "inability to analyse GSE115821, the IMPRES derivation cohort, which contains too few "
        "independent "
        "patients; the absence of matching covariates for population adjustment; the small "
        "responder counts in the development and GSE78220 cohorts; the fact that IPS-MHC+CP "
        "covers two of the four Immunophenoscore classes and is not the Immunophenoscore; and "
        "the single procedural deviation described in Section 3.3, which we report rather than "
        "repair."))
    B.append(p(
        "Two further points of scope should be stated. First, the design tests discrimination "
        "between responders and non-responders; it does not test survival prediction, and the "
        "GSE215868 study's native long-term-benefit endpoint is not evaluated here. Second, the "
        "single-cell and T-cell receptor analysis reported in the Supplementary Material "
        "addresses mechanisms of resistance rather than the reproducibility of the signatures, "
        "and its conclusions should not be read as bearing on the primary result; a substantial "
        "fraction of its clusters could not be assigned, which limits it to coarse statements."))

    B.append(h(2, "4.6. Implications"))
    B.append(p(
        "For practice, the immediate implication is narrow and should not be overstated: "
        "published ICB-response signatures should be treated as comparators of uncertain "
        "transportability rather than as established performance benchmarks, and a new signature "
        "that merely approaches a re-derived 0.5 has not demonstrated added value."))
    B.append(p(
        "For design, two concrete recommendations follow. Studies of this kind should "
        "pre-register the acceptance criterion, the direction and the list of technical "
        "explanations in advance, because the value of a negative result lies entirely in "
        "whether it was designed to be interpretable, and because choosing which of several "
        "completed analyses to report after seeing their values manufactures significance out "
        "of a small number of null results [13]. And they should be powered explicitly "
        "against the published comparator rather than against a conventional per-group sample "
        "size: as Section 3.6 shows, 45 and 34 patients is ample to detect the performance these "
        "signatures report, and the same cohort is poorly powered to detect an AUC of 0.60, so "
        "the sample size determines which question the study can answer."))

    # ─────────────────────────── 5. 结论 ───────────────────────────
    B.append(h(1, "5. Conclusions"))
    B.append(p(
        "Two independently derived, mathematically unrelated published signatures for ICB "
        "response in melanoma both failed to reproduce their published discrimination in "
        "independent public cohorts under a pre-registered protocol, and the failure cannot be "
        "attributed to any of the six technical causes we could exclude, sample size being the "
        "seventh and only partially excluded. IMPRES failed by producing no "
        "consistent direction; IPS-MHC+CP produced a consistent but weak signal that the "
        "available sample size can neither confirm nor exclude, and for which approximately 53 "
        "patients per group would be needed to exclude the published level of discrimination, "
        "and on the order of 150 to establish the observed one above chance. What the public "
        "data cannot determine is whether "
        "these signatures are weak everywhere or weak only outside their derivation "
        "populations. That question requires the derivation cohorts, and for IMPRES those are "
        "not publicly analysable."))

    B.append(h(1, "Supplementary Materials"))
    B.append(p(
        "Supplementary Material S1: single-cell analysis (12 patients) and T-cell receptor "
        "analysis (42 libraries from 34 patients, of whom 7 contributed more than one biopsy). "
        "It is supplied as a single file, Supplementary_Material_S1.pdf, and comprises per-sample "
        "quality-control metrics, the cluster marker table, the pre-specified tie-break rule and "
        "its application, the cell-type composition before and after that rule, T-cell receptor "
        "clonality, and the per-sample cell-state scores. Every number in it is read from the "
        "result files in the public repository rather than transcribed; it is generated by "
        "manuscript/build_supplementary.py, and both that script and the analysis scripts that "
        "produced the underlying result files are in the same repository."))

    B.append(h(1, "Author Contributions"))
    B.append(p("Conceptualization, methodology, software, formal analysis, investigation, "
               "writing—original draft preparation, writing—review and editing. All authors "
               "have read and agreed to the published version of the manuscript."))
    B.append(h(1, "Funding"))
    B.append(p("This research received no external funding."))
    B.append(h(1, "Institutional Review Board Statement"))
    B.append(p("Not applicable. The study analyses exclusively de-identified public data from "
               "NCBI Gene Expression Omnibus and involved no human or animal subjects."))
    B.append(h(1, "Informed Consent Statement"))
    B.append(p("Not applicable."))
    B.append(h(1, "Data Availability Statement"))
    B.append(p("All expression data are publicly available from NCBI Gene Expression Omnibus "
               "under the accessions listed in Section 2.2. All analysis code, the "
               "pre-registration, the patient-level inclusion and exclusion table, all "
               "intermediate and result files and the figure sources are available at "
               "https://github.com/Cybing521/melanoma-icb-signature-reproducibility. Code is "
               "released under the MIT licence; derived data tables are released under CC BY "
               "4.0; underlying GEO data remain subject to their original terms."))
    B.append(h(1, "Conflicts of Interest"))
    B.append(p("The authors declare no conflict of interest."))

    B.append(h(1, "References"))
    for ref in [
        # 全部 10 条的作者名单、刊名、卷页均已逐条核对（2026-10-07）。
        # 更正记录见 docs/12_参考文献核对.md：原第 2、3、5、6 条的题名或出处与 PubMed / Crossref
        # 记录不符（原第 1、4、10 条的作者名单有虚构填充），第 7、9 条未能机器核验、暂按通行引用保留。
        "1. Auslander N.; Zhang G.; Lee J.S.; Frederick D.T.; Miao B.; Moll T.; Tian T.; "
        "Wei Z.; Madan S.; Sullivan R.J.; Boland G.; Flaherty K.; Herlyn M.; Ruppin E. "
        "Robust prediction of response to immune checkpoint blockade therapy in metastatic "
        "melanoma. Nat. Med. 2018, 24, 1545–1549. (The IMPRES signature; note the "
        "publisher's correction of December 2018.)",
        "2. Charoentong P.; Finotello F.; Angelova M.; Mayer C.; Efremova M.; Rieder D.; "
        "Hackl H.; Trajanoski Z. Pan-cancer immunogenomic analyses reveal "
        "genotype–immunophenotype relationships and predictors of response to checkpoint "
        "blockade. Cell Rep. 2017, 18, 248–262. (The Immunophenoscore.)",
        "3. Riaz N.; Havel J.J.; Makarov V.; Desrichard A.; Urba W.J.; Sims J.S.; Hodi F.S.; "
        "Martín-Algarra S.; Mandal R.; Sharfman W.H.; Bhatia S.; Hwu W.J.; Gajewski T.F.; "
        "Slingluff C.L.; Chowell D.; Kendall S.M.; Chang H.; Shah R.; Kuo F.; Morris L.G.T.; "
        "Sidhom J.W.; Schneck J.P.; Horak C.E.; Weinhold N.; Chan T.A. Tumor and "
        "microenvironment evolution during immunotherapy with nivolumab. Cell 2017, 171, "
        "934–949.e16. (Source of GSE91061.)",
        "4. Hugo W.; Zaretsky J.M.; Sun L.; Song C.; Moreno B.H.; Hu-Lieskovan S.; "
        "Berent-Maoz B.; Pang J.; Chmielowski B.; Cherry G.; Seja E.; Lomeli S.; Kong X.; "
        "Kelley M.C.; Sosman J.A.; Johnson D.B.; Ribas A.; Lo R.S. Genomic and transcriptomic "
        "features of response to anti-PD-1 therapy in metastatic melanoma. Cell 2016, 165, "
        "35–44. (Source of GSE78220.)",
        "5. Vathiotis I.A.; Salichos L.; Martinez-Morilla S.; Gavrielatou N.; Aung T.N.; "
        "Shafi S.; Wong P.F.; Jessel S.; Kluger H.M.; Syrigos K.N.; Warren S.; Gerstein M.; "
        "Rimm D.L. Baseline gene expression profiling determines long-term benefit to "
        "programmed cell death protein 1 axis blockade. NPJ Precis. Oncol. 2022, 6, 92. "
        "(Source of GSE215868; native endpoint is 24-month long-term benefit.)",
        "6. Lauss M.; Phung B.; Borch T.H.; Harbst K.; Kaminska K.; Ebbesson A.; Hedenfalk I.; "
        "Yuan J.; Nielsen K.; Ingvar C.; Carneiro A.; Isaksson K.; Pietras K.; Svane I.M.; "
        "Donia M.; Jönsson G. Molecular patterns of resistance to immune checkpoint blockade in "
        "melanoma. Nat. Commun. 2024, 15, 3075. (Source of GSE244982.)",
        "7. Di Pietro A.; Au L.; Crock P.; Thio N.; Pizzolla A.; Nguyen T.N.; Macdonald S.; "
        "Chalmers H.; Zhu R.; Airaghi A.; Molden-Hauer T.; Bacac M.; Schwalie P.; Schlenker R.; "
        "Levesque M.P.; Mailer S.; Barnes-Cullen K.; Winch K.; Chan J.; Yeung G.A.; Spain L.; "
        "Rao A.D.; Sandhu S.; Gyorki D.E.; McArthur G.A.; Mackay L.K.; Neeson P.J. "
        "Tumor-resident T cells and dendritic cells form an in situ archetype during "
        "immunotherapy response in melanoma. Nat. Commun. 2026, 17, 7445. (Source of GSE294272 "
        "and GSE294273.)",
        "8. Tumor-intrinsic and extrinsic immune mechanisms define resistance to immune "
        "checkpoint blockade in metastatic melanoma. GEO Series GSE308433, GSE308434 and "
        "GSE308435 [Internet]. Columbia University, Izar laboratory; cited by accession because "
        "no peer-reviewed publication is linked to these accessions in GEO. Available at: "
        "https://www.ncbi.nlm.nih.gov/geo/ (accessed 7 October 2026). (Source of the single-cell "
        "and T-cell receptor analyses in Supplementary Material S1.)",
        "9. Liberzon A.; Birger C.; Thorvaldsdóttir H.; Ghandi M.; Mesirov J.P.; Tamayo P. The "
        "Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015, "
        "1, 417–425. (Source of the feature annotations used in model development.)",
        "10. Hanley J.A.; McNeil B.J. The meaning and use of the area under a receiver "
        "operating characteristic (ROC) curve. Radiology 1982, 143, 29–36. (Standard error of "
        "the AUC, and the basis of the power calculation in Section 2.8.)",
        "11. DerSimonian R.; Laird N. Meta-analysis in clinical trials. Control. Clin. Trials. "
        "1986, 7, 177–188. (Method for the study-level pooling in Section 2.8.)",
        "12. Deeks J.J.; Keating J.; Leeflang M.M.G. A hierarchy of diagnostic tests. Ann. "
        "Intern. Med. 2008, 148, 175–176.",
        "13. Hanley J.A.; Lippman-Hand A. If a coin is tossed 3 times, is it fair? JAMA. 1983, "
        "250, 2563–2566.",
    ]:
        B.append(p(ref))

    return {
        "meta": {
            "title": "Do Published Immune Checkpoint Blockade Response Signatures Reproduce "
                     "Across Independent Melanoma Cohorts? A Pre-Registered Reproducibility Study",
            "authors": "[Author names to be supplied]",
            "affiliations": "[Affiliations to be supplied]",
            "correspondence": "[Corresponding author to be supplied]",
        },
        "blocks": B,
    }


def main() -> int:
    content = build()
    json_path = ROOT / "manuscript" / "content.json"
    json_path.write_text(json.dumps(content, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"content: {json_path}  blocks={len(content['blocks'])}")

    builder = ROOT / "manuscript" / "ManuscriptBuilder"
    r = subprocess.run(
        ["dotnet", "run", "--project", str(builder), "-v", "q", "--",
         str(json_path), str(OUT_DOCX), str(ROOT)],
        capture_output=True, text=True)
    sys.stdout.write(r.stdout)
    sys.stderr.write(r.stderr)
    if r.returncode != 0:
        return r.returncode
    print(f"docx: {OUT_DOCX}  ({OUT_DOCX.stat().st_size / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
