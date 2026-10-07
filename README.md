# Are published immune-checkpoint-blockade response signatures reproducible? A pre-registered negative study in melanoma

**Short title:** Reproducibility of published ICB-response signatures in melanoma

A pre-registered, multi-cohort reproducibility study of published immune checkpoint
blockade (ICB) response signatures in melanoma, using exclusively public data.

---

## What this repository is

This is the full analysis code, pre-registration, and results for a study that set out
to build an ICB-response classifier for melanoma and instead found that **two published,
independently-derived response signatures do not reproduce their published performance
in independent cohorts.**

The repository is organised so that every claim in the manuscript can be traced to a
script, a pre-registration section written *before* the corresponding computation, and
a result file.

### The central finding

Two signatures were tested, chosen because they use **mathematically different
constructions** and come from **independent research groups**:

| | Signature 1 | Signature 2 |
|---|---|---|
| Name | IMPRES (Auslander et al., *Nat Med* 2018) | `IPS-MHC+CP` (reduced Immunophenoscore) |
| Construction | pairwise ordering of 15 immune-checkpoint genes | weighted z-scores of 20 MHC + checkpoint genes |
| Cohorts | GSE91061, GSE78220, GSE215868 | GSE91061, GSE78220, GSE215868 |
| Best AUC (largest cohort, n=79) | 0.505 [0.380, 0.631] | 0.591 [0.463, 0.714] |
| Direction consistent across cohorts | **no** | **yes** (0.604 / 0.583 / 0.591) |
| Gene list obtainable from open sources | **no** (rebuilt from author GitHub) | **yes** |
| Coverage by a routine clinical panel | 13/15 (86.7%) | 20/20 (100%) |

Neither reproduces its published performance. The two fail in **different ways**, which
is the substantive result: Signature 1 is directionally erratic, Signature 2 gives a
weak but consistently-signed signal that this design can neither confirm nor exclude.

A second finding concerns **accessibility rather than accuracy**: the complete gene list
of the field's second-most-cited ICB-response signature could not be retrieved from any
open source through six documented routes.

### What we explicitly ruled out

We pre-specified, before running anything, that a negative result would not be
attributed to the usual technical excuses. Each was tested and excluded:

| Excuse | How it was excluded | Artefact |
|---|---|---|
| "Your implementation is buggy" | Monotone-transformation invariance check: 4 normalisations × 2 directions produced **8/8 identical** per-sample scores | `qc/GSE215868_impres_invariant_check.tsv` |
| "Your normalisation is wrong" | Same check — IMPRES is provably invariant to any monotone transform | `docs/03` §17.4.1 |
| "The platform is unsuitable" | Feature-informativeness compared across 5 cohorts / 3 platforms; the targeted panel was **less** degenerate than the full-transcriptome cohorts | `qc/GSE215868_panel_diagnostic.tsv` |
| "The sample is too small" | Largest cohort n=79 (4.5× the development set); CI upper bound 0.631 excludes the published range from 0.70 | `docs/07` |

---

## Repository layout

```
docs/      Pre-registration and results, in execution order
             03  Endpoint definitions and analysis pre-registration  ← read first
                 §17  GSE215868 IMPRES validation protocol  (written pre-download)
                 §19  Cluster-assignment tie-break rule   (written pre-result)
                 §20  Feature-set deviation record        (self-reported violation)
                 §22  Second-signature protocol            (written pre-computation)
             07  GSE215868 IMPRES validation results
             08  Pre-registered primary analysis, re-run
             09  Second-signature (IPS-MHC+CP) results
scripts/   Numbered in execution order, 00–24
figures/   Manuscript figures, vector PDF + PNG (see below)
meta/      Sample-level metadata and the patient-level inclusion/exclusion table
qc/        All intermediate and result artefacts (TSV + markdown reports)
expr/      Gene × GSM expression matrices (five cohorts)
01_方案/    Study design document (Chinese)
PROJECT.md Status panel
todolist.md  Task ledger with acceptance criteria
```

---

## Figures

Manuscript figures, drawn to MDPI specifications (≤175 mm width, Arial 8–10 pt,
lowercase bold parenthesised panel labels). Vector PDF is the delivery format; PNG is
provided for preview. 600 dpi TIFFs are **not** committed (≈145 MB) — regenerate with
the figure script.

| Figure | File | Printed width | Content |
|---|---|---|---|
| 1 | `figures/Figure1_StudyDesign.pdf` | 163.6 mm | Cohort design, patient counts, and where each pre-registered decision applies |
| 2 | `figures/Figure2_SignaturePerformance.pdf` | 172.9 mm | AUC of both signatures across all three validation cohorts with 95 % CIs |
| 3 | `figures/Figure3_RuledOutExplanations.pdf` | 151.4 mm | Each pre-specified technical excuse and the artefact that excludes it |
| 4 | `figures/Figure4_NestedCV.pdf` | 158.4 mm | Nested-CV performance under three AUC conventions, 20 repeats |

```bash
python scripts/24_figures_mdpi.py   # writes PDF + PNG + 600 dpi TIFF
```

---

## Data provenance

All expression data are public. Nothing requires controlled access.

| Accession | Role | Platform | Samples → patients |
|---|---|---|---|
| GSE91061 | development | RNA-seq | 109 → 65 (nivolumab monotherapy; Riaz et al. *Cell* 2017) |
| GSE78220 | external validation | microarray | 28 → 26 patients |
| GSE215868 | large-scale validation | NanoString IO 360 | 105 → 105 |
| GSE244982 | direction consistency | RNA-seq | 41 (all post-progression) |
| GSE294272 | reference validation | RNA-seq | 29 → 26 |
| GSE308433/434/435 | paired tumour–TCR analysis | snRNA-seq / scTCR / lpWGS | 42 → 34 |

Fetched directly from NCBI GEO. On macOS/Linux, use `curl --noproxy '*'` for NCBI domains.

---

## Reproducing the analyses

```bash
# 1. metadata and data gate
python scripts/00_fetch_geo_metadata.py

# 2. downloads (expression matrices, TCR, WGS)
python scripts/02_download_supl.py
python scripts/04_download_izar.py

# 3. patient-level inclusion/exclusion table (single source of truth for all n)
python scripts/05_build_patient_level.py

# 4. expression matrices + gene coverage gate
python scripts/06_build_expression_matrices.py
python scripts/11_check_feature_coverage.py

# 5. IMPRES baseline
python scripts/14_impres_baseline.py --cohort GSE91061 --orientation g1_low

# 6. nested CV modelling
python scripts/16_nested_cv_model.py --outdir .
#    the pre-registered primary panel only:
python scripts/16_nested_cv_model.py --outdir . --tag dev_cv_mainpanel \
  --features antigen_presentation,ifng_response,g2m_checkpoint,emt,hypoxia,tnf_nfkb

# 7. the two signature validations
python scripts/17_build_gse215868.py
python scripts/18_impres_gse215868.py
python scripts/23_ips_mhccp.py

# 8. diagnostics
python scripts/19_panel_diagnostic.py
python scripts/20_task14_external_validation.py
python scripts/22_feature_accessibility.py
```

Random seed is **20261007** throughout. Each result file records the seed used.

---

## Honest limitations

These are stated in the manuscript, not buried:

1. **The IMPRES gene pairs are reconstructed, not copied.** Supp. Table 2 of the
   original paper is not open access; four retrieval routes were tried and failed. The
   15 pairs were rebuilt from the authors' GitHub. The reconstruction has an independent
   validity check (all 15 pairs contain ≥1 direct ICB target), but this remains the
   largest attack surface in the study.
2. **`IPS-MHC+CP` is not the full Immunophenoscore.** Two of its four classes
   (`EC`, `SC`) have gene lists that could not be obtained openly.
3. **The development set is 33 patients with 10 responders.** A feature-count cap of ≤10
   was pre-registered for this reason.
4. **One self-reported procedural deviation.** The feature set used for the main model
   differed from the pre-registered set (3 added, 4 removed). It is recorded in
   `docs/03` §20, the pre-registered primary panel was re-run, and both are reported.
   The violation record is retained deliberately.
5. **A weak signal is not a positive finding.** The `IPS-MHC+CP` consistency across
   cohorts is reported as an unconfirmed observation requiring far larger cohorts, not
   as partial replication.

---

## Citation

If you use this repository, please cite the associated manuscript. A `CITATION.cff`
file is provided.

## License

Code: MIT (see `LICENSE`). Derived data tables in `meta/` and `qc/`: CC BY 4.0.
Underlying GEO data remain subject to their original terms.