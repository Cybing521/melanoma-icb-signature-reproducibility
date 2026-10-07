# Supplementary Material S1

## Single-cell and T-cell receptor analyses

This supplement reports two analyses performed alongside the reproducibility study in the main text. They address a different question — mechanisms of resistance and the cellular composition of the microenvironment — and they are reported separately because they do not bear on whether the two published signatures reproduce. Every number below is read directly from the result files in the public repository; the scripts that produced those files are named in each subsection.

The two components differ in scale and must not be described together as one longitudinal cohort. The single-cell component covers 12 patients with one sample each. The T-cell receptor component covers 42 libraries drawn from 34 patients, of whom only 7 contributed more than one biopsy; it therefore supports a cross-sectional comparison with a small paired subset, and no within-patient trajectory analysis is reported.

## S1.1 Single-cell component: data and quality control

Source series GSE294273 (Nat. Commun. 2026, 17, 7445; main text reference 7). Script: `scripts/05`–`scripts/07`, `scripts/13`, `scripts/21`, `scripts/22`.

Of 88,715 nuclei passing initial loading, 85,167 (96.0%) survived quality control, resolving into 41 clusters. Per-sample quality-control metrics are in `qc/scrna_qc_per_sample.tsv` and are reproduced below.

| Sample | Patient | Treatment | Nuclei after QC | Median genes | Median UMI | Median MT% | MT% p95 |
|---|---|---|---|---|---|---|---|
| GSM8901614 | HM002 | Untreated | 9,284 | 1346.0 | 2549.5 | 4.12 | 14.27 |
| GSM8901615 | HM004 | ICI-resistant | 7,475 | 1762.0 | 4474.0 | 3.08 | 8.83 |
| GSM8901616 | HM005 | Untreated | 8,919 | 1581.0 | 3829.0 | 2.92 | 9.71 |
| GSM8901617 | HM006 | ICI-resistant | 10,253 | 2003.0 | 5956.0 | 2.58 | 7.68 |
| GSM8901618 | HM009 | ICI-responder | 6,502 | 1792.0 | 4211.0 | 2.43 | 11.19 |
| GSM8901619 | HM012 | Untreated | 8,714 | 1631.0 | 3810.5 | 3.26 | 15.01 |
| GSM8901620 | HM014 | Untreated | 6,323 | 1659.0 | 4117.0 | 3.53 | 9.25 |
| GSM8901621 | HM019 | ICI-resistant | 1,749 | 1929.0 | 5687.0 | 2.84 | 9.85 |
| GSM8901622 | HM026 | ICI-responder | 3,960 | 1760.0 | 4581.0 | 2.77 | 8.0 |
| GSM8901623 | HM028 | ICI-resistant | 10,376 | 2535.0 | 7894.5 | 1.47 | 6.04 |
| GSM8901624 | HM035 | ICI-responder | 6,620 | 1963.0 | 5371.5 | 1.9 | 7.88 |
| GSM8901625 | HM038 | ICI-responder | 4,992 | 1700.5 | 5209.5 | 2.69 | 8.97 |

## S1.2 Cluster markers

Top-30 differentially expressed genes per cluster by Wilcoxon rank-sum test, from `qc/scrna_cluster_markers.tsv` (1,230 rows). The top 5 genes of each cluster are listed; the complete table is in the repository.

| Cluster | Top 5 markers |
|---|---|
| 0 | CTLA4, IL32, TIGIT, B2M, CD2 |
| 1 | RPL28, RPS12, RPL10, RPL11, RPL13 |
| 2 | JUNB, CD69, ZFP36L2, LTB, TXNIP |
| 3 | IL7R, SRGN, SARAF, CCR7, ANXA1 |
| 4 | RPL34, RPS3A, RPS27A, RPLP2, RPS14 |
| 5 | MALAT1, MT-CO1, MT-CO2, PTPRC, TNFAIP3 |
| 6 | MALAT1, FYN, AC016831.7, IL7R, HSPH1 |
| 7 | NKG7, CCL5, CST7, GZMA, LAG3 |
| 8 | CD74, HLA-DPB1, HLA-DPA1, LSP1, HLA-DRA |
| 9 | CXCL13, CD7, CCL5, HLA-C, UBC |
| 10 | CCL5, NKG7, CD8A, CST7, LAG3 |
| 11 | CCL5, GZMK, CST7, DUSP2, ZFP36L2 |
| 12 | SERPINE2, MT2A, APOD, FTL, VIM |
| 13 | THBS2, FKBP10, RCN3, ITGB8, CALU |
| 14 | MT-ATP6, MT2A, MT-CYB, MT-ND4L, MT-ND2 |
| 15 | FTL, FTH1, C1QB, C1QA, C1QC |
| 16 | TYROBP, PSAP, PLAUR, SOD2, FCER1G |
| 17 | IL3RA, IRF4, SLC15A4, IRF8, IRF7 |
| 18 | CD79A, CD74, HLA-DRA, HLA-DPA1, EZR |
| 19 | MS4A1, CD83, MALAT1, BANK1, AC016831.7 |
| 20 | GAPDH, VIM, ACTB, PPIA, CD63 |
| 21 | MS4A1, CD79A, BANK1, CD83, CD37 |
| 22 | SARAF, RPL13, RPS27A, RPS26, RPS12 |
| 23 | CTSW, GNLY, AREG, KLRD1, TYROBP |
| 24 | HMGB2, MKI67, STMN1, TUBA1B, TMPO |
| 25 | POU2AF1, CD79A, MZB1, SEL1L3, UBE2J1 |
| 26 | CD79A, MS4A1, CD74, BANK1, HLA-DRA |
| 27 | CD79A, CD37, CD74, MS4A1, HLA-DRA |
| 28 | IGKV3-20, IGHV4-34, CD79A, CD37, MS4A1 |
| 29 | IGFBP7, RAMP2, PECAM1, CRIP2, VWF |
| 30 | TRIM33, MDM4, EDNRB, MBNL2, EPS8 |
| 31 | MALAT1, PDE3B, PTPRC, XIST, BACH2 |
| 32 | CD8B, RPS3A, RPL34, RPS12, EEF1A1 |
| 33 | CA8, MCAM, L1CAM, FN1, PMEL |
| 34 | IGFBP2, AKAP12, POSTN, AEBP1, STXBP1 |
| 35 | SEMA5A, GDF15, MCAM, AEBP1, AKAP12 |
| 36 | MGP, IFITM3, HSPA1A, SNRPE, CDK4 |
| 37 | A2M, SFRP1, COL11A1, S100A1, ISYNA1 |
| 38 | GAPDH, S100A1, RPL8, LGALS1, IFI27 |
| 39 | MFSD12, PMEL, PAEP, LY6K, GNG11 |
| 40 | MGP, BANCR, PAGE5, SNRPE, BCAS2 |

## S1.3 The pre-specified tie-break rule and its application

The tie-break rule was written into the pre-registration before any tie-break result existed. Clustering and the marker table are reused unchanged from the locked run; the resolution was not revisited. For each cluster the top-30 markers are tested for hypergeometric enrichment within each candidate reference profile, using the full set of detected genes as the population and α = 0.05; a cluster that fails to reach significance for any candidate is labelled `Unassigned` rather than being forced into the nearest profile.

Applied to 41 clusters, the rule changed the assignment of **11** clusters and left **6** clusters unassigned, together **15.1%** of all nuclei.

### S1.3.1 Clusters whose assignment changed

| Cluster | Original | New | Best candidate | Markers hit | Hypergeometric p |
|---|---|---|---|---|---|
| 4 | T_CD4 | Unassigned | T_CD4 | 0/30 | 1.00e+00 |
| 5 | T_CD8 | Unassigned | T_CD8 | 0/30 | 1.00e+00 |
| 9 | T_CD8 | NK | NK | 2/30 | 1.62e-05 |
| 13 | Melanocyte | Unassigned | Melanocyte | 0/30 | 1.00e+00 |
| 14 | NK | Melanocyte | Melanocyte | 1/30 | 6.25e-03 |
| 16 | Macrophage | DC | DC | 2/30 | 1.62e-05 |
| 24 | T_CD8 | Unassigned | T_CD8 | 0/30 | 1.00e+00 |
| 25 | Plasma | B | B | 4/30 | 1.82e-11 |
| 34 | Melanocyte | Unassigned | Melanocyte | 0/30 | 1.00e+00 |
| 37 | Melanocyte | Fibroblast | Fibroblast | 1/30 | 6.25e-03 |
| 38 | Melanocyte | Unassigned | Melanocyte | 0/30 | 1.00e+00 |

### S1.3.2 The unassigned clusters

Abstention here does not mean the cluster lacks biological meaning; it means the marker set available does not reach it. This is the only evidence a reader needs to judge whether the marker set should be extended, and it is given deliberately.

| Cluster | Previously assigned | Top 5 markers at stake |
|---|---|---|
| 4 | T_CD4 | RPL34, RPS3A, RPS27A, RPLP2, RPS14 |
| 5 | T_CD8 | MALAT1, MT-CO1, MT-CO2, PTPRC, TNFAIP3 |
| 13 | Melanocyte | THBS2, FKBP10, RCN3, ITGB8, CALU |
| 24 | T_CD8 | HMGB2, MKI67, STMN1, TUBA1B, TMPO |
| 34 | Melanocyte | IGFBP2, AKAP12, POSTN, AEBP1, STXBP1 |
| 38 | Melanocyte | GAPDH, S100A1, RPL8, LGALS1, IFI27 |

One implementation correction is recorded here rather than repaired silently: the first version of the tie-break script set the hypergeometric population to the genes appearing within the top-30 lists rather than to the full set of detected genes. That population had itself been enriched by differential-expression filtering, which artificially suppressed enrichment and caused roughly forty per cent of cells to abstain. The population was restored to the full detected gene set. This was a bug fix, not a change of rule made after seeing results: the rule itself was not altered.

## S1.4 Cell-type composition before and after the tie-break

From `qc/scrna_reference_profile_v2.tsv`. 'Before' is the first-pass assignment, 'after' the pre-specified tie-break.

| Cell type | Before | After | Change |
|---|---|---|---|
| T_CD4 | 26.1% | 21.5% | -4.6 pp |
| Melanocyte | 25.7% | 18.9% | -6.8 pp |
| Unassigned | 0.0% | 15.1% | +15.1 pp |
| T_CD8 | 20.9% | 13.2% | -7.7 pp |
| B | 11.4% | 12.3% | +1.0 pp |
| T_reg | 6.0% | 6.0% | +0.0 pp |
| DC | 0.3% | 4.3% | +4.0 pp |
| Fibroblast | 0.0% | 3.8% | +3.8 pp |
| NK | 2.3% | 2.4% | +0.1 pp |
| Macrophage | 5.4% | 1.4% | -4.0 pp |
| pDC | 0.8% | 0.8% | +0.0 pp |
| Endothelial | 0.2% | 0.2% | +0.0 pp |
| Plasma | 1.0% | 0.0% | -1.0 pp |

## S1.5 T-cell receptor clonality

Source series GSE308435 (main text reference 8). Script: `scripts/08`. 42 libraries, all reported; none was excluded.

Median clonality is 0.031 (range 0.000–0.281). Per-library values are in `qc/GSE308435_tcr_clonality.tsv`.

The series is described upstream as 42 sequential biopsies from 34 patients. That count is a tally of libraries, not a paired design: only 7 of the 34 patients contributed more than one biopsy, and a single untreated anchor exists across the whole cohort. Any statement about within-patient change therefore rests on those few patients and none is made here.

## S1.6 Cell-state scores

Per-sample dispersion of five immune state programmes across 42 samples, from `qc/GSE308434_state_scores.tsv`. Each programme is summarised by its standard deviation across nuclei, its 90th-minus-10th percentile range, and the percentage of nuclei in the upper decile.

| Sample | Nuclei after QC | QC kept % | antigen_presentation | ifng_response | proliferation | emt | hypoxia | tnf_nfkb | ifna_response | inflammatory | tgf_beta | il2_stat5 | glycolysis | oxphos | mtorc1 | pi3k_akt_mtor | apoptosis | angiogenesis | p53 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSM9245398_F01_on | 11,053 | 98.11 | 8.4 | 12.3 | 11.1 | 8.6 | 25.2 | 17.2 | 16.8 | 12.5 | 28.2 | 18.8 | 27.0 | 26.4 | 25.7 | 25.8 | 24.0 | 15.4 | 25.6 |
| GSM9245399_F01_pre | 8,258 | 97.6 | 26.8 | 23.1 | 12.3 | 16.6 | 29.4 | 22.6 | 24.0 | 19.9 | 31.1 | 26.0 | 29.0 | 26.7 | 29.0 | 29.7 | 29.5 | 21.9 | 29.4 |
| GSM9245400_F02_on | 5,954 | 97.88 | 11.5 | 20.1 | 17.2 | 12.1 | 26.5 | 22.0 | 22.9 | 15.1 | 31.9 | 28.8 | 30.8 | 29.8 | 30.9 | 34.4 | 28.8 | 13.0 | 29.5 |
| GSM9245401_F02_pre | 10,095 | 96.88 | 8.2 | 16.3 | 15.7 | 10.3 | 24.7 | 20.5 | 21.8 | 11.0 | 29.0 | 23.5 | 28.8 | 26.9 | 29.5 | 30.6 | 26.7 | 12.8 | 29.0 |
| GSM9245402_F03_post1_on2 | 12,658 | 81.51 | 13.4 | 19.5 | 17.2 | 12.6 | 32.3 | 28.6 | 21.6 | 22.1 | 32.9 | 34.3 | 35.0 | 29.1 | 33.2 | 36.0 | 32.4 | 18.0 | 34.4 |
| GSM9245403_F03_post1_pre2 | 9,216 | 99.26 | 13.1 | 19.7 | 13.2 | 9.7 | 26.1 | 21.4 | 23.1 | 16.9 | 27.1 | 26.0 | 29.7 | 29.7 | 29.3 | 30.8 | 26.9 | 17.0 | 29.6 |
| GSM9245404_F04_pre | 10,780 | 97.11 | 15.0 | 21.1 | 27.5 | 9.0 | 27.7 | 20.8 | 22.7 | 17.4 | 28.6 | 27.7 | 29.9 | 27.0 | 29.6 | 31.1 | 28.2 | 16.0 | 28.8 |
| GSM9245405_F05_pre | 11,033 | 97.04 | 18.1 | 20.5 | 26.1 | 13.9 | 29.2 | 24.2 | 22.6 | 20.0 | 31.9 | 30.1 | 32.5 | 29.8 | 32.0 | 33.3 | 31.2 | 17.2 | 30.4 |
| GSM9245406_F06_post1_pre2 | 3,439 | 98.01 | 25.4 | 28.6 | 28.8 | 27.6 | 30.8 | 30.1 | 30.5 | 26.9 | 32.2 | 30.6 | 32.0 | 32.9 | 32.2 | 32.6 | 31.6 | 26.9 | 31.1 |
| GSM9245407_F07_post1_pre2 | 2,630 | 97.81 | 29.4 | 27.9 | 16.6 | 18.4 | 23.7 | 24.8 | 26.7 | 27.1 | 27.9 | 29.5 | 26.1 | 25.7 | 25.6 | 30.7 | 27.1 | 22.2 | 26.6 |
| GSM9245408_F08_post | 13,827 | 100.0 | 16.2 | 23.7 | 16.9 | 14.5 | 31.6 | 31.3 | 24.3 | 20.9 | 32.0 | 31.8 | 32.7 | 31.8 | 32.6 | 32.8 | 32.2 | 19.9 | 32.3 |
| GSM9245409_F09_post | 9,916 | 99.37 | 14.1 | 21.0 | 24.4 | 18.4 | 28.5 | 23.5 | 21.4 | 16.3 | 28.6 | 26.1 | 29.2 | 27.0 | 29.0 | 29.8 | 28.0 | 22.4 | 29.1 |
| GSM9245410_F10_post | 5,082 | 99.69 | 25.7 | 27.3 | 17.4 | 29.2 | 27.8 | 24.8 | 26.5 | 25.3 | 30.4 | 29.1 | 28.1 | 25.8 | 25.0 | 29.7 | 28.6 | 27.3 | 26.5 |
| GSM9245411_F12_post | 10,870 | 98.77 | 13.3 | 24.8 | 24.6 | 25.4 | 30.5 | 23.8 | 25.8 | 25.6 | 34.5 | 33.7 | 33.9 | 30.8 | 32.9 | 36.2 | 33.8 | 24.1 | 34.2 |
| GSM9245412_F12_pre | 5,530 | 96.54 | 28.5 | 26.1 | 19.1 | 18.7 | 26.1 | 24.4 | 24.8 | 24.8 | 25.7 | 31.0 | 30.3 | 27.2 | 26.9 | 34.2 | 28.4 | 21.4 | 30.3 |
| GSM9245413_F15_pre | 7,105 | 98.83 | 19.5 | 21.1 | 13.2 | 12.1 | 29.4 | 24.8 | 22.2 | 20.0 | 29.3 | 27.6 | 31.5 | 32.5 | 32.4 | 32.2 | 30.5 | 17.3 | 32.0 |
| GSM9245414_F16_post1_pre2 | 7,353 | 99.06 | 10.1 | 24.1 | 14.6 | 23.0 | 27.5 | 24.6 | 26.0 | 19.5 | 31.7 | 28.6 | 30.9 | 27.1 | 29.2 | 31.2 | 29.3 | 24.1 | 31.1 |
| GSM9245415_F16_pre | 9,576 | 89.53 | 13.2 | 20.0 | 21.2 | 18.0 | 25.3 | 20.0 | 20.0 | 18.2 | 29.9 | 25.9 | 30.5 | 30.1 | 30.0 | 32.8 | 28.1 | 18.3 | 29.2 |
| GSM9245416_F17_post | 11,682 | 99.98 | 19.1 | 21.8 | 17.2 | 13.3 | 30.6 | 32.8 | 21.2 | 19.3 | 31.6 | 31.3 | 30.3 | 31.1 | 33.6 | 32.7 | 33.3 | 18.2 | 32.6 |
| GSM9245417_F18_post | 3,531 | 100.0 | 31.6 | 33.5 | 34.8 | 33.2 | 35.4 | 34.6 | 33.2 | 32.7 | 33.3 | 34.7 | 34.7 | 34.7 | 35.5 | 34.9 | 34.3 | 28.7 | 34.4 |
| GSM9245418_F20_post1_pre2 | 12,430 | 99.94 | 15.3 | 19.2 | 15.5 | 13.1 | 35.8 | 25.1 | 23.4 | 17.4 | 34.8 | 33.1 | 36.6 | 35.3 | 35.6 | 36.9 | 35.3 | 20.3 | 37.5 |
| GSM9245419_F22_post | 13,824 | 99.4 | 18.8 | 20.9 | 17.8 | 23.8 | 28.7 | 23.4 | 23.2 | 22.7 | 30.9 | 28.7 | 29.7 | 27.5 | 28.2 | 31.9 | 28.4 | 27.5 | 29.2 |
| GSM9245420_F23_post | 689 | 97.59 | 15.8 | 18.7 | 28.4 | 22.5 | 29.2 | 23.1 | 21.8 | 18.3 | 33.4 | 29.5 | 32.9 | 28.7 | 31.5 | 33.5 | 28.0 | 20.8 | 29.8 |
| GSM9245421_F25_post | 4,752 | 98.71 | 29.5 | 28.3 | 27.9 | 11.3 | 22.8 | 20.4 | 28.0 | 22.6 | 24.9 | 27.0 | 27.5 | 29.1 | 30.8 | 32.5 | 27.1 | 17.3 | 27.0 |
| GSM9245422_F26_post | 4,844 | 99.36 | 27.7 | 28.2 | 18.5 | 16.5 | 28.4 | 28.6 | 26.3 | 25.2 | 26.8 | 31.0 | 29.0 | 25.4 | 29.0 | 32.2 | 29.4 | 25.6 | 27.2 |
| GSM9245423_F27_post1_pre2 | 2,096 | 92.78 | 20.3 | 22.4 | 15.8 | 24.6 | 27.4 | 19.2 | 24.2 | 17.0 | 33.1 | 27.7 | 31.2 | 27.5 | 31.2 | 34.4 | 29.2 | 21.6 | 29.8 |
| GSM9245424_F28_post1_pre2 | 11,275 | 100.0 | 11.3 | 21.8 | 17.3 | 9.1 | 31.1 | 26.1 | 24.9 | 17.2 | 31.0 | 28.5 | 30.8 | 29.7 | 30.9 | 32.0 | 30.1 | 19.7 | 30.1 |
| GSM9245425_F29_post1_pre2 | 10,049 | 93.1 | 16.8 | 25.6 | 28.0 | 12.7 | 33.0 | 29.9 | 26.6 | 22.9 | 33.0 | 31.4 | 32.8 | 29.4 | 32.4 | 34.5 | 32.4 | 18.3 | 32.3 |
| GSM9245426_F30_post | 10,063 | 92.54 | 23.8 | 25.4 | 20.1 | 22.4 | 30.7 | 26.3 | 27.0 | 25.7 | 32.8 | 30.9 | 31.9 | 29.0 | 30.7 | 33.6 | 30.3 | 24.7 | 31.9 |
| GSM9245427_F31_post | 18,759 | 99.99 | 25.4 | 31.1 | 27.7 | 28.9 | 32.6 | 30.9 | 31.8 | 29.3 | 32.9 | 33.0 | 32.2 | 32.1 | 32.3 | 33.7 | 32.7 | 27.2 | 31.9 |
| GSM9245428_R204_pre | 16,418 | 99.99 | 14.5 | 24.0 | 13.1 | 14.5 | 31.1 | 28.2 | 26.7 | 22.8 | 35.2 | 31.4 | 32.9 | 33.6 | 31.8 | 34.3 | 33.2 | 16.1 | 32.9 |
| GSM9245429_R294_on | 10,311 | 99.73 | 27.2 | 24.2 | 15.4 | 18.4 | 29.9 | 23.5 | 25.2 | 20.9 | 32.4 | 30.7 | 32.8 | 31.5 | 31.8 | 33.7 | 32.2 | 16.3 | 30.9 |
| GSM9245430_R308_pre | 13,690 | 99.25 | 26.8 | 26.0 | 19.0 | 20.2 | 28.0 | 24.3 | 26.2 | 24.7 | 29.9 | 31.8 | 26.2 | 23.7 | 25.7 | 32.1 | 30.2 | 20.7 | 28.9 |
| GSM9245431_R310_on1 | 11,699 | 100.0 | 29.1 | 21.3 | 13.7 | 9.8 | 28.4 | 17.7 | 23.1 | 17.1 | 31.3 | 27.2 | 33.0 | 28.2 | 31.7 | 34.4 | 28.1 | 17.1 | 30.9 |
| GSM9245432_R310_on2 | 13,870 | 99.99 | 29.4 | 23.1 | 17.7 | 17.5 | 27.1 | 25.1 | 25.9 | 20.3 | 29.3 | 31.0 | 26.5 | 25.7 | 27.7 | 31.9 | 30.0 | 20.6 | 28.4 |
| GSM9245433_R310_pre | 13,406 | 100.0 | 11.7 | 20.2 | 15.4 | 8.8 | 27.2 | 19.1 | 25.6 | 15.9 | 29.5 | 26.9 | 30.6 | 29.9 | 30.9 | 31.4 | 27.3 | 16.2 | 30.1 |
| GSM9245434_R319_on | 13,681 | 98.05 | 15.9 | 21.9 | 15.6 | 11.0 | 25.7 | 20.6 | 23.3 | 17.1 | 27.9 | 25.8 | 29.6 | 28.2 | 29.9 | 30.5 | 28.0 | 17.6 | 29.8 |
| GSM9245435_R319_pre | 8,777 | 99.5 | 20.4 | 22.5 | 18.7 | 15.1 | 25.7 | 20.2 | 23.5 | 19.2 | 28.1 | 25.4 | 27.9 | 26.6 | 26.8 | 28.8 | 27.6 | 19.8 | 28.4 |
| GSM9245436_R328_on | 9,625 | 99.96 | 13.5 | 18.3 | 16.5 | 11.6 | 27.5 | 19.5 | 20.9 | 16.3 | 27.6 | 25.6 | 29.4 | 28.7 | 28.2 | 29.7 | 27.2 | 17.1 | 27.1 |
| GSM9245437_R329_on | 2,699 | 98.68 | 24.6 | 26.8 | 18.8 | 22.7 | 29.9 | 27.6 | 27.4 | 25.6 | 30.6 | 28.7 | 29.2 | 28.7 | 28.4 | 31.4 | 30.0 | 21.2 | 29.9 |
| GSM9245438_R334_pre | 14,954 | 99.46 | 15.8 | 20.3 | 20.4 | 14.9 | 28.7 | 25.6 | 20.8 | 20.2 | 32.8 | 31.4 | 32.1 | 30.6 | 31.8 | 33.8 | 31.1 | 17.3 | 32.0 |
| GSM9245439_R354_pre | 19,224 | 98.49 | 16.3 | 20.4 | 26.8 | 11.1 | 26.4 | 23.4 | 23.1 | 19.8 | 30.6 | 28.8 | 31.4 | 28.3 | 30.3 | 31.8 | 29.4 | 18.9 | 29.7 |

## S1.7 What this supplement does and does not support

Because a substantial fraction of nuclei could not be assigned with confidence, the deconvolution performed here supports coarse statements only — the broad division of the microenvironment into lymphoid, myeloid, melanocytic and stromal compartments. No claim of fine-grained cell-type decomposition is made, and none of these results is used to support or qualify the reproducibility conclusion in the main text.

All result files are in the public repository at the address given in the Data Availability statement, together with the scripts that produced them.
