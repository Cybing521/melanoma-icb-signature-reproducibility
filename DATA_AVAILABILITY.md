# Data availability

All primary data are public and were obtained from NCBI Gene Expression Omnibus.
No controlled access is required. Nothing in this repository redistributes the
original GEO files; the scripts fetch them by accession number.

| Accession | URL | What was used |
|---|---|---|
| GSE91061 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061 | `GSE91061_series_matrix.txt.gz` |
| GSE78220 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220 | series matrix + supplementary |
| GSE215868 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215868 | `GSE215868_RAW.tar` (105 NanoString RPT) + series matrix |
| GSE244982 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244982 | supplementary expression table |
| GSE294272 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294272 | supplementary (29 samples) |
| GSE294273 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294273 | scRNA-seq matrices (reference profile) |
| GSE308433 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308433 | lpWGS FACETS output |
| GSE308434 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308434 | per-sample `_sn_counts.csv.gz` |
| GSE308435 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308435 | scTCR libraries |
| GSE115821 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115821 | evaluated, **excluded** (8 patients) |

## External reference data (not redistributed here)

| Resource | URL | Purpose |
|---|---|---|
| MSigDB Hallmark 2024.1.0 | https://www.gsea-msigdb.org/gsea/msigdb/ | signature gene sets |
| NCBI `gene_info` | https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz.gz | Entrez → Symbol mapping (GSE91061 rows are Entrez IDs) |
| IMPRES author code | https://github.com/noamaus/IMPRES-codes | reconstruction of the 15 feature pairs |

## Derived data

Expression matrices rebuilt from GEO are in `expr/`; patient-level metadata and
inclusion/exclusion flags are in `meta/`; all result tables are in `qc/`. These are
derived from CC0/CC-BY public sources and are released under CC BY 4.0.
