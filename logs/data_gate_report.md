# GEO 数据门禁报告

生成时间：服务器执行，核查 6 个队列。

| 队列 | 角色 | 方案登记样本量 | series_matrix 实得 | 状态 | 响应标签字段 | 患者标识字段 |
|---|---|---:|---:|---|---|---|
| GSE308434 | 纵向 ICB 单核转录组 | 42 | 42 | ok | — | — |
| GSE294273 | 淋巴结 scRNA-seq | 12 | 12 | ok | — | — |
| GSE294272 | 淋巴结 bulk RNA-seq | 29 | 29 | ok | — | — |
| GSE78220 | 治疗前 bulk ICB | 28 | 28 | ok | anti-pd-1 response | patient id |
| GSE91061 | ICB bulk（含两瘤种） | 109 | 109 | ok | response | — |
| GSE244982 | 进展后 bulk ICB | 41 | 41 | ok | — | — |

## GSE78220 标签取值
- `anti-pd-1 response`: Complete Response, Partial Response, Progressive Disease

## GSE91061 标签取值
- `response`: PD, PRCR, SD, UNK

