# GSE215868 解析自检

生成脚本 `scripts/17_build_gse215868.py`。

## 一、元数据与 RPT 的一致性

- series matrix 样本数：105
- RPT 文件数：105
- 只有元数据无 RPT：无
- 只有 RPT 无元数据：无

## 二、表达矩阵

- 基因数（内源+内参）：770
- 列数（GSM）：105
- 写出：`expr/GSE215868_symbol_matrix.tsv.gz`

## 三、患者级表

- 写出：`meta/GSE215868_samples.tsv`，105 行 × 14 列
- 分析单位：**患者**。105 个 GSM 与 105 个唯一临床组合一一对应（门禁已核）

## 四、关键字段分布（实测，非抄预登记）

- **bor**：
    - `PD`：34
    - `SD`：25
    - `PR`：24
    - `CR`：21
    - `UK`：1
- **response_group**：
    - `PRCR`：45
    - `PD`：34
    - `SD`：25
    - `UK`：1
- **itx_raw**：
    - `IPI+NIVO`：47
    - `PEMBRO`：31
    - `NIVO`：25
    - `NIVO+EXPERIMENTAL`：1
    - `Pembro`：1
- **itx**：
    - `IPI+NIVO`：47
    - `PEMBRO`：32
    - `NIVO`：25
    - `NIVO+EXPERIMENTAL`：1
- **prior_icb**：
    - `NO`：83
    - `YES`：22
- **os_index**：
    - `0`：62
    - `1`：43

- 唯一 (age, bor, itx, pfs_days, os_days) 组合数：105 / 105

- 标题去重后唯一值数：105 / 105
- 标题可解析为板位编号的比例：105 / 105
- 板位（`ITx<板>_<lane>`）分布：{'ITx1_1': 12, 'ITx1_2': 12, 'ITx1_4': 12, 'ITx1_5': 12, 'ITx2_1': 12, 'ITx2_2': 12, 'ITx2_3': 12, 'ITx1_3': 11, 'ITx2_4': 10}
- 孔号重复情况：{'1': 9, '2': 8, '3': 9, '4': 9, '5': 9, '6': 9, '7': 9, '8': 9, '9': 9, '10': 9, '11': 8, '12': 8}
- **板位编号不是患者编号**：孔号在板内唯一，且 105 个样本的 (age, bor, itx, pfs_days, os_days) 组合两两互不相同 → 一 GSM 一患者
