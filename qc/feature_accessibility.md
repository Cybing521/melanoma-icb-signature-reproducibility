# 候选特征的临床可及性分级（任务 11）

生成脚本 `scripts/22_feature_accessibility.py`。
**经验代理**：NanoString PanCancer IO 360 面板——真实临床研究在用的 IO 面板。面板共 784 个条目，其中 6 个阳性对照 + 8 个阴性对照不参与打分，**实测基因 770 个**，本项目已下载其清单（`scripts/17`）。

分级规则见脚本头，**定义性判据、与结果无关**——本项目主线结论为阴性，没有把等级调高的动机。

## 一、分布

| 等级 | 特征数 |
|---|---:|
| C 难可及 | 42 |
| B 条件可及 | 9 |
| D 非常规可及 | 2 |
| A 常规可及 | 1 |

## 二、逐特征

| 特征 | 类型 | 基因数 | 面板覆盖 | 等级 | 测量手段 |
|---|---|---:|---:|---|---|
| `IMPRES_15pairs` | 基因对打分 | 15 | 86.7% | **A 常规可及** | 靶向转录组面板 / 甚至 RT-PCR 即可 |
| `allograft_rejection` | Hallmark 签名 | 200 | 55.0% | **B 条件可及** | bulk RNA-seq / 靶向转录组面板 |
| `il6_jak_stat3_signaling` | Hallmark 签名 | 87 | 54.0% | **B 条件可及** | bulk RNA-seq / 靶向转录组面板 |
| `interferon_alpha_response` | Hallmark 签名 | 97 | 43.3% | **B 条件可及** | bulk RNA-seq / 靶向转录组面板 |
| `interferon_gamma_response` | Hallmark 签名 | 200 | 43.0% | **B 条件可及** | bulk RNA-seq / 靶向转录组面板 |
| `inflammatory_response` | Hallmark 签名 | 200 | 40.0% | **B 条件可及** | bulk RNA-seq / 靶向转录组面板 |
| `angiogenesis` | Hallmark 签名 | 36 | 38.9% | **B 条件可及** | bulk RNA-seq / 靶向转录组面板 |
| `tnfa_signaling_via_nfkb` | Hallmark 签名 | 200 | 36.0% | **B 条件可及** | bulk RNA-seq / 靶向转录组面板 |
| `wnt_beta_catenin_signaling` | Hallmark 签名 | 42 | 35.7% | **B 条件可及** | bulk RNA-seq / 靶向转录组面板 |
| `notch_signaling` | Hallmark 签名 | 32 | 34.4% | **B 条件可及** | bulk RNA-seq / 靶向转录组面板 |
| `il2_stat5_signaling` | Hallmark 签名 | 199 | 28.1% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `apoptosis` | Hallmark 签名 | 161 | 26.7% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `epithelial_mesenchymal_transition` | Hallmark 签名 | 200 | 26.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `pi3k_akt_mtor_signaling` | Hallmark 签名 | 105 | 20.9% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `hedgehog_signaling` | Hallmark 签名 | 36 | 19.4% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `kras_signaling_up` | Hallmark 签名 | 200 | 19.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `tgf_beta_signaling` | Hallmark 签名 | 54 | 18.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `hypoxia` | Hallmark 签名 | 200 | 18.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `complement` | Hallmark 签名 | 200 | 18.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `apical_junction` | Hallmark 签名 | 200 | 15.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `glycolysis` | Hallmark 签名 | 200 | 14.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `p53_pathway` | Hallmark 签名 | 200 | 14.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `coagulation` | Hallmark 签名 | 138 | 13.8% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `apical_surface` | Hallmark 签名 | 44 | 13.6% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `e2f_targets` | Hallmark 签名 | 200 | 13.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `uv_response_up` | Hallmark 签名 | 158 | 13.3% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `uv_response_dn` | Hallmark 签名 | 144 | 13.2% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `mtorc1_signaling` | Hallmark 签名 | 200 | 12.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `estrogen_response_late` | Hallmark 签名 | 200 | 11.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `cholesterol_homeostasis` | Hallmark 签名 | 74 | 10.8% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `g2m_checkpoint` | Hallmark 签名 | 200 | 9.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `xenobiotic_metabolism` | Hallmark 签名 | 200 | 8.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `kras_signaling_dn` | Hallmark 签名 | 200 | 8.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `estrogen_response_early` | Hallmark 签名 | 200 | 7.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `unfolded_protein_response` | Hallmark 签名 | 113 | 7.1% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `androgen_response` | Hallmark 签名 | 101 | 6.9% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `dna_repair` | Hallmark 签名 | 150 | 6.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `myogenesis` | Hallmark 签名 | 200 | 6.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `mitotic_spindle` | Hallmark 签名 | 199 | 5.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `adipogenesis` | Hallmark 签名 | 200 | 5.5% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `myc_targets_v2` | Hallmark 签名 | 58 | 5.2% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `heme_metabolism` | Hallmark 签名 | 200 | 5.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `pancreas_beta_cells` | Hallmark 签名 | 40 | 5.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `fatty_acid_metabolism` | Hallmark 签名 | 158 | 4.4% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `reactive_oxygen_species_pathway` | Hallmark 签名 | 49 | 4.1% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `myc_targets_v1` | Hallmark 签名 | 200 | 4.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `peroxisome` | Hallmark 签名 | 104 | 3.9% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `spermatogenesis` | Hallmark 签名 | 135 | 3.7% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `protein_secretion` | Hallmark 签名 | 96 | 3.1% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `oxidative_phosphorylation` | Hallmark 签名 | 200 | 3.0% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `bile_acid_metabolism` | Hallmark 签名 | 112 | 2.7% | **C 难可及** | bulk RNA-seq / 靶向转录组面板 |
| `GSE308433_state_sd` | 非表达类 | 0 | 0.0% | **C 难可及** | 单核 RNA-seq 状态分散度 |
| `GSE308433_cn_dosage` | 非表达类 | 0 | 0.0% | **D 非常规可及** | 低深度 WGS (lpWGS) 拷贝数 |
| `GSE308435_tcr_clonality` | 非表达类 | 0 | 0.0% | **D 非常规可及** | scTCR 文库测序 |

## 三、对交付物的直接影响

- **A 级 1 个**：常规 FFPE + 商业面板即可拿到，落地无额外门槛
- **B 级 9 个**：需定制面板或更好的 FFPE RNA-seq
- **C 级 42 个**：需全转录组且对降解敏感，**不应写进模型卡的「常规可测」清单**
- **D 级 2 个**：平台层面不可及，**必须从任何声称临床可部署的清单中剔除**

### 一条必须写进结论的观察

**IMPRES 的 15 个基因在本项目所有候选里临床可及性最高**（面板覆盖 86.7%，13/15 个基因在常规 IO 面板上）。

但 `docs/07` 的结论是：**它在 79 例预登记锁定的独立验证上AUC 0.505 / 0.483，与随机无法区分**。

这两条事实并列，正好构成项目最有说服力的一处对比：

> **一个最容易测、测起来最便宜的签名，恰恰是最没有判别力的那个。**
> 它的失败不能归咎于「测不准」「测不到」「样本量不够」——这三项在本项目里都已逐一排除（不变量校验 8/8、面板覆盖 86.7%、n=79）。

**同时必须写明反向的边界**：可及性分级**不能**用来预测判别力。本表给的是测量成本，不是诊断价值；把它读成「C 级特征更值得做」是误读。

