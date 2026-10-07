# IMPRES 基线的获取、实现与核对

> 核对日期 2026-10-07。对应任务 13。对应原始文献：Auslander et al., *Nature Medicine* 2018, PMID 30127394。

## 一、结论摘要

| 队列（治疗前） | 方向约定 | PRCR n | PD n | AUC | 自助法 95% CI |
|---|---|---:|---:|---:|---|
| GSE91061（建模集 33 例） | `g1_low`（主口径） | 10 | 23 | **0.359** | [0.172, 0.552] |
| GSE91061 | `g2_low`（镜像） | 10 | 23 | 0.659 | [0.452, 0.843] |
| GSE78220（验证集 26 例） | `g1_low`（主口径） | 14 | 12 | **0.298** | [0.122, 0.521] |
| GSE78220 | `g2_low`（镜像） | 14 | 12 | 0.583 | [0.351, 0.786] |

**四个 AUC 的置信区间全部包含 0.5。** 无论采用哪个方向约定，IMPRES 在这两个队列的治疗前样本上都**不能显著区分应答与进展**。

**这与原文报告不符**：原文称 IMPRES 在其考察的 ICB 应答队列上 AUC 为 0.77–0.96，且 GSE91061 正是其验证队列之一。本项目在**治疗前子集**上未能复现该性能。可能原因见第五节，**不据此改动算法或换数据**。

## 二、15 个特征对的来源：不凭记忆写

论文正文只给构造方法，具体的 15 对基因在 **Supp. Table 2**。取不到：

| 路径 | 结果 |
|---|---|
| `europepmc.org/articles/PMC6693632?pdf=render` | HTTP 403 |
| Europe PMC `fullTextXML` / `supplementaryFiles` | `Article with id PMC6693632 is not open access one` |
| PMC 网页直下 `NIHMS1032642-supplement-Supp.docx` | 返回 "Preparing to download" 的 JS cookie 门，无直链 |
| `ftp.ncbi.nlm.nih.gov/pub/pmc/oa_package/...` | HTTP 404（该文非 OA） |

**改走作者 GitHub 仓库** `https://github.com/noamaus/IMPRES-codes`（与 GSE115821 核对时同一条替代路径）：

| 文件 | 内容 |
|---|---|
| `ADDITIONAL_CODES/Additional feature sets/FEATS.mat` | 15 个特征索引（Matlab uint16，1-based）：`31 61 128 148 169 237 493 549 567 575 603 606 619 650 710` |
| `ADDITIONAL_CODES/Additional feature sets/CPall.mat` | 28 个 IC 基因（字母序）：BTLA…TNFSF9 |
| `getIMPRESRAT.m` / `classifyImmuneCOMP.m` | 索引 → 基因对的枚举顺序与比值方向 |

### 索引 → 基因对

`C(28,2)=378` 对分两个区块，各自按 `i<j` 升序枚举：

* 区块 1（idx 1–378）：`g1=gene[j], g2=gene[i]`
* 区块 2（idx 379–756）：`g1=gene[i], g2=gene[j]`

还原出的 15 对：

| # | idx | 区块 | g1 | g2 |
|---:|---:|---|---|---|
| 1 | 31 | 1 | CD274 | C10orf54 |
| 2 | 61 | 1 | CD86 | CD200 |
| 3 | 128 | 1 | CD40 | CD274 |
| 4 | 148 | 1 | CD28 | CD276 |
| 5 | 169 | 1 | CD40 | CD28 |
| 6 | 237 | 1 | TNFRSF14 | CD86 |
| 7 | 493 | 2 | CD27 | PDCD1 |
| 8 | 549 | 2 | CD28 | CD86 |
| 9 | 567 | 2 | CD40 | CD80 |
| 10 | 575 | 2 | CD40 | PDCD1 |
| 11 | 603 | 2 | CD80 | TNFSF9 |
| 12 | 606 | 2 | CD86 | HAVCR2 |
| 13 | 619 | 2 | CD86 | TNFSF4 |
| 14 | 650 | 2 | CTLA4 | TNFSF4 |
| 15 | 710 | 2 | PDCD1 | TNFSF4 |

### 映射正确性的独立验证

原文写明候选对必须"至少含一个基因属于 6 个直接 ICB 靶点（CTLA4、CD28、CD80、CD86、PD-1、PD-L1）"。

**还原出的 15 对全部满足这一条**（涉及的靶点基因：CD274、CD86、CD28、CD80、PDCD1、CTLA4）。若索引偏移一位或反了，这个性质会被破坏。这是索引映射正确、不依赖 AUC 高低的独立判据。

## 三、评分实现

论文 Methods：

```
F_{i,j}(x) = 1  若 exp_i(x) < exp_j(x)，否则 0
IMPRES(x)  = Σ_{k=1..15} F(x)        取值 0–15
```

实现（`scripts/14_impres_baseline.py`）：log2(x+1) → 跨样本分位数归一化 → 逐对比较取满足计数。四个队列全部算出，取值分布在 1–12 之间，与 0–15 的定义相容。

**别名处理**：GSE91061 与 GSE244982 用旧名 `C10orf54` 不存在，按 HGNC 官方改名映射到 `VSTM1`；GSE78220 与 GSE294272 直接命中。这是已证实的改名，不是"找不到就随便换一个"。

## 四、方向约定：两个都报，不挑好看的

**作者仓库里的两个 .m 文件本身方向就不一致**：

| 文件 | 区块 1 的比值 | 区块 2 的比值 |
|---|---|---|
| `getIMPRESRAT.m` | `GE(i)/GE(j)` | `GE(j)/GE(i)` |
| `classifyImmuneCOMP.m` | `GE(j)/GE(i)` | `GE(i)/GE(j)` |

论文正文只写 `F_{i,j}=1 当 exp_i < exp_j`，**没说 Supp. Table 2 里每对的书写方向**，而该文件不可获取。

因此脚本默认 `--orientation g1_low`（与 `ratio = exp_{g1}/exp_{g2}` 一致），并提供 `g2_low` 镜像方向。**两个方向都报，不把好看的那个当主结果。** 本项目的主口径固定为 `g1_low`。

结果的方向性事实：两个队列都是 `g1_low < 0.5`、镜像 `> 0.5`。这**只说明二者之一与作者实际用法不符，不能据此反推哪个才对**——判定标准必须来自原文的 Supp. Table 2，而它取不到。

## 五、为什么不显著：三种可能，均未排除

1. **方向约定**（最可能）：`g1_low` 在两个队列都给出反向判别，`g2_low` 给出正向但都不显著。
2. **时点不同**：原文报告的 0.77–0.96 可能包含治疗中样本。本项目按方案第二节的临床决策点，只用**治疗前**样本；GSE91061 治疗前只有 33 例（PRCR 10）。
3. **归一化与量纲**：原文用 quantile-normalized 表达。本项目各队列量纲不同（raw 计数 / FPKM / 作者标准化值），统一 log2 后再做分位数归一化，可能与作者原始处理有差异。

**处置**：按预登记，不因不显著而更换算法、换特征或换队列。IMPRES 作为"已发表 SOTA 预测器"的基线地位不变，但**在本项目的治疗前口径下未能复现原文性能**这一事实必须写进结论，不能只报自己模型的对比数字。

## 六、对模块 4 的影响

模块 4 原计划"说明新模型是否比 IMPRES 增加信息"。现在这个对比的前提变了：

- IMPRES 在同一套治疗前样本上 AUC ≈ 0.30–0.36（主口径）；
- 任何新模型若 AUC 明显高于此，**只能说明优于本项目实现的 IMPRES**，不能写成"优于已发表的 IMPRES"；
- 引用时必须同时给出"原文报告 0.77–0.96"与"本项目治疗前口径下未复现"两条，不得只取其一。

## 七、产物

| 文件 | 内容 |
|---|---|
| `scripts/impres_pairs.py` | 索引 → 基因对的还原，含自校验 |
| `scripts/14_impres_baseline.py` | IMPRES 计算，含 `--orientation` 与 `--no-qnorm` 敏感性开关 |
| `qc/IMPRES_scores_<cohort>.tsv` | 四个队列的主口径分数 |
| `qc/IMPRES_scores_<cohort>_mirror.tsv` | GSE91061 / GSE78220 的镜像方向分数 |
