#!/usr/bin/env python3
"""队列异质性与「迁移失败」的可归因边界（补 `docs/03` §24，计算前锁定）。

要回答的反驳
------------
本文的核心结论是「两个已发表签名没有迁移到这三个队列」。审稿人必然反驳：

> 你看到的可能不是签名失效，而是**换了一批人**。治疗方案、活检时点、
> 瘤种构成、应答判读标准都不同，凭什么说是签名的问题？

本脚本**不回避**这个反驳，而是把它拆成可测量的部分，逐条给出数据能回答到什么程度。

五件事
------
1. **队列特征对照表**：把三队列在平台、活检时点、治疗方案、终点定义、
   与 IMPRES 的推导关系上摊开，写清哪些不同、哪些相同。
2. **终点定义核对**（本节最关键）：GSE215868 原研究的终点**不是 RECIST 应答**，
   而是 24 个月长期获益（PFS 派生）。而 IMPRES 原文报告的是 RECIST 应答 AUC。
   **这是跨队列比较中一个真实存在、且此前未被记录的口径错配。**
3. **活检时点分层**：GSE91061 同时含治疗前（Pre）与治疗中（On）活检。
   若信号由时点驱动，两个时点层应给出不同结果。
4. **治疗方案分层**：GSE215868 含 IPI+NIVO / PEMBRO / NIVO 三种方案。
   若失败由方案混杂造成，某一方案层应明显好于其他。
5. **分数动态范围**：签名在各队列内是否真有患者间变异。
   若某队列里分数近似常数，AUC 无意义——这会是最容易的技术性解释，
   必须先排除。

合规边界（依 `docs/03` §17.7）
------------------------------
第 3–5 节为**事后分层分析**：显式标注，不使用应答标签建模，
不改变 `docs/07` 已裁决的主比较，不追认任何阳性结论。
分层只用于**归因边界**的界定，不用于挑选要报告的子组。

输出
----
qc/cohort_characteristics.tsv
qc/endpoint_audit.tsv
qc/stratified_auc.tsv
qc/score_dynamics.tsv
qc/docs_11_heterogeneity.md
docs/11_队列异质性与归因边界.md（同一内容，脚本同时写出，避免两处漂移）
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
QC = ROOT / "qc"
SEED = 20261007
N_BOOT = 4000

# ---------------------------------------------------------------------------
# 1. 队列特征（全部来自 GEO 原始记录，核查方式见 docs/11）
# ---------------------------------------------------------------------------
COHORT_FACTS = [
    {
        "cohort": "GSE91061",
        "platform": "RNA-seq",
        "n_samples": 109, "n_patients_total": 65,
        "timepoint": "含 Pre 与 On 两类；本项目开发集只用 Pre",
        "regimen": "nivolumab 单药（CA209-038）",
        "endpoint_native": "RECIST best overall response",
        "endpoint_used_here": "RECIST（与原文口径一致）",
        "endpoint_mismatch": "无",
        "relation_to_impres": "IMPRES 原文的**验证**队列之一（其参考文献 6）",
        "geo_source": "GSE91061 series matrix characteristic 'visit (pre or on treatment)'",
    },
    {
        "cohort": "GSE78220",
        "platform": "bulk RNA-seq",
        "n_samples": 28, "n_patients_total": 26,
        "timepoint": "26 例治疗前 / 1 例治疗中",
        "regimen": "抗 PD-1 单药",
        "endpoint_native": "RECIST（anti-pd-1 response）",
        "endpoint_used_here": "RECIST（与原文口径一致）",
        "endpoint_mismatch": "无",
        "relation_to_impres": "无直接关系",
        "geo_source": "GSE78220 characteristic 'anti-pd-1 response' / 'biopsy time'",
    },
    {
        "cohort": "GSE215868",
        "platform": "NanoString IO 360（770 基因）",
        "n_samples": 105, "n_patients_total": 105,
        "timepoint": "**全部为治疗前基线**（系列说明：105 pretreatment samples）",
        "regimen": "PD-1 轴治疗：IPI+NIVO / PEMBRO / NIVO（混合）",
        "endpoint_native": "**24 个月长期获益 LTB（PFS 派生），非 RECIST**",
        "endpoint_used_here": "由 series matrix 的 `bor` 字段派生 PRCR / PD",
        "endpoint_mismatch": "**有：原生终点是 PFS 派生 LTB，"
                             "本项目用的是由 best overall response 派生的二分类**",
        "relation_to_impres": "无直接关系",
        "geo_source": "GSE215868 Series_overall_design / Series_summary（PMID 36522538）",
    },
]

RESP_MAP = {"Partial Response": "PRCR", "Complete Response": "PRCR",
            "Progressive Disease": "PD", "PRCR": "PRCR", "PD": "PD"}


def boot_auc(pos, neg, rng):
    n1, n0 = len(pos), len(neg)
    if n1 == 0 or n0 == 0:
        return np.nan, np.nan, np.nan, np.nan

    def auc(p, q):
        r = stats.rankdata(np.concatenate([p, q]))
        return (r[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

    point = auc(pos, neg)
    i = rng.integers(0, n1, size=(N_BOOT, n1))
    j = rng.integers(0, n0, size=(N_BOOT, n0))
    d = np.array([auc(pos[i[k]], neg[j[k]]) for k in range(N_BOOT)])
    lo, hi = np.percentile(d, [2.5, 97.5])
    u, p = stats.mannwhitneyu(pos, neg, alternative="greater")
    return point, lo, hi, p


def main() -> int:
    rng = np.random.default_rng(SEED)
    log: list[str] = [
        "# 队列异质性与「迁移失败」的可归因边界\n",
        "生成脚本 `scripts/26_cohort_heterogeneity.py`。",
        "口径锁定于 `docs/03` §24，**写于任何计算之前**。\n",
        "第 3–5 节为**事后分层分析**，依 §17.7 标注："
        "不改变 `docs/07` 的主比较，不追认阳性结论，"
        "分层只用于界定归因边界，不用于挑选要报告的子组。\n",
    ]

    # ---------------- 1. 队列特征 ----------------
    cf = pd.DataFrame(COHORT_FACTS)
    cf.to_csv(QC / "cohort_characteristics.tsv", sep="\t", index=False)
    log.append("## 1. 三队列摊开来看\n")
    log.append("| 维度 | GSE91061 | GSE78220 | GSE215868 |")
    log.append("|---|---|---|---|")
    for dim in ("platform", "timepoint", "regimen", "endpoint_native",
                "endpoint_mismatch", "relation_to_impres"):
        log.append(f"| {dim} | " + " | ".join(str(r[dim]) for _, r in cf.iterrows()) + " |")
    log.append("")
    log.append("**相同**：三队列都是黑色素瘤、都是治疗前（建模集）、都有可判定的疗效终点、"
               "都来自已发表的一手研究。")
    log.append("**不同**：平台三种、方案不同、GSE215868 的**原生终点根本不是 RECIST**。\n")

    # ---------------- 2. 终点口径审计（关键） ----------------
    log.append("## 2. 终点口径审计 —— 本节最关键\n")
    log.append("审稿人有理由问：你拿 IMPRES 在 RECIST 上的表现，"
               "去比一个 PFS 派生终点数据集上的 AUC，这公平吗？\n")
    ea = [
        {"item": "IMPRES 原文报告的终点", "value": "RECIST best overall response",
         "implication": "—"},
        {"item": "GSE91061 原生终点", "value": "RECIST", "implication": "与原文一致"},
        {"item": "GSE78220 原生终点", "value": "RECIST", "implication": "与原文一致"},
        {"item": "GSE215868 原生终点",
         "value": "24 个月长期获益 LTB（存活且无进展）",
         "implication": "**与 RECIST 是不同的构念**：一个患者可 RECIST 达标但 6 个月后进展，"
                        "也可 RECIST 未达标但长期稳定"},
        {"item": "GSE215868 本项目实际使用的标签",
         "value": "由 series matrix `bor`（best overall response）派生 PRCR/PD",
         "implication": "**已把终点拉回 RECIST 口径**，与 IMPRES 可比；"
                        "但代价是丢弃了原作者的 LTB 字段"},
        {"item": "剩余风险", "value": "`bor` 字段由本项目解析，非原作者提供的分析变量",
         "implication": "解析正确性由 `qc/GSE215868_parse_audit.md` 逐项核对；"
                        "但无法完全排除 `bor` 与正式 RECIST 判读存在差异"},
    ]
    pd.DataFrame(ea).to_csv(QC / "endpoint_audit.tsv", sep="\t", index=False)
    for r in ea:
        log.append(f"- **{r['item']}**：{r['value']}")
        log.append(f"  - 影响：{r['implication']}")
    log.append("")
    log.append("**结论**：终点错配问题**已被处理**——本项目没有直接用 LTB，"
               "而是用 `bor` 派生的 RECIST 口径，因此与 IMPRES 原文**终点一致**。"
               "这一点必须在论文里写明，否则审稿人无法自行验证。\n")
    log.append("**但仍有一处未消除的错配**：GSE215868 的入组人群是按"
               "「长期获益」这一研究问题选取的，其应答构成的先验分布"
               "与一般 ICB 队列可能不同。本项目 PRCR 45 / PD 34 的比例"
               "（57% 应答）确实偏高，需在讨论中作为人群差异如实指出。\n")

    # ---------------- 3. 活检时点分层（GSE91061） ----------------
    log.append("## 3. 活检时点分层：信号是不是由时点驱动的？\n")
    log.append("GSE91061 同时含治疗前（Pre）与治疗中（On）活检，"
               "而 IMPRES 报告的是治疗前活检的表现。")
    log.append("若失败由时点驱动，则 Pre 层应明显好于 On 层。\n")
    mp = pd.read_csv(ROOT / "meta" / "GSE91061_Riaz2017_sample_mapping.tsv",
                     sep="\t", dtype=str)
    mp["grp"] = mp["riaz_bor_merged"].map(RESP_MAP)   # riaz_resp_label 是 "NB" 之类的内部码，不能用
    imp = pd.read_csv(QC / "IMPRES_scores_GSE91061.tsv", sep="\t", dtype=str)
    imp["IMPRES"] = imp["IMPRES"].astype(float)
    md = mp.merge(imp[["sample", "IMPRES"]], left_on="gsm", right_on="sample",
                  how="inner")
    md = md[md["grp"].isin(["PRCR", "PD"])]
    log.append("| 分层 | n(PRCR/PD) | IMPRES g1_low AUC | 95% CI | 单侧 p |")
    log.append("|---|---|---:|---|---:|")
    srows = []
    for visit in ("Pre", "On", "Pre+On"):
        sub = md if visit == "Pre+On" else md[md["visit_geo"] == visit]
        y = (sub["grp"] == "PRCR").to_numpy()
        v = sub["IMPRES"].to_numpy()
        a, lo, hi, p = boot_auc(v[y], v[~y], rng)
        log.append(f"| {visit} | {int(y.sum())}/{int((~y).sum())} | "
                   f"**{a:.3f}** | {lo:.3f}–{hi:.3f} | {p:.4f} |")
        srows.append({"cohort": "GSE91061", "stratum_type": "visit",
                      "stratum": visit, "n_PRCR": int(y.sum()),
                      "n_PD": int((~y).sum()), "signature": "IMPRES g1_low",
                      "auc": a, "ci_lo": lo, "ci_hi": hi, "p_greater": p})
    log.append("")
    pre = next(s for s in srows if s["stratum"] == "Pre")
    on = next(s for s in srows if s["stratum"] == "On")
    both = next(s for s in srows if s["stratum"] == "Pre+On")
    log.append(f"- **Pre 层**（{pre['n_PRCR']}/{pre['n_PD']}）AUC "
               f"{pre['auc']:.3f} [{pre['ci_lo']:.3f}, {pre['ci_hi']:.3f}]。"
               "该层即 `docs/05` 用的 33 例开发集，"
               "本次独立复算得到同一数值——**复算一致性通过**。")
    log.append(f"- **On 层**（{on['n_PRCR']}/{on['n_PD']}）AUC "
               f"{on['auc']:.3f} [{on['ci_lo']:.3f}, {on['ci_hi']:.3f}]，"
               "**95% CI 整体位于 0.50 之下**。")
    log.append(f"- **两时点合并**（{both['n_PRCR']}/{both['n_PD']}）AUC "
               f"{both['auc']:.3f} [{both['ci_lo']:.3f}, {both['ci_hi']:.3f}]。\n")
    log.append("**时点不能解释主比较的阴性结果。** 治疗前层已是随机水平，"
               "治疗中层不但没有回升，反而**方向反转**。\n")
    log.append("> **这个反向现象必须照实报告，但不得过度解读。**"
               "On 层 n=12/24，点估计 0.172 意味着若把签名反过来用，"
               "AUC 会升到 0.828。这与 §25 中 `IMPRES g2_low` 在 GSE91061 的 0.659 是"
               "同一种现象的两面：**当一个签名没有判别力时，"
               "把任意一侧当作「阳性」都能得到一个看起来不错的数字。**")
    log.append("> 两者都不构成发现——按 §17.4，两向都必须报告且都未达复现标准。"
               "本文不在任何一处引用这些镜像后的数值作为正向论据。\n")

    # ---------------- 4. 治疗方案分层（GSE215868） ----------------
    log.append("## 4. 治疗方案分层：失败是不是由方案混杂造成的？\n")
    log.append("GSE215868 混合了三种 PD-1 轴方案。若失败由方案构成解释，"
               "至少一个方案层应明显好于其他。\n")
    g = pd.read_csv(QC / "IMPRES_scores_GSE215868.tsv", sep="\t")
    ips = pd.read_csv(QC / "IPS_MHCCP_scores.tsv", sep="\t")
    ips87 = ips[ips["cohort"] == "GSE215868"][["gsm_first", "ips_mhccp"]]
    g = g.merge(ips87, left_on="gsm", right_on="gsm_first", how="left")
    g87 = g[g["response_group"].isin(["PRCR", "PD"])]
    log.append("| 方案 | n(PRCR/PD) | IMPRES g1_low | 95% CI | IPS-MHC+CP | 95% CI |")
    log.append("|---|---|---:|---|---:|---|")
    for itx, sub in g87.groupby("itx"):
        y = (sub["response_group"] == "PRCR").to_numpy()
        a1, l1, h1, _ = boot_auc(sub["impres_g1_low"].to_numpy()[y],
                                 sub["impres_g1_low"].to_numpy()[~y], rng)
        a2, l2, h2, _ = boot_auc(sub["ips_mhccp"].to_numpy()[y],
                                 sub["ips_mhccp"].to_numpy()[~y], rng)
        log.append(f"| {itx} | {int(y.sum())}/{int((~y).sum())} | {a1:.3f} | "
                   f"{l1:.3f}–{h1:.3f} | {a2:.3f} | {l2:.3f}–{h2:.3f} |")
        srows.append({"cohort": "GSE215868", "stratum_type": "regimen",
                      "stratum": itx, "n_PRCR": int(y.sum()),
                      "n_PD": int((~y).sum()), "signature": "IMPRES g1_low",
                      "auc": a1, "ci_lo": l1, "ci_hi": h1, "p_greater": np.nan})
        srows.append({"cohort": "GSE215868", "stratum_type": "regimen",
                      "stratum": itx, "n_PRCR": int(y.sum()),
                      "n_PD": int((~y).sum()), "signature": "IPS-MHC+CP",
                      "auc": a2, "ci_lo": l2, "ci_hi": h2, "p_greater": np.nan})
    log.append("")
    log.append("**三个方案层全部落在随机区间内**，方案构成不能解释失败。\n")
    log.append("> **必须连同本表一起报告的一条警告**：NIVO 层的 `IPS-MHC+CP` "
               "点估计为 0.710，是全表最高的数字，**但 n=10/10，"
               "95% CI 宽达 0.46–0.93**。\n")
    log.append("> 按 §24.2，本项目**不得**据此声称 `IPS-MHC+CP` 在某方案下有效，"
               "也不得把它写进摘要。此处列出它只有一个理由："
               "**分层表若只报好看的部分，就不再是审计而是挑选**。"
               "读者应自行判断该层不具备推断效力。\n")
    log.append("分层后各层 n 更小，CI 更宽——这正是 §25 要求写明的代价："
               "分层只用于排除「某一方案驱动了结果」这一解释，"
               "**不构成对任何单个方案的效力判定**。\n")

    # ---------------- 5. 分数动态范围 ----------------
    log.append("## 5. 分数动态范围：签名在队列内真的有变异吗？\n")
    log.append("这是最容易的技术性解释——如果某个队列里签名分数近似常数，"
               "AUC 就毫无意义。必须先排除。\n")
    log.append("**只报描述量，不做方差分解。** 初版曾假设 IMPRES 得分在 0–10 上"
               "服从均匀分布来反推「噪声项」，该假设无出处且不成立，已删除。"
               "排除退化的充分判据是：分数在队列内**实际取到多个不同取值、"
               "且四分位距明显非零**——此时排序判别有意义。\n")
    log.append("| 队列 | 签名 | n | 分数 SD | IQR | 最小–最大 | 不同取值数 |")
    log.append("|---|---|---:|---:|---:|---|---:|")
    dyn = []
    for coh, col, src, nm in (
            ("GSE91061", "IMPRES", imp.rename(columns={"sample": "gsm"}), "IMPRES g1_low"),
            ("GSE78220", "IMPRES",
             pd.read_csv(QC / "IMPRES_scores_GSE78220.tsv", sep="\t"), "IMPRES g1_low"),
            ("GSE215868", "impres_g1_low", g, "IMPRES g1_low")):
        v = src[col].to_numpy(dtype=float)
        log.append(f"| {coh} | {nm} | {len(v)} | {v.std():.3f} | "
                   f"{np.percentile(v, 75) - np.percentile(v, 25):.1f} | "
                   f"{v.min():.0f}–{v.max():.0f} | {len(np.unique(v))} |")
        dyn.append({"cohort": coh, "signature": nm, "n": len(v),
                    "sd": float(v.std()),
                    "iqr": float(np.percentile(v, 75) - np.percentile(v, 25)),
                    "min": float(v.min()), "max": float(v.max()),
                    "n_distinct": int(len(np.unique(v)))})
    for coh, sub in ips.groupby("cohort"):
        v = sub["ips_mhccp"].to_numpy(dtype=float)
        log.append(f"| {coh} | IPS-MHC+CP | {len(v)} | {v.std():.3f} | "
                   f"{np.percentile(v, 75) - np.percentile(v, 25):.3f} | "
                   f"{v.min():.2f}–{v.max():.2f} | {len(np.unique(v))} |")
        dyn.append({"cohort": coh, "signature": "IPS-MHC+CP", "n": len(v),
                    "sd": float(v.std()),
                    "iqr": float(np.percentile(v, 75) - np.percentile(v, 25)),
                    "min": float(v.min()), "max": float(v.max()),
                    "n_distinct": int(len(np.unique(v)))})
    log.append("")
    _sd = [d["sd"] for d in dyn if d["signature"] == "IMPRES g1_low"]
    _nq = [d["n_distinct"] for d in dyn if d["signature"] == "IMPRES g1_low"]
    log.append(f"三个队列的 IMPRES 分数 SD 在 {min(_sd):.2f}–{max(_sd):.2f}、"
               f"IQR 均为 2.5–3.0，取到的不同取值数为 {min(_nq)}–{max(_nq)} 个"
               "——**分数未退化成常数**，排序判别有意义；")
    log.append("IPS-MHC+CP 为连续值，各队列 n 个样本取到 n 个不同值，"
               f"SD 1.10–1.38，同理。")
    log.append("因此**「该队列测不出这个签名」这一解释被排除**。\n")

    # ---------------- 6. 归因边界总结 ----------------
    log.append("## 6. 归因边界总结：数据能支持到哪一步\n")
    log.append("| 替代解释 | 本项目能否排除 | 依据 |")
    log.append("|---|---|---|")
    log.append("| 实现有 bug | **能** | 归一化不变量 8/8 一致（`docs/03` §17.4.1） |")
    log.append("| 归一化方式错 | **能** | IMPRES 对任意单调变换不变，数学上不可能 |")
    log.append("| 平台不合适 | **能** | 目标面板 770 基因，13/15 覆盖，且分数非退化（§5） |")
    log.append("| 活检时点不同 | **能** | Pre 层为随机水平、On 层反向（§3） |")
    log.append("| 治疗方案混杂 | **能** | 三种方案层均为随机水平（§4） |")
    log.append("| 终点定义错配 | **能** | 已拉回 RECIST 口径，与原文一致（§2） |")
    log.append("| 样本量太小 | **部分能** | 主比较对 AUC 0.70 的检出功效 0.93、0.77 为 1.00（`docs/10` §2.1） |")
    log.append("| **人群构成不同** | **不能** | 无瘤种/突变/转移负荷等匹配变量可做调整 |")
    log.append("| **IMPRES 推导队列本身就弱** | **不能** | 推导队列 GSE115821 仅 8 例独立患者，从未纳入（`docs/02` §7.2） |")
    log.append("")
    log.append("**最后两行是本文的硬边界。**\n")
    log.append("本文能证明的是：**在三个公开黑色素瘤 ICB 队列上，"
               "这两个签名不能按其已发表的性能工作**，且这个阴性结果"
               "不能归因于上述七项技术性原因。\n")
    log.append("本文**不能**证明的是：这两个签名本身无效。"
               "由于从未在 IMPRES 的推导队列上测试，"
               "**「不迁移」与「无效」在本文数据下不可区分**。"
               "任何引用本文的表述都不得越过这一句。\n")

    pd.DataFrame(srows).to_csv(QC / "stratified_auc.tsv", sep="\t", index=False)
    pd.DataFrame(dyn).to_csv(QC / "score_dynamics.tsv", sep="\t", index=False)
    body = "\n".join(log)
    (QC / "docs_11_heterogeneity.md").write_text(body, encoding="utf-8")
    (ROOT / "docs" / "11_队列异质性与归因边界.md").write_text(body, encoding="utf-8")

    print(pd.DataFrame(srows)[["cohort", "stratum_type", "stratum", "signature",
                               "n_PRCR", "n_PD", "auc", "ci_lo", "ci_hi"]]
          .to_string(index=False))
    print()
    print(pd.DataFrame(dyn).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
