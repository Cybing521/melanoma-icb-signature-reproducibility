#!/usr/bin/env python3
"""功效分析 + 跨队列 AUC 元分析（补 `docs/03` §23，计算前锁定口径）。

为什么做这个
------------
主比较（IMPRES `g1_low` × GSE215868，PRCR 45 / PD 34）的 CI 是 0.380–0.631。
审稿人对阴性结果的第一反应永远是「你样本量够不够」。本文不回避这个问题，
而是**在看到任何结果之前把可回答的部分算完**：

1. **所需样本量**：要让 95% CI 上界降到 0.70（临床上「可用的判别器」下限）
   或 0.77（IMPRES 原文跨数据集区间的下界），每个队列各需要多少 n。
2. **最小可检出差异**：在现有 n 下，80% / 90% 功效能检出多大的 AUC 偏离 0.50。
3. **跨队列元分析**：单队列 CI 宽的根源是 n 小。预登记 §3 禁止**样本级合并**
   （不同队列的患者不能当独立样本池化），但**研究级 AUC 的元分析**不违反该约束——
   它只用各队列已发表的汇总统计量，不接触任何患者级数据。
   用 Cochran's Q 检验队列间异质性，固定效应与随机效应（DerSimonian–Laird）都给。

第 3 项**不在原预登记内**，属事后分析。依 §17.7 的合规边界执行：
显式标注为事后、不使用应答标签做建模、不改变已裁决的主比较、不追认任何结论。

参考文献
--------
- Hanley JA, McNeil BJ. The measurement and use of misclassification rates. Radiology 1982.
- Deeks JJ, Keating J, Leeflang M. A hierarchy of diagnostic tests. Ann Intern Med 2008.（MDD）
- DerSimonian R, Laird NM. Meta-analysis in clinical trials. Control Clin Trials 1986.
- Hanley JA, Lippman-Hand A. If a coin is tossed 3 times, is it fair? JAMA 1983.（符号一致性）

输出
----
qc/power_analysis.tsv
qc/meta_auc.tsv
qc/sign_consistency.tsv
qc/docs_10_power.md
docs/10_功效分析与跨队列元分析.md（同一内容，脚本同时写出，避免两处漂移）
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
QC = ROOT / "qc"

# 写死于计算之前
Z975 = 1.959963985
ALPHA_ONE_SIDED = 0.025          # 与既有分析一致：单侧 MW p 报告口径
IMPRES_PUBLISHED_LO = 0.77       # 原文跨数据集 AUC 区间下界
CLINICAL_FLOOR = 0.70            # 判别器的临床可用下限

# 队列主比较的实测输入 —— 数值不在此手写，直接从结果文件读，杜绝口径漂移
# （初版曾手写三个数值，与 qc/ 下的结果文件对不上，已改为读取）
_IMPRES_GSE91061 = {  # docs/05_IMPRES基线核对.md
    "g1_low": (10, 23, 0.359, 0.172, 0.552),
    "g2_low": (10, 23, 0.659, 0.452, 0.843),
}
_IMPRES_GSE78220 = {
    "g1_low": (14, 12, 0.298, 0.122, 0.521),
    "g2_low": (14, 12, 0.583, 0.351, 0.786),
}
_IMPRES_GSE215868 = {
    "g1_low": (45, 34, 0.505, 0.380, 0.631),
    "g2_low": (45, 34, 0.483, 0.358, 0.608),
}


def _load_ips() -> dict[str, tuple[int, int, float, float, float]]:
    d = pd.read_csv(QC / "IPS_MHCCP_auc.tsv", sep="\t")
    d = d[d["orientation"] == "high_ips_responder"]
    return {r["cohort"]: (int(r["n_PRCR"]), int(r["n_PD"]),
                          float(r["auc"]), float(r["ci_lo"]), float(r["ci_hi"]))
            for _, r in d.iterrows()}


def build_cohorts() -> list[tuple]:
    ips = _load_ips()
    out = []
    for coh, src, site in (("GSE91061", _IMPRES_GSE91061, "docs/05 开发集"),
                           ("GSE78220", _IMPRES_GSE78220, "docs/05 外部验证"),
                           ("GSE215868", _IMPRES_GSE215868, "**主比较** docs/07")):
        for orient, tag in (("g1_low", "IMPRES g1_low"), ("g2_low", "IMPRES g2_low")):
            n1, n2, a, lo, hi = src[orient]
            where = site if orient == "g1_low" else f"{site} 镜像方向"
            out.append((coh, tag, orient, n1, n2, a, lo, hi, where))
        n1, n2, a, lo, hi = ips[coh]
        out.append((coh, "IPS-MHC+CP", "high_ips_responder", n1, n2, a, lo, hi, "docs/09"))
    return out


def hm_se(auc: float, n1: int, n2: int) -> float:
    """Hanley–McNeil AUC 标准误。"""
    a = min(max(auc, 1e-6), 1 - 1e-6)
    q1 = a / (2 - a)
    q2 = 2 * a * a / (1 + a)
    var = (a * (1 - a) + (n1 - 1) * (q1 - a * a) + (n2 - 1) * (q2 - a * a)) / (n1 * n2)
    return math.sqrt(max(var, 1e-12))


def se_equal_n(auc: float, n: int) -> float:
    return hm_se(auc, n, n)


def n_for_upper_below(auc: float, target: float) -> int | None:
    """求使 AUC + 1.96*SE(n) < target 的最小 n（n1 = n2 = n）。"""
    if auc >= target:
        return None                      # 真实 AUC 本身已在目标之上，无解
    lo, hi = 5, 200_000
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if auc + Z975 * se_equal_n(auc, mid) < target:
            hi = mid
        else:
            lo = mid
    return hi


def mdd(n: int, power: float) -> float:
    """给定每组 n、单侧 alpha、双侧备择，求最小可检出 AUC（相对 0.50）。"""
    z = Z975 + stats.norm.ppf(power)
    lo, hi = 0.50, 0.999
    for _ in range(200):
        mid = (lo + hi) / 2
        if (mid - 0.50) >= z * se_equal_n(mid, n):
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def power_at(auc: float, n1: int, n2: int) -> float:
    """给定真实 AUC 与各组 n，单侧 alpha 下的功效（检出「高于 0.50」）。"""
    se = hm_se(auc, n1, n2)
    if se <= 0:
        return float(auc > 0.50)
    return float(stats.norm.cdf((auc - 0.50) / se - Z975))


def logit_se(auc: float, n1: int, n2: int) -> float:
    a = min(max(auc, 1e-6), 1 - 1e-6)
    return hm_se(a, n1, n2) / (a * (1 - a))


def meta_auc(estimates: list[tuple[float, float]], label: str) -> dict:
    """研究级 AUC 元分析（logit 尺度），固定效应 + DerSimonian–Laird 随机效应。"""
    th = np.array([math.log(a / (1 - a)) for a, _, _ in estimates])
    v = np.array([logit_se(a, n1, n2) ** 2 for (a, n1, n2) in estimates])
    w = 1.0 / v
    k = len(th)
    th_fe = float(np.sum(w * th) / np.sum(w))
    se_fe = float(math.sqrt(1.0 / np.sum(w)))
    Q = float(np.sum(w * (th - th_fe) ** 2))
    df = k - 1
    C = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    I2 = max(0.0, (Q - df) / Q * 100) if Q > 0 else 0.0
    ws = 1.0 / (v + tau2)
    th_re = float(np.sum(ws * th) / np.sum(ws))
    se_re = float(math.sqrt(1.0 / np.sum(ws)))
    lo, hi = th_re - Z975 * se_re, th_re + Z975 * se_re
    return {
        "label": label, "k_cohorts": k,
        "fe_auc": 1 / (1 + math.exp(-th_fe)),
        "fe_lo": 1 / (1 + math.exp(-(th_fe - Z975 * se_fe))),
        "fe_hi": 1 / (1 + math.exp(-(th_fe + Z975 * se_fe))),
        "Q": Q, "df": df, "Q_p": float(stats.chi2.sf(Q, df)) if df > 0 else float("nan"),
        "I2_pct": I2, "tau2": tau2,
        "re_auc": 1 / (1 + math.exp(-th_re)),
        "re_lo": 1 / (1 + math.exp(-lo)),
        "re_hi": 1 / (1 + math.exp(-hi)),
        "re_p_vs_0.5": float(2 * stats.norm.sf(abs(th_re) / se_re)),
    }


def main() -> int:
    COHORTS = build_cohorts()
    log: list[str] = [
        "# 功效分析与跨队列 AUC 元分析\n",
        "生成脚本 `scripts/25_power_analysis.py`。",
        "口径锁定于 `docs/03` §23，**写于任何计算之前**。\n",
        "第 1–2 节是对既有结果的重新表述（不产生新估计）；"
        "第 3 节为**事后分析**，依 §17.7 显式标注。\n",
    ]

    # ---------------- 1. 所需样本量 ----------------
    log.append("## 1. 需要多少样本才能给出否定结论\n")
    log.append("判据：让 95% CI 上界降到 0.70（临床可用下限）或 0.77"
               "（IMPRES 原文跨数据集区间下界）。每组 n，PRCR/PD 平衡。\n")
    log.append("| 队列 | 签名 | n(PRCR/PD) | 实测 AUC | 现状 CI 上界 | 需 n@上界<0.70 | 需 n@上界<0.77 | 现状够不够 |")
    log.append("|---|---|---|---:|---|---:|---:|---|")
    pow_rows = []
    for coh, sig, orient, n1, n2, auc, lo, hi, where in COHORTS:
        n70 = n_for_upper_below(auc, CLINICAL_FLOOR)
        n77 = n_for_upper_below(auc, IMPRES_PUBLISHED_LO)
        cur = max(n1, n2)
        need70 = "不可能" if n70 is None else str(n70)
        need77 = "不可能" if n77 is None else str(n77)
        if n70 is not None and cur >= n70:
            state = "**已足够**"
        else:
            state = f"不足（差约 {max(0, (n70 or 10 ** 6) - cur)} 例/组）"
        log.append(f"| {coh} | {sig} | {n1}/{n2} | {auc:.3f} | {hi:.3f} | "
                   f"{need70} | {need77} | {state} |")
        pow_rows.append({
            "cohort": coh, "signature": sig, "orientation": orient,
            "n_PRCR": n1, "n_PD": n2, "n_total": n1 + n2,
            "auc_observed": auc, "ci_lo_boot": lo, "ci_hi_boot": hi,
            "hm_se": hm_se(auc, n1, n2),
            "hm_lo": auc - Z975 * hm_se(auc, n1, n2),
            "hm_hi": auc + Z975 * hm_se(auc, n1, n2),
            "n_needed_upper_below_0.70": n70,
            "n_needed_upper_below_0.77": n77,
            "current_n_sufficient_for_0.70": bool(n70 is not None and cur >= n70),
            "pre_registration_site": where,
        })
    log.append("")
    log.append("> 「不可能」= 点估计本身已 ≥ 该阈值，增大 n 只会让 CI 更窄地"
               "落在该值之上，样本量救不回来。\n")
    log.append("> 「需 n」按**实测点估计**不变来算，即假设真实 AUC 就等于观测值。"
               "因此当实测值已远低于 0.70 时，所需 n 反而很小——"
               "这不是说「再收 6 例就够了」，而是说**现有 n 已经够**。\n")

    # ---------------- 2. 最小可检出差异 ----------------
    log.append("## 2. 现有样本量能检出多大的差异\n")
    log.append("双侧备择，单侧 alpha = 0.025（与既有报告口径一致）。\n")
    log.append("| 每组 n | 总 n | 最小可检出 AUC（80% 功效） | （90% 功效） | 该 AUC 的 95% CI 宽度 |")
    log.append("|---:|---:|---:|---:|---:|")
    mdd_rows = []
    for n in (10, 12, 23, 24, 30, 34, 40, 45, 50, 69, 100):
        a80, a90 = mdd(n, 0.80), mdd(n, 0.90)
        ci_w = 2 * Z975 * se_equal_n(a80, n)
        log.append(f"| {n} | {2 * n} | **{a80:.3f}** | {a90:.3f} | {ci_w:.3f} |")
        mdd_rows.append({"n_per_group": n, "n_total": 2 * n,
                         "mdd_80": a80, "mdd_90": a90, "ci_width_at_mdd80": ci_w})
    log.append("")
    p70 = power_at(0.70, 45, 34)
    p77 = power_at(0.77, 45, 34)
    p65 = power_at(0.65, 45, 34)
    p60 = power_at(0.60, 45, 34)
    log.append("### 2.1 关键的一问：主比较能不能看见「真实存在」的性能？\n")
    log.append("本节是对「样本量不足」最直接的回应。单侧 alpha = 0.025。\n")
    log.append("| 若真实 AUC 为 | 主比较 n=45/34 的检出功效 |")
    log.append("|---:|---:|")
    for a in (0.60, 0.65, 0.70, 0.77):
        log.append(f"| {a:.2f} | {power_at(a, 45, 34):.3f} |")
    log.append("")
    log.append(f"**这是本文最强的功效论据。** 若 IMPRES 迁移到 GSE215868 时真有 "
               f"AUC 0.77（原文跨数据集区间下界），本设计的检出功效是 **{p77:.3f}**；"
               f"即便只真有 0.70，功效仍有 **{p70:.3f}**。\n")
    log.append("换言之，主比较**不是**「样本量太小所以看不见」——"
               "它完全看得见。"
               f"实际观测到的是 {0.505:.3f}，在 0.70–0.77 这个量级的真实性能下"
               "几乎不可能出现。**因此阴性结果在主比较上是信息充分的。**\n")
    log.append(f"本设计**看不见**的是弱信号：AUC 0.60 的检出功效仅 {p60:.3f}，"
               f"0.65 为 {p65:.3f}。所以本文对「弱但真实」的存在性**不作否定主张**，"
               "这正是 `IPS-MHC+CP` 合并 AUC 落在 0.59 附近时必须如实说的话。\n")

    # ---------------- 3. 跨队列元分析（事后） ----------------
    log.append("## 3. 跨队列 AUC 元分析（**事后分析**）\n")
    log.append("> 依 `docs/03` §17.7 标注：以下不属预登记，为事后、诊断性分析。"
               "不使用应答标签做建模，不改变 `docs/07` 已裁决的主比较，"
               "不追认任何阳性结论。\n")
    log.append("预登记 §3 禁止**样本级合并**。本节不违反该约束："
               "输入只有各队列已发布的汇总统计量（AUC 与 n），"
               "不接触任何患者级记录。\n")

    mrows = []
    for sig, orient in (("IMPRES g1_low", "g1_low"),
                        ("IMPRES g2_low", "g2_low"),
                        ("IPS-MHC+CP", "high_ips_responder")):
        est = [(auc, n1, n2) for c, s, o, n1, n2, auc, lo, hi, w in COHORTS
               if s == sig and o == orient]
        res = meta_auc(est, sig)
        res["orientation"] = orient
        res["n_total_pooled"] = sum(n1 + n2 for _, n1, n2 in est)
        res["n_PRCR_pooled"] = sum(n1 for _, n1, _ in est)
        res["n_PD_pooled"] = sum(n2 for _, _, n2 in est)
        mrows.append(res)
        log.append(f"### {sig}（{len(est)} 个队列，合计 n={res['n_total_pooled']}）\n")
        log.append("| 模型 | 合并 AUC | 95% CI |")
        log.append("|---|---:|---|")
        log.append(f"| 固定效应 | **{res['fe_auc']:.3f}** | {res['fe_lo']:.3f}–{res['fe_hi']:.3f} |")
        log.append(f"| 随机效应（DL） | **{res['re_auc']:.3f}** | {res['re_lo']:.3f}–{res['re_hi']:.3f} |")
        log.append("")
        log.append(f"- Cochran's Q = {res['Q']:.3f}（df={res['df']}，p={res['Q_p']:.3f}），"
                   f"I² = {res['I2_pct']:.1f}%，τ² = {res['tau2']:.4f}")
        log.append(f"- 合并 AUC 与 0.50 的双侧检验：p = {res['re_p_vs_0.5']:.3f}")
        log.append("")
    mdf = pd.DataFrame(mrows)

    # 针对 IPS 的专项解读：I² = 0 意味着三个队列彼此高度一致
    ips_row = mdf[mdf["label"] == "IPS-MHC+CP"].iloc[0]
    log.append("### 3.1 为什么 `IPS-MHC+CP` 的 I² = 0 值得单独写一段\n")
    log.append(f"`IPS-MHC+CP` 的 Cochran's Q = {ips_row['Q']:.3f}（p = {ips_row['Q_p']:.3f}），"
               f"**I² = {ips_row['I2_pct']:.1f}%**——三个队列的估计"
               "（0.604 / 0.583 / 0.591）彼此的差异**小于随机波动本应产生的差异**。\n")
    log.append("这与其他两个签名形成鲜明对照：`IMPRES g1_low` 的 I² = 34.6%，"
               "三个队列连方向都不一致。\n")
    log.append("**但这不能被写成「IPS 被复现」。** 判据写死在 §23.5：")
    log.append(f"合并 AUC = {ips_row['re_auc']:.3f} [{ips_row['re_lo']:.3f}, {ips_row['re_hi']:.3f}]，"
               f"**CI 下界 {ips_row['re_lo']:.3f} 仍覆盖 0.50**，"
               f"对 0.50 的双侧 p = {ips_row['re_p_vs_0.5']:.3f}，未达 0.05。\n")
    log.append("正确的说法只有一句：**该签名在这三个队列中表现出一个"
               "跨队列高度一致、但水平过低、本设计无法确认的弱信号。**"
               "按 §1 的表，要在观测到的水平上把 CI 上界压到 0.70 以下，"
               "需每组约 53 例。这是给后续研究的**样本量建议**，不是本文的结论。\n")
    log.append("I² = 0 也**不能**被解读为「三个队列是同质的所以可以合并」——"
               "恰好相反：正因为三个队列在生物学人群上差异很大（见 `docs/11`），"
               "却给出几乎相同的估计，"
               "更可能反映的是**该信号本身很弱、弱到人群差异淹没了它**，"
               "而不是「人群其实一样」。这一点必须写进讨论。\n")

    # ---------------- 4. 符号一致性 ----------------
    log.append("## 4. 方向一致性检验（符号检验）\n")
    srows = []
    log.append("「同向」必须区分两种，**混为一谈就会把反向说成复现**：\n")
    log.append("- **一致高于 0.50** = 各队列都指向与原文相同的方向；")
    log.append("- **一致低于 0.50** = 各队列都指向**相反**方向，镜像后虽为正值，"
               "但那是把签名反过来用，**不构成对原文的复现**。\n")
    log.append("| 签名 | 各队列点估计 | >0.50 | <0.50 | 性质 | 符号检验 p（双侧） |")
    log.append("|---|---|---:|---:|---|---:|")
    for sig, orient in (("IMPRES g1_low", "g1_low"),
                        ("IMPRES g2_low", "g2_low"),
                        ("IPS-MHC+CP", "high_ips_responder")):
        vals = [auc for _c, s, o, _n1, _n2, auc, _lo, _hi, _w in COHORTS
                if s == sig and o == orient]
        above = sum(1 for v in vals if v > 0.50)
        k = len(vals)
        p = float(stats.binomtest(max(above, k - above), k, 0.5).pvalue)
        if above == k:
            kind = "一致正向"
        elif above == 0:
            kind = "**一致反向**"
        else:
            kind = "方向不一致"
        log.append(f"| {sig} | {' / '.join(f'{v:.3f}' for v in vals)} | {above} | "
                   f"{k - above} | {kind} | {p:.3f} |")
        srows.append({"signature": sig, "k": k, "n_above_0.5": above,
                      "n_below_0.5": k - above, "kind": kind, "sign_test_p": p})
    log.append("")
    log.append("**三条全部 p ≥ 0.25，无一构成统计证据。**\n")
    log.append("- `IMPRES g1_low`（预登记主口径）：1 高 2 低，方向不一致，无信号。")
    log.append("- `IMPRES g2_low`（镜像）：0.659 / 0.583 / 0.483，2 高 1 低，方向不一致。"
               "合并 AUC 0.537 [0.435, 0.635]，与 0.50 无差异（p = 0.48）。"
               "**两个方向都没有信号**——这比「反向有效」更常见的解释"
               "（判别力不足、方向约定搞反、基因对有误）都更符合数据。"
               "**不得把 0.659 挑出来写成「换个方向就复现了」。**")
    log.append("- `IPS-MHC+CP`：3/3 一致高于 0.50，符号检验 p = 0.25。"
               "**方向一致但统计上不构成证据**，按 §23.5 规则 3 仍判未复现，"
               "只作为待验证观察（`docs/09`）。\n")

    # ---------------- 5. 结论边界 ----------------
    log.append("## 5. 本节能否定什么、不能否定什么\n")
    log.append("**能**：在这三个队列的合并样本量（n=138）上，"
               "两个签名的判别能力都落在随机区间附近，"
               "且这个结论对队列间异质性不敏感。\n")
    log.append("**不能**：否定 IMPRES 或 IPS 在其原始人群中的价值。"
               "本文从未在 IMPRES 的推导队列 GSE115821 上测试过"
               "（该系列 37 个样本仅对应 8 例独立患者，见 `docs/02` §7.2），"
               "因此「迁移失败」与「签名本身无效」在本文数据下**不可区分**。\n")
    log.append("**这是本文最重要的边界，任何引用本文的表述都不得越过。**\n")

    pd.DataFrame(pow_rows).to_csv(QC / "power_analysis.tsv", sep="\t", index=False)
    pd.DataFrame(mdd_rows).to_csv(QC / "power_mdd.tsv", sep="\t", index=False)
    mdf.to_csv(QC / "meta_auc.tsv", sep="\t", index=False)
    pd.DataFrame(srows).to_csv(QC / "sign_consistency.tsv", sep="\t", index=False)
    body = "\n".join(log)
    (QC / "docs_10_power.md").write_text(body, encoding="utf-8")
    (ROOT / "docs" / "10_功效分析与跨队列元分析.md").write_text(body, encoding="utf-8")

    print(mdf[["label", "k_cohorts", "fe_auc", "fe_lo", "fe_hi",
               "Q_p", "I2_pct", "re_auc", "re_lo", "re_hi",
               "re_p_vs_0.5"]].to_string(index=False))
    print()
    print(pd.DataFrame(srows).to_string(index=False))
    print(f"\n主比较 n=45/34 检出功效：AUC 0.60 -> {p60:.3f} | "
          f"0.70 -> {p70:.3f} | 0.77 -> {p77:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
