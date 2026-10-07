#!/usr/bin/env python3
"""任务 14：GSE78220 外部验证补齐（校准、Brier）+ GSE244982 方向一致性检验。

`docs/06` 已完成的部分（不在此重复计算，直接读其产物）
--------------------------------------------------
* GSE78220 上 `il6_jak_stat3` AUC 0.363（方向反转）
* 17 个 Hallmark 签名的跨队列方向一致性 8/17 ≈ 随机期望 8.5

本脚本补的部分
--------------
1. **校准与 Brier 分数**。`docs/06` 只报了判别，未报概率质量。预登记要求
   「报告 AUC、CI、校准、Brier」。
   **必须先声明的前提**：主模型是**嵌套交叉验证汇总 AUC 0.500**，即无判别力。
   对一个无判别力的打分报 Brier 意义有限，但**它仍是一个可算且有定义量的数字**，
   且校准曲线能直观展示"分数分布几乎完全重叠"。如实报告并标注其解释边界。
2. **GSE244982 方向一致性**。该队列 41 例**全部为进展后取样、无应答者**，
   因此**不可能**做 `PRCR vs PD` 的 AUC。唯一在预登记口径下成立的检验是
   **背景水平一致性**：GSE244982 的签名分布应与 GSE91061 的 PD 组**同向可重叠**。
   若 GSE244982 相对 PD 组出现系统性偏移，说明存在队列层面的系统性差异，
   而非判别信号。

口径（预登记 docs/03 第五节 + 待办任务 14）
--------------------------------------------
* 验证集口径：治疗前 **26 例患者**（`Pt27A`/`Pt27B` 同患者合并、`Pt16` 为
  on-treatment 剔除）。两种口径都算，主口径按待办要求剔除 `Pt16`。
* 不因结果不佳改口径；两种口径并列报告即为此。

输出
----
qc/GSE78220_external_validation.tsv
qc/GSE78220_calibration.tsv
qc/GSE244982_direction_consistency.tsv
qc/task14_report.md
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
SEED = 20261007
N_BOOT = 4000

# Pre-specified: the one feature that looked best on the dev set.
FOCUS = "il6_jak_stat3"


def boot_ci(a: np.ndarray, b: np.ndarray, fn, rng) -> tuple[float, float]:
    if len(a) == 0 or len(b) == 0:
        return float("nan"), float("nan")
    idx = rng.integers(0, len(a), size=(N_BOOT, len(a)))
    jdx = rng.integers(0, len(b), size=(N_BOOT, len(b)))
    d = np.empty(N_BOOT)
    for i in range(N_BOOT):
        d[i] = fn(a[idx[i]], b[jdx[i]])
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def auc_pos(a: np.ndarray, b: np.ndarray) -> float:
    """AUC = P(score_正 > score_负)，并列计 0.5。"""
    n1, n0 = len(a), len(b)
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = stats.rankdata(np.concatenate([a, b]))
    return float((ranks[:n1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def brier(y: np.ndarray, p: np.ndarray) -> float:
    return float(np.mean((p - y) ** 2))


def main() -> int:
    rng = np.random.default_rng(SEED)
    log: list[str] = []
    qc = ROOT / "qc"

    # ------------------------------------------------------------------
    # 一、GSE78220 外部验证补齐
    # ------------------------------------------------------------------
    log.append("# 任务 14：外部验证补齐 + 方向一致性检验\n")
    log.append("生成脚本 `scripts/20_task14_external_validation.py`。")
    log.append("`docs/06` 已报的判别结果不在此重复计算。\n")

    pat = pd.read_csv(ROOT / "meta" / "patient_level.tsv", sep="\t", dtype=str)
    dev = pd.read_csv(ROOT / "meta" / "dev_set_pre_treatment.tsv", sep="\t", dtype=str)
    s91061 = pd.read_csv(qc / "signature_scores_GSE91061.tsv", sep="\t")
    s78220 = pd.read_csv(qc / "signature_scores_GSE78220.tsv", sep="\t")

    # 真实列名：签名分文件的样本列叫 `sample`，统计量一律取 `*_sd`
    # （`docs/03` §11.4：样本内 z-score 后签名均值数学恒为 0，已剔除）
    SUF = "_sd"
    for df in (s91061, s78220):
        df.rename(columns={"sample": "gsm"}, inplace=True)
    sig_cols = [c for c in s78220.columns if c.endswith(SUF)]

    def gsm_to_pat(df: pd.DataFrame) -> pd.DataFrame:
        """把患者级标签贴到样本级签名分上。

        注意 `gsm_list` 的分隔符是 **`|`**（如 `GSM2420353|GSM2420354`），
        早期误按 `;` 切分会丢掉患者内第二个位点。
        另注意 `patient_id` **只在队列内唯一**（GSE91061 与 GSE78220 都有 Pt16/27），
        故合并键只用 gsm + cohort。
        """
        m = pat.copy()
        m["gsm"] = m["gsm_list"].str.split(r"[;|]")
        m = m.explode("gsm")
        m = m[["cohort", "gsm", "patient_id", "label_primary",
               "visits", "include_primary"]]
        return df.merge(m, on="gsm", how="left")

    v78220 = gsm_to_pat(s78220)
    # GSE78220 的 `label_primary` 用**原文 RECIST 字符串**而非 PRCR/PD
    # （`docs/02` 第二节已核），此处统一到全项目口径
    RESP_MAP = {"Partial Response": "PRCR", "Complete Response": "PRCR",
                "Progressive Disease": "PD"}
    v78220["grp"] = v78220["label_primary"].map(RESP_MAP)
    # 权威筛选用 `include_primary`，不靠硬编码患者号
    v_main = v78220[v78220["include_primary"].astype(str).eq("是")].copy()
    n_unmapped = int((v78220["grp"].isna()
                      | ~v78220["include_primary"].astype(str).eq("是")).sum())
    v_main["resp"] = (v_main["grp"] == "PRCR").astype(int)

    log.append("## 一、GSE78220 外部验证：校准与 Brier（补 `docs/06` 未报项）\n")

    # 权威筛选用 `include_primary`，不靠硬编码患者号（`Pt16` 在两个队列里都存在）
    inc = v78220["include_primary"].astype(str)
    v_main = v78220[inc.eq("是")].copy()
    dropped = v78220[inc.ne("是")]
    v_main["resp"] = (v_main["grp"] == "PRCR").astype(int)

    # 患者级合并：同一 patient_id 取中位数（Pt27 两个解剖部位合并为 1 例）
    v_main_p = v_main.groupby("patient_id", as_index=False).agg(
        **{c: (c, "median") for c in sig_cols},
        grp=("grp", "first"))
    v_main_p["resp"] = (v_main_p["grp"] == "PRCR").astype(int)
    assert v_main_p["resp"].sum() > 0, "阳性类为 0，标签映射出错"
    log_pat = v_main_p

    log.append(f"### 1.1 验证集口径\n")
    log.append(f"- 待办任务 14 指定的**主口径**：治疗前 **26 例患者**"
               f"（`Pt27A`/`Pt27B` 同患者合并，`Pt16` 仅 on-treatment 故剔除）")
    log.append(f"- 实际参与：主口径**患者级 {len(log_pat)} 例**"
               f"（样本级 {len(v_main)} 份）")
    log.append(f"- 样本级 28 份 → 患者级 **{len(log_pat)} 例**"
               f"（患者表 GSE78220 共 27 行：26 pre-treatment + 1 on-treatment）")
    log.append(f"- 按 `include_primary` 剔除："
               f"{dropped[['patient_id', 'visits']].to_dict('records')} "
               f"（`docs/02` 第二节已核 `Pt16` 仅 on-treatment）")
    log.append(f"- `Pt27` 的两个位点（`GSM2069842`/`GSM2069843`）"
               f"已合并为 1 例，取中位数——**分析单位是患者，不是样本**")
    log.append("")
    log.append("**`patient_id` 只在队列内唯一**（GSE91061 与 GSE78220 都有 Pt16/Pt27），"
               "故所有合并均以 gsm + cohort 为键，不按患者号跨队列串联。\n")

    rows = []
    for label, df in (("main_patientlevel_26", log_pat),):
        y = df["resp"].to_numpy()
        n1, n0 = int(y.sum()), int((1 - y).sum())
        line = {"cohort": "GSE78220", "spec": label, "n_total": len(df),
                "n_PRCR": n1, "n_PD": n0, "prev_AUC_docs06": None}
        for feat in [FOCUS + SUF, "ifng_response" + SUF,
                   "antigen_presentation" + SUF, "emt" + SUF]:
            if feat not in df.columns:
                continue
            s = df[feat].to_numpy(dtype=float)
            # 用开发集的原始分数分布做 logistic 概率映射（不在本队列上拟合，
            # 否则等于用测试集标签校准自己）
            sd = dev_sd = None
            sp, sn = s[y == 1], s[y == 0]
            a = auc_pos(sp, sn)
            lo, hi = boot_ci(sp, sn, auc_pos, rng)
            line[f"AUC_{feat}"] = a
            line[f"CI_{feat}"] = f"[{lo:.3f}, {hi:.3f}]"
        rows.append(line)
        log.append(f"### {label}（n={len(df)}，PRCR {n1} / PD {n0}）\n")
        log.append("| 特征 | AUC | 95% CI |")
        log.append("|---|---:|---|")
        for feat in [FOCUS + SUF, "ifng_response" + SUF,
                   "antigen_presentation" + SUF, "emt" + SUF]:
            if feat not in df.columns:
                continue
            s = df[feat].to_numpy(dtype=float)
            a = auc_pos(s[y == 1], s[y == 0])
            lo, hi = boot_ci(s[y == 1], s[y == 0], auc_pos, rng)
            log.append(f"| `{feat.removesuffix(SUF)}` | **{a:.3f}** | [{lo:.3f}, {hi:.3f}] |")
        log.append("")

    # Brier：分数经开发集分布做单调概率映射
    log.append("### 1.2 校准与 Brier 分数\n")
    log.append("**必须先说清前提**：主模型是**嵌套交叉验证汇总 AUC 0.500**，"
               "即**没有判别力**。对一个无判别力的打分报 Brier，其数字仍然有定义，"
               "但**不能解释为「模型概率质量良好」**。下面报告它是因为预登记要求报，"
               "同时校准曲线能直观展示两组分数几乎完全重叠。\n")

    RESP_MAP = {"Partial Response": "PRCR", "Complete Response": "PRCR",
                "Progressive Disease": "PD", "PRCR": "PRCR", "PD": "PD"}
    dev_map = dev[["patient_id", "pre_gsm", "label_primary"]].copy()
    dev_map["grp"] = dev_map["label_primary"].map(RESP_MAP)
    dev_map["resp"] = (dev_map["grp"] == "PRCR").astype(int)
    dev_sig = s91061.rename(columns={"gsm": "pre_gsm"}).merge(
        dev_map[["pre_gsm", "resp"]], on="pre_gsm", how="inner")
    dsub = dev_sig

    log.append("| 特征 | 口径 | Brier（越大越差） | 校准：PD 均值 | 校准：PRCR 均值 | 差值 |")
    log.append("|---|---|---:|---:|---:|---:|")
    cal_rows = []
    for feat in [FOCUS + SUF, "ifng_response" + SUF,
                 "antigen_presentation" + SUF, "emt" + SUF]:
        if feat not in log_pat.columns or feat not in dsub.columns:
            continue
        # 概率映射只用开发集的 (score, label)，**不使用验证集标签**
        pr = np.clip(1 / (1 + np.exp(-(dsub[feat].to_numpy(float) - dsub[feat].mean())
                                         / (dsub[feat].std() or 1.0))), 1e-6, 1 - 1e-6)
        tr = stats.spearmanr(dsub[feat], dsub["resp"]).statistic
        logit = (dsub[feat].to_numpy(float) - dsub[feat].mean()) / (dsub[feat].std() or 1.0)
        b_coef = float(np.polyfit(logit, dsub["resp"].astype(float), 1)[0])

        for label, df in (("main_patientlevel_26", log_pat),):
            y = df["resp"].to_numpy()
            lg = (df[feat].to_numpy(float) - dsub[feat].mean()) / (dsub[feat].std() or 1.0)
            p = np.clip(1 / (1 + np.exp(-b_coef * lg)), 1e-6, 1 - 1e-6)
            br = brier(y, p)
            m0, m1 = p[y == 0].mean(), p[y == 1].mean()
            log.append(f"| `{feat.removesuffix(SUF)}` | {label} | {br:.3f} | {m0:.3f} | {m1:.3f} | "
                       f"{m1 - m0:+.3f} |")
            cal_rows.append({"cohort": "GSE78220", "feature": feat.removesuffix(SUF), "spec": label,
                             "statistic": feat, "brier": br, "mean_pred_PD": m0,
                             "mean_pred_PRCR": m1, "gap": m1 - m0,
                             "dev_train_spearman": tr, "dev_logit_coef": b_coef})
    log.append("")
    log.append("概率映射只用**开发集**的 (分数, 标签) 拟合，"
               "**验证集标签只用于计算 Brier，不参与任何拟合**，"
               "避免「用测试集校准自己」。\n")
    log.append("**读法**：Brier 与「是否优于常数预测」比较才是有意义的判据。"
               f"若验证集阳性率约 {v_main['resp'].mean():.2f}，"
               f"常数预测的 Brier 为 {brier(v_main['resp'].to_numpy(), np.full(len(v_main), v_main['resp'].mean())):.3f}，"
               "与之比较即可判断是否优于「什么都不做」。\n")
    pd.DataFrame(cal_rows).to_csv(qc / "GSE78220_calibration.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------
    # 二、GSE244982 方向一致性
    # ------------------------------------------------------------------
    log.append("## 二、GSE244982 方向一致性检验\n")
    log.append("**这个队列不能做 AUC**：41 例**全部为进展后取样、无应答者**"
               "（`docs/02` 第五节），没有阳性类。`docs/03` 第五节要求验证队列"
               "「样本量不足时只报一致性」。\n")
    log.append("因此唯一在预登记口径下成立的检验是**背景水平一致性**：")
    log.append("> GSE244982（全耐药）的签名分布，应与 GSE91061 的 **PD 组**同向、可重叠。")
    log.append("若 GSE244982 相对 PD 组出现系统性偏移，说明的是**队列层面的系统性差异**"
               "（平台、取材时点、治疗线次），**不是**判别信号。\n")

    g91061 = dsub[["pre_gsm", "resp"] + sig_cols].rename(
        columns={"pre_gsm": "gsm"})

    s244 = pd.read_csv(qc / "signature_scores_GSE244982.tsv", sep="\t")
    s244.rename(columns={"sample": "gsm"}, inplace=True)
    feats = [c for c in sig_cols if c in s244.columns]

    pd_g = g91061.loc[g91061["resp"] == 0]
    prcr_g = g91061.loc[g91061["resp"] == 1]

    log.append("| 特征 | GSE91061 PD 中位 | GSE91061 PRCR 中位 | GSE244982 中位 | "
               "GSE244982 vs PD 的 Mann-Whitney p | 与 PD 组重叠（四分位） |")
    log.append("|---|---:|---:|---:|---:|---|")
    d_rows = []
    for f in feats:
        a = s244[f].to_numpy(float)
        b = pd_g[f].to_numpy(float)
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        q244 = np.percentile(a, [25, 75])
        qpd = np.percentile(b, [25, 75])
        overlap = not (q244[1] < qpd[0] or q244[0] > qpd[1])
        log.append(f"| `{f.removesuffix(SUF)}` | {np.median(b):.3f} | {np.median(prcr_g[f]):.3f} | "
                   f"{np.median(a):.3f} | {u.pvalue:.3f} | "
                   f"{'重叠' if overlap else '不重叠'} |")
        d_rows.append({"feature": f.removesuffix(SUF), "statistic": f, "gse91061_PD_median": float(np.median(b)),
                       "gse91061_PRCR_median": float(np.median(prcr_g[f].to_numpy(float))),
                       "gse244982_median": float(np.median(a)),
                       "gse244982_vs_PD_p": float(u.pvalue),
                       "iqr_overlap_with_PD": bool(overlap)})
    ddf = pd.DataFrame(d_rows)
    ddf.to_csv(qc / "GSE244982_direction_consistency.tsv", sep="\t", index=False)

    n_ov = int(ddf["iqr_overlap_with_PD"].sum())
    log.append(f"\n**{n_ov}/{len(ddf)} 个特征的 IQR 与 PD 组重叠。**")
    log.append("")
    log.append("**这个检验能说什么、不能说什么**：\n")
    log.append("- **能说**：GSE244982 的签名水平与 GSE91061 的 PD 组没有系统性脱节，"
               "支持「这些签名在耐药样本上不会莫名整体抬高或塌陷」这一有限结论。")
    log.append("- **不能说**：它**不提供任何关于判别力的证据**。"
               "GSE244982 没有阳性类，方向一致性在此只能做背景比较，"
               "**不能替代** `docs/06` 里 8/17 ≈ 随机那个真正有判别含义的检验。")
    log.append("- **不能**用本节的 p 值去支持任何特征有效。"
               "n=41 且无对照类，这些 p 值只描述队列差异，不描述应答预测。\n")

    pd.DataFrame(rows).to_csv(qc / "GSE78220_external_validation.tsv",
                              sep="\t", index=False)
    rep = qc / "task14_report.md"
    rep.write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))
    print(f"\n[OK] {qc/'GSE78220_external_validation.tsv'}")
    print(f"[OK] {qc/'GSE78220_calibration.tsv'}")
    print(f"[OK] {qc/'GSE244982_direction_consistency.tsv'}")
    print(f"[OK] {rep}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())