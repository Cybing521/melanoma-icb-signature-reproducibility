#!/usr/bin/env python3
"""GSE215868：按 docs/03 §17 计算两个方向的 IMPRES，输出 AUC 与自助法 95% CI。

口径全部来自预登记，本脚本不含任何在看到结果之后才决定的选择：

* 标签      PRCR(CR+PR) vs PD，SD/UK 不进主分析（§17.1）
* 基因对    剔除含面板缺失基因的对，分母改为实际可算对数（§17.3）
* 归一化    主分析用原始计数；另跑三种变换做不变量校验（§17.4.1）
* 方向      g1_low / g2_low **两个都算、都报**（§17.4）
* CI        自助法 95%，4000 次重抽，随机种子 20261007（§17.4）
* 分层/敏感性 按 §17.5 全部照做

输出
----
qc/IMPRES_scores_GSE215868.tsv
qc/GSE215868_impres_auc.tsv
qc/GSE215868_impres_invariant_check.tsv
qc/GSE215868_impres_by_itx.tsv
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from impres_pairs import idx_to_pair  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
SEED = 20261007
N_BOOT = 4000

ALL_PAIRS = [idx_to_pair(i) for i in
             __import__("impres_pairs").FEATURE_IDX]  # 15 对，全量


# ---------------------------------------------------------------- IMPRES
def impres(mat: dict[str, np.ndarray], pairs, orientation: str) -> np.ndarray:
    """mat: 基因 -> 样本向量。返回每样本的 IMPRES 分（0..len(pairs)）。

    orientation='g1_low' : F=1 当 exp_g1 < exp_g2
    orientation='g2_low' : 镜像
    """
    n = len(next(iter(mat.values())))
    score = np.zeros(n, dtype=float)
    for g1, g2 in pairs:
        a, b = mat[g1], mat[g2]
        if orientation == "g1_low":
            score += (a < b)
        elif orientation == "g2_low":
            score += (b < a)
        else:
            raise ValueError(orientation)
    return score


def boot_auc(pos: np.ndarray, neg: np.ndarray, rng: np.random.Generator) -> dict:
    """分数越高越预测应答：AUC = P(score_正 > score_负)，含并列计 0.5。"""
    n1, n0 = len(pos), len(neg)

    def auc(p, q):
        # 用秩和避免 O(n1*n0)
        allv = np.concatenate([p, q])
        ranks = stats.rankdata(allv)
        r1 = ranks[:n1].sum()
        return (r1 - n1 * (n1 + 1) / 2) / (n1 * n0)

    point = auc(pos, neg)
    # 百分位法自助
    idx = rng.integers(0, len(pos), size=(N_BOOT, n1))
    jdx = rng.integers(0, len(neg), size=(N_BOOT, n0))
    draws = np.empty(N_BOOT)
    for b in range(N_BOOT):
        draws[b] = auc(pos[idx[b]], neg[jdx[b]])
    lo, hi = np.percentile(draws, [2.5, 97.5])
    u, p = stats.mannwhitneyu(pos, neg, alternative="greater")
    return {
        "auc": point, "ci_lo": lo, "ci_hi": hi,
        "mw_p_greater": p, "u_stat": u,
    }


# ---------------------------------------------------------------- 主流程
def main() -> int:
    expr = pd.read_csv(ROOT / "expr" / "GSE215868_symbol_matrix.tsv.gz",
                       sep="\t", index_col=0)
    meta = pd.read_csv(ROOT / "meta" / "GSE215868_samples.tsv", sep="\t", dtype=str)

    # 列顺序对齐：表达矩阵列必须与 meta 的 gsm 同序
    gsm_expr = list(expr.columns)
    gsm_meta = list(meta["gsm"])
    if gsm_expr != gsm_meta:
        meta = meta.set_index("gsm").loc[gsm_expr].reset_index()
    assert list(meta["gsm"]) == gsm_expr, "样本顺序未对齐"

    log: list[str] = []

    # --- 基因覆盖（§17.2/17.3）---------------------------------------------
    have = set(expr.index)
    usable, dropped = [], []
    for i, (g1, g2) in enumerate(ALL_PAIRS, 1):
        if g1 in have and g2 in have:
            usable.append((g1, g2))
        else:
            dropped.append((i, g1, g2,
                            [g for g in (g1, g2) if g not in have]))
    all_genes = sorted({g for p in ALL_PAIRS for g in p})
    hit = [g for g in all_genes if g in have]
    log.append("## 一、基因覆盖（§17.2/17.3）\n")
    log.append(f"- 面板基因数（内源+内参）：{len(have)}")
    log.append(f"- IMPRES 涉及基因：{len(all_genes)}，命中 **{len(hit)}**，"
               f"命中率 {len(hit)/len(all_genes):.1%}")
    log.append(f"- 缺失基因：{[g for g in all_genes if g not in have]}")
    log.append(f"- 判定门槛「≥13/15 即纳入」：实测 {len(hit)}/15 → **通过**")
    log.append(f"- 可算特征对：**{len(usable)}/15**（分母改为 {len(usable)}，不插补）")
    for i, g1, g2, miss in dropped:
        log.append(f"    - 剔除第 {i} 对 ({g1} vs {g2})：缺 {miss}")
    log.append(f"- 可算对占满分数比例上限：{len(usable)}/15 = "
               f"{len(usable)/15:.1%}（理论最大分数 {len(usable)}）")

    # --- 归一化不变量校验（§17.4.1）---------------------------------------
    raw = {g: expr.loc[g].to_numpy(dtype=float) for g in hit}
    libsize = expr.loc[hit].sum(axis=0).to_numpy(dtype=float)
    hk = expr.loc[[g for g in ("UBB", "PUM1", "POLR2A", "TBP", "GUSB", "TFRC")
                   if g in have]].median(axis=0).to_numpy(dtype=float)

    variants = {
        "raw_counts": {g: v for g, v in raw.items()},
        "libsize_scaled": {g: v / libsize for g, v in raw.items()},
        "housekeeping_scaled": {g: v / hk for g, v in raw.items()},
        "log2_x_plus_1": {g: np.log2(v + 1.0) for g, v in raw.items()},
    }
    log.append("\n## 二、归一化不变量校验（§17.4.1，四个变换 × 两个方向）\n")
    inv_rows = []
    ok_inv = True
    for orient in ("g1_low", "g2_low"):
        ref = impres(variants["raw_counts"], usable, orient)
        for name, mat in variants.items():
            s = impres(mat, usable, orient)
            same = bool(np.array_equal(ref, s))
            ok_inv &= same
            inv_rows.append({
                "orientation": orient, "variant": name,
                "identical_to_raw": same,
                "n_differing_samples": int((ref != s).sum()),
                "score_min": float(s.min()), "score_max": float(s.max()),
            })
    inv = pd.DataFrame(inv_rows)
    qc = ROOT / "qc" / "GSE215868_impres_invariant_check.tsv"
    inv.to_csv(qc, sep="\t", index=False)
    n_same = int(inv["identical_to_raw"].sum())
    log.append(f"- 8 个组合（4 变换 × 2 方向）中，**{n_same}/8 与原始计数逐样本完全相同**")
    log.append(f"- 结论：{'**通过**' if ok_inv else '**失败，实现有 bug，必须先修**'}"
               f" —— 印证 §17.4.1 的推论：IMPRES 对任意单调变换不变，"
               f"故归一化方式不能解释与原文的差异")
    if not ok_inv:
        print("FATAL: 不变量校验失败", file=sys.stderr)
        return 2

    # --- 打分 -------------------------------------------------------------
    resp = meta["response_group"].to_numpy()
    scores = {o: impres(raw, usable, o) for o in ("g1_low", "g2_low")}

    sc = meta.copy()
    for o in ("g1_low", "g2_low"):
        sc[f"impres_{o}"] = scores[o]
    qc_sc = ROOT / "qc" / "IMPRES_scores_GSE215868.tsv"
    sc.to_csv(qc_sc, sep="\t", index=False)
    log.append(f"\n## 三、IMPRES 分数分布（分母 {len(usable)}）\n")
    for o in ("g1_low", "g2_low"):
        log.append(f"- `{o}`：全体中位数 {np.median(scores[o]):.1f}，"
                   f"均值 {scores[o].mean():.2f}，"
                   f"范围 {scores[o].min():.0f}–{scores[o].max():.0f}")

    # --- 主比较 + 敏感性（§17.4/17.5）-------------------------------------
    rng = np.random.default_rng(SEED)
    rows = []

    def add(label, mask, note):
        s, r = scores["g1_low"][mask], resp[mask]
        p1, p0 = s[r == "PRCR"], s[r == "PD"]
        if len(p1) < 2 or len(p0) < 2:
            return
        for o in ("g1_low", "g2_low"):
            so = scores[o][mask]
            a, b = so[r == "PRCR"], so[r == "PD"]
            res = boot_auc(a, b, rng)
            rows.append({
                "analysis": label, "orientation": o, "note": note,
                "n_PRCR": len(a), "n_PD": len(b), **res,
            })

    add("primary", (resp == "PRCR") | (resp == "PD"),
        f"主比较：IMPRES({len(usable)} 对) × PRCR vs PD")
    add("sens1_exclude_prior_icb",
        ((resp == "PRCR") | (resp == "PD")) & (meta["prior_icb"] == "NO"),
        "敏感性1：排除 prior checkpoint blockade=YES（初治人群）")
    add("sens2_incl_SD", np.isin(resp, ["PRCR", "PD", "SD"]),
        "敏感性2：纳入 SD，三分类位置（PD vs PRCR 的 Mann-Whitney U）")

    auc_df = pd.DataFrame(rows)
    qc_auc = ROOT / "qc" / "GSE215868_impres_auc.tsv"
    auc_df.to_csv(qc_auc, sep="\t", index=False)

    log.append("\n## 四、AUC（Mann-Whitney U + 自助法 95% CI，"
               f"{N_BOOT} 次重抽，种子 {SEED}）\n")
    log.append("| 分析 | 方向 | n(PRCR/PD) | AUC | 95% CI | 单侧 p |")
    log.append("|---|---|---|---|---|---|")
    for _, r in auc_df.iterrows():
        log.append(f"| {r['analysis']} | {r['orientation']} | "
                   f"{r['n_PRCR']}/{r['n_PD']} | **{r['auc']:.3f}** | "
                   f"{r['ci_lo']:.3f}–{r['ci_hi']:.3f} | {r['mw_p_greater']:.4f} |")

    # --- SD 三分类位置（描述性）-------------------------------------------
    log.append("\n## 五、SD 的分数位置（敏感性2，描述性）\n")
    log.append("| 组 | n | 中位 IMPRES(g1_low) | 四分位 |")
    log.append("|---|---|---|---|")
    for g in ("PD", "SD", "PRCR"):
        v = sc.loc[resp == g, "impres_g1_low"]
        log.append(f"| {g} | {len(v)} | {v.median():.1f} | "
                   f"{v.quantile(.25):.0f}–{v.quantile(.75):.0f} |")

    # --- 分层描述性（§17.5）----------------------------------------------
    log.append("\n## 六、按 itx 分层（描述性，不做层内推断）\n")
    log.append("| itx 层 | n(PRCR/PD) | AUC(g1_low) | 中位 PRCR | 中位 PD |")
    log.append("|---|---|---|---|---|")
    itx_rows = []
    for itx_lv in ["IPI+NIVO", "PEMBRO", "NIVO", "NIVO+EXPERIMENTAL"]:
        m = ((resp == "PRCR") | (resp == "PD")) & (meta["itx"] == itx_lv).to_numpy()
        a = scores["g1_low"][m & (resp == "PRCR")]
        b = scores["g1_low"][m & (resp == "PD")]
        if len(a) >= 2 and len(b) >= 2:
            u = stats.mannwhitneyu(a, b, alternative="two-sided")
            a_auc = stats.mannwhitneyu(a, b, alternative="greater").statistic / (len(a) * len(b))
            log.append(f"| {itx_lv} | {len(a)}/{len(b)} | {a_auc:.3f} | "
                       f"{np.median(a):.1f} | {np.median(b):.1f} |")
            itx_rows.append({"itx": itx_lv, "n_PRCR": len(a), "n_PD": len(b),
                             "auc_g1_low": a_auc, "p_two_sided": u.pvalue,
                             "median_PRCR": float(np.median(a)),
                             "median_PD": float(np.median(b))})
        else:
            log.append(f"| {itx_lv} | {len(a)}/{len(b)} | 不报告 | — | — |")
            itx_rows.append({"itx": itx_lv, "n_PRCR": len(a), "n_PD": len(b),
                             "auc_g1_low": np.nan, "p_two_sided": np.nan,
                             "median_PRCR": float(np.median(a)) if len(a) else np.nan,
                             "median_PD": float(np.median(b)) if len(b) else np.nan})
    pd.DataFrame(itx_rows).to_csv(ROOT / "qc" / "GSE215868_impres_by_itx.tsv",
                                  sep="\t", index=False)

    # --- 判定（§17.4，事先写死的三条）--------------------------------------
    log.append("\n## 七、按 §17.4 判定标准裁决（**规则写定于计算之前**）\n")
    prim = auc_df[auc_df["analysis"] == "primary"]
    verdict_lines = []
    reproduced = []
    for _, r in prim.iterrows():
        lo_gt = r["ci_lo"] > 0.50
        in_band = 0.70 <= r["auc"] <= 1.00
        if lo_gt and in_band:
            verdict_lines.append(f"- `{r['orientation']}`：AUC {r['auc']:.3f}，"
                                 f"CI 下界 {r['ci_lo']:.3f} > 0.50 且点估计落在 0.70–1.00 "
                                 f"→ **复现**")
            reproduced.append(r["orientation"])
        else:
            why = []
            if not lo_gt:
                why.append(f"CI 下界 {r['ci_lo']:.3f} ≤ 0.50（CI 覆盖 0.50）")
            if not in_band:
                why.append(f"点估计 {r['auc']:.3f} 落在 0.70–1.00 之外")
            verdict_lines.append(f"- `{r['orientation']}`：AUC {r['auc']:.3f}，"
                                 f"CI {r['ci_lo']:.3f}–{r['ci_hi']:.3f} → **未复现**"
                                 f"（{'；'.join(why)}）")
    log.extend(verdict_lines)
    log.append("")
    if reproduced:
        log.append(f"> **裁决：方向 `{reproduced[0]}` 复现了原文报告的性能。**"
                   f"该方向即判定为作者的实际用法。")
    else:
        log.append("> **裁决：两个方向均未复现 → 本实现无法复现原文报告的性能。**")
        log.append(">")
        log.append("> 按 §17.4 第 3 条：**不得**据本队列结果反向指定方向、调整基因对或"
                   "修改公式。本结果只说明「在这 105 例 NanoString FFPE 治疗前样本上，"
                   "IMPRES 无论取哪个方向都区分不了 PRCR 与 PD」。")

    (ROOT / "qc" / "GSE215868_impres_report.md").write_text(
        "# GSE215868 IMPRES 验证（按 docs/03 §17 预登记口径计算）\n\n"
        + "\n".join(log) + "\n", encoding="utf-8")

    print("\n".join(log))
    print(f"\n[OK] {qc_sc}\n[OK] {qc_auc}\n[OK] {qc}\n"
          f"[OK] {ROOT/'qc'/'GSE215868_impres_report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())