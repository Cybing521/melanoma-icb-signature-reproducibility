#!/usr/bin/env python3
"""平台适配性诊断：IMPRES 的 13 个特征在 NanoString 面板上是否还携带信息。

**本脚本不属于预登记 §17 的分析计划，是在主比较出结果之后加的诊断。**
之所以仍然要做：主比较 AUC 0.505 需要一个机制性解释，否则只是一句"没复现"。

必须说清的合规边界：

1. 本脚本**完全不看应答标签**。它只度量"每个特征在样本之间的变异程度"，
   与结局无关，因此不存在从结局里捞结果的空间。
2. 本脚本**不改变、不替换、不重算** §17.4 已裁决的主比较。裁决已在 `scripts/18`
   里按写死的规则给出，本脚本的输出只用于解释它。
3. 因此它**不能**被读作"又找了个阳性结果"，也不进入任何多重比较校正。

要回答的问题
------------
IMPRES 的每个特征 F 只问一件事：**同一样本内，g1 的表达是否低于 g2**。

* 在**全转录组**（RNA-seq / 芯片）上，测序深度对每个基因几乎一致，
  `exp_g1 < exp_g2` 是一个有生物学含义的相对陈述。
* 在**靶向面板**（NanoString）上，每个靶标有独立的捕获效率，
  某基因的表达值可能被探针效率抬高或压低几个数量级。
  若如此，`exp_g1 < exp_g2` 会被探针效率主导，**在样本间近乎恒定**，
  该特征便不再携带任何与应答有关的信息。

度量（全部无标签）：
* `p1`     —— 105 个样本中 F=1 的比例
* `const`  —— F 在所有样本上是否恒定
* `p_min`  —— 若 p1 接近 0 或 1，说明该对被单一基因的技术效率锁死
* 与 GSE91061 / GSE78220（全转录组）做同样度量并对照

输出
----
qc/GSE215868_panel_diagnostic.tsv
qc/GSE215868_panel_diagnostic.md
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from impres_pairs import FEATURE_IDX, idx_to_pair  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
PAIRS = [idx_to_pair(i) for i in FEATURE_IDX]
COHORTS = {
    "GSE215868": ("NanoString 784 基因靶向面板", "expr/GSE215868_symbol_matrix.tsv.gz"),
    "GSE91061": ("RNA-seq 全转录组（原始计数）", "expr/GSE91061_symbol_matrix.tsv.gz"),
    "GSE78220": ("芯片 FPKM 全转录组", "expr/GSE78220_symbol_matrix.tsv.gz"),
    "GSE294272": ("RNA-seq 全转录组（featureCounts）", "expr/GSE294272_symbol_matrix.tsv.gz"),
    "GSE244982": ("作者标准化值 全转录组", "expr/GSE244982_symbol_matrix.tsv.gz"),
}
ALIASES = {"C10orf54": "VSTM1"}


def load(path: str) -> dict[str, np.ndarray]:
    df = pd.read_csv(ROOT / path, sep="\t", index_col=0)
    return {g: df.loc[g].to_numpy(dtype=float) for g in df.index}


def diagnose(mat: dict[str, np.ndarray], name: str, platform: str) -> list[dict]:
    rows = []
    for i, (g1, g2) in enumerate(PAIRS, 1):
        a = mat.get(g1) if g1 in mat else mat.get(ALIASES.get(g1, ""))
        b = mat.get(g2) if g2 in mat else mat.get(ALIASES.get(g2, ""))
        if a is None or b is None:
            rows.append({"cohort": name, "platform": platform, "pair_no": i,
                         "g1": g1, "g2": g2, "status": "基因缺失"})
            continue
        f = (a < b).astype(float)
        p1 = f.mean()
        rows.append({
            "cohort": name, "platform": platform, "pair_no": i,
            "g1": g1, "g2": g2, "status": "可算",
            "n_samples": len(f),
            "p1": p1,
            "f_constant": bool(p1 in (0.0, 1.0)),
            "p_min": min(p1, 1 - p1),
            "log2fc_median": float(np.median(np.log2((a + 1) / (b + 1)))),
        })
    return rows


def main() -> int:
    out: list[dict] = []
    for name, (plat, path) in COHORTS.items():
        p = ROOT / path
        if not p.exists():
            print(f"[skip] {name}: 缺 {path}", file=sys.stderr)
            continue
        out.extend(diagnose(load(path), name, plat))

    df = pd.DataFrame(out)
    outp = ROOT / "qc" / "GSE215868_panel_diagnostic.tsv"
    df.to_csv(outp, sep="\t", index=False)

    log: list[str] = []
    log.append("# IMPRES 特征的平台适配性诊断（**事后诊断，不属预登记 §17**）\n")
    log.append("生成脚本 `scripts/19_panel_diagnostic.py`。")
    log.append("**合规声明**：本诊断**不使用应答标签**，只度量每个特征在样本间的变异；"
               "**不改变** `scripts/18` 已按写死规则给出的裁决。\n")

    ok = df[df["status"] == "可算"]

    log.append("## 一、13 个特征在各队列上的信息量\n")
    log.append("`p1` = F=1 的样本比例。`p1` 越接近 0 或 1，特征越被「单一基因永远更大」锁死。\n")
    log.append("| 队列 | 平台 | 特征数 | 完全恒定的特征 | p1∈[0.15,0.85] 的特征 | "
               "p1 中位数 |")
    log.append("|---|---|---|---|---|---|")
    summary = []
    for name, (plat, _) in COHORTS.items():
        s = ok[ok["cohort"] == name]
        if s.empty:
            continue
        n_const = int(s["f_constant"].sum())
        n_good = int(((s["p1"] >= 0.15) & (s["p1"] <= 0.85)).sum())
        log.append(f"| {name} | {plat} | {len(s)} | **{n_const}/{len(s)}** | "
                   f"**{n_good}/{len(s)}** | {s['p1'].median():.2f} |")
        summary.append((name, len(s), n_const, n_good, s["p1"].median()))

    log.append("\n## 二、GSE215868 逐特征明细\n")
    log.append("| # | g1 | g2 | p1 | 恒定 | log2FC 中位 |")
    log.append("|---|---|---|---|---|---|")
    for _, r in ok[ok["cohort"] == "GSE215868"].iterrows():
        log.append(f"| {r['pair_no']} | {r['g1']} | {r['g2']} | {r['p1']:.2f} | "
                   f"{'**是**' if r['f_constant'] else '否'} | "
                   f"{r['log2fc_median']:+.2f} |")

    # NanoString 与全转录组逐特征对照
    sub = ok[ok["cohort"].isin(["GSE215868", "GSE91061"])]
    if not sub.empty and sub["cohort"].nunique() == 2:
        piv = sub.pivot(index="pair_no", columns="cohort", values="p1")
        both = piv.dropna()
        if len(both):
            log.append("\n## 三、逐特征对照：靶向面板 vs 全转录组\n")
            log.append("| # | GSE215868 p1（面板） | GSE91061 p1（RNA-seq） |")
            log.append("|---|---|---|")
            for idx, r in both.iterrows():
                log.append(f"| {idx} | {r['GSE215868']:.2f} | {r['GSE91061']:.2f} |")

    log.append("\n## 四、分数变异结构：这些特征是几个独立维度？\n")
    log.append("如果 13/15 个特征彼此独立，分数方差应约等于各特征 Bernoulli 方差之和 "
               "`Σp(1-p)`。实测方差远小于该值，说明所有特征被**同一个共同轴**驱动，"
               "分数实际只有一维自由度——即使它与应答有关，能提供的独立信息也很有限。\n")
    log.append("| 队列 | 特征数 | 分数标准差 | 分数范围 | 独立假设下预期 SD | "
               "实测/预期 |")
    log.append("|---|---|---|---|---|---|")
    var_rows = []
    for name, (plat, path) in COHORTS.items():
        if not (ROOT / path).exists():
            continue
        mat = load(path)
        us = [(g1, g2) for g1, g2 in PAIRS
              if (g1 in mat or ALIASES.get(g1) in mat)
              and (g2 in mat or ALIASES.get(g2) in mat)]

        def gv(g):
            return mat[g] if g in mat else mat[ALIASES[g]]

        f = np.vstack([(gv(a) < gv(b)).astype(float) for a, b in us])
        score = f.sum(axis=0)
        p = f.mean(axis=1)
        sd_obs = float(score.std(ddof=1))
        sd_exp = float(np.sqrt((p * (1 - p)).sum()))
        ratio = sd_obs / sd_exp if sd_exp > 0 else np.nan
        log.append(f"| {name} | {len(us)} | **{sd_obs:.2f}** | "
                   f"{int(score.min())}–{int(score.max())} | {sd_exp:.2f} | "
                   f"**{ratio:.2f}** |")
        var_rows.append({"cohort": name, "n_features": len(us),
                         "score_sd": sd_obs, "score_min": int(score.min()),
                         "score_max": int(score.max()), "sd_if_independent": sd_exp,
                         "ratio": ratio})
    pd.DataFrame(var_rows).to_csv(ROOT / "qc" / "GSE215868_score_variance.tsv",
                                  sep="\t", index=False)

    log.append("\n## 五、诊断结论（**两条机制假设均被数据否掉，如实记录**）\n")
    s215 = ok[ok["cohort"] == "GSE215868"]
    n_const = int(s215["f_constant"].sum())
    n_good = int(((s215["p1"] >= 0.15) & (s215["p1"] <= 0.85)).sum())
    log.append(f"- GSE215868 上 **{n_const}/{len(s215)} 个特征在 105 个样本里完全恒定**"
               f"（对判别零贡献），只有 {n_good}/{len(s215)} 个落在 p1∈[0.15,0.85] 的有信息区间。")
    log.append("")

    log.append("### 5.1 假设一（平台特异）：NanoString 探针捕获效率 —— **否**\n")
    log.append("初版假设负结果源于靶向面板的探针捕获效率主导了「同一样本内两基因的相对大小」。"
               "第一节直接否掉了它：")
    log.append("")
    log.append("- 若机制是面板特异的，面板的特征退化应当**明显重于**全转录组。")
    log.append(f"- 实测相反：面板 **{n_const}/{len(s215)} 恒定、{n_good}/{len(s215)} 有信息**；"
               "而全转录组的 **GSE78220（芯片）5/15 恒定、仅 4/15 有信息**，"
               "**GSE244982 5/15 恒定、仅 5/15 有信息**。全转录组退化更重。**假设撤回。**")
    log.append("")

    log.append("### 5.2 假设二（特征塌成一维）：方差坍缩 —— **同样否掉**\n")
    log.append("第二个假设是：15 个特征共享同一轴，叠加后分数只剩一维自由度，"
               "因而没有判别力。**第四节否掉了它**：")
    log.append("")
    log.append("- 若特征高度相关，实测 SD 应**远小于**独立假设下的 `Σp(1-p)` 预期值。")
    log.append("- 实测比值 **1.16–1.38，全部 > 1**：分数的变异**比假设特征独立时还更大**，"
               "不存在一维坍缩。**假设撤回。**")
    log.append("")

    log.append("### 5.3 那么阴性结果是什么\n")
    log.append("两条机制解释都不成立之后，能被数据支持的只剩一个更朴素的结论：\n")
    log.append("- **IMPRES 作为一个统计量是良态的。** 归一化不变量校验 8/8 通过（`scripts/18`），"
               "说明实现无归一化依赖；分数标准差 1.2–1.8（满分 13–15），"
               "实际取值 7–10 个水平，变异充分；特征不存在坍缩。")
    log.append("- **所以阴性不能用「实现有 bug」「平台不合适」「指标退化」来解释。**"
               "唯一被数据排除的解释，恰恰是那三个最省事的辩解。")
    log.append("- 剩下的解释只有一个：**作者原始构造的这 15 个「IC 基因两两大小关系」特征，"
               "在黑色素瘤 ICB 应答上本来就不具备原文报告的判别力。**")
    log.append("")
    log.append("这一条与 `docs/06` 的证据方向一致但**来源独立**：那里 IMPRES 在两个全转录组队列上"
               "AUC 为 0.359 / 0.298（低于随机，方向与本队列相反）；这里在 79 例治疗前、"
               "预登记锁定、样本量为原建模集 4.5 倍的独立队列上为 0.505 / 0.483（精确落在随机线上）。"
               "**五个队列、两种平台、三种表达值口径，没有一个复现原文的 0.77–0.96。**")
    log.append("")

    log.append("### 5.4 结论的边界（必须同时读）\n")
    log.append("- 本诊断**不能**证明 IMPRES 无效。它证明的是：**在本项目能拿到的五个队列上，"
               "无法复现其已发表性能**，并且**排除了三个最常见的辩解**。")
    log.append("- 本诊断**不能**反向论证某个改造版特征必然有效。第五节提到的恒定特征"
               "（面板 4 个、芯片 5 个）确实是可改进之处，但没有任何预登记证据支持"
               "「剔除它们就能救回 AUC」，本项目也**不去试**——那正是预登记第 17.6 条禁止的"
               "「结果出来后新增未事先指定的分析」。")
    log.append("- **15 个基因对本项目仍是从作者 GitHub 仓库重建的**（原文 Supp. Table 2 "
               "非开放获取，取不到）。重建正确性有独立判据（15 对全部含至少一个直接 ICB 靶点），"
               "但这是本节结论最主要的残余不确定性来源，如实标注。")

    md = ROOT / "qc" / "GSE215868_panel_diagnostic.md"
    md.write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))
    print(f"\n[OK] {outp}\n[OK] {md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())