#!/usr/bin/env python3
"""模块 2 配对分析：把肿瘤内在状态分与同一份活检的 TCR 克隆性指标对上。

配对键
------
GSE308434（snRNA）、GSE308435（TCR）、GSE308433（lpWGS）三个系列**连 GSM 号都不同**
（同一个 `F01_on` 在 434 里是 GSM9245398、在 435 里是 GSM9245440、在 433 里是 GSM9245359）。
唯一跨系列共享的是去掉 `GSM<数字>_` 前缀后的**文库代号**，因此一律按 library 后缀配对。

⚠️ 曾踩的坑：`scripts/08_tcr_clonality.py` 早期版本只取 GSM + 访视，产出
`GSM9245440_on`，把文库代号整个丢掉；本脚本按整串取交集，第一步就会报
"样本名没有交集"。现已改为同时输出 `sample` 与 `library` 两列。

统计口径（事先约定，见《终点口径与分析预登记》第十一节）
------------------------------------------------------
* TCR 侧过滤：仅保留 ≥200 个参与比对 T 细胞的文库（阈值已锁定）→ 预计 18 个文库。
* 统计量：Spearman 秩相关（样本量小、指标非正态时比 Pearson 稳）。
* 状态侧统计量用 `*_sd`（样本内状态分离散度）。**不能用 `*_mean`**：状态分由样本内
  基因 z-score 求均值得到，其样本均值在数学上恒等于 0，恒零量没有信息。
* **主比较 3 对**（事先指定，不做挑选）：
    1. ifng_response × clonality
    2. antigen_presentation × clonality
    3. emt × clonality
  其余全部为探索性，用 Benjamini–Hochberg 控制 FDR，且不得用于支持主要结论。
* **不做的事**：不因相关性不显著而更换指标、阈值或改用配对样本子集重跑。

输入：qc/GSE308434_state_scores.tsv、qc/GSE308435_tcr_clonality.tsv
输出：qc/module2_pairing.tsv（配对明细）、qc/module2_correlations.tsv（相关性）
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

GSM_PREFIX = re.compile(r"^GSM\d+_")
MIN_TCR_CELLS = 200
PRIMARY_PAIRS = [
    ("ifng_response", "clonality"),
    ("antigen_presentation", "clonality"),
    ("emt", "clonality"),
]
CLONALITY_METRICS = ["clonality", "top1_frac", "top10_frac", "clones_for_50pct"]
STATE_METRICS = [
    "ifng_response", "antigen_presentation", "proliferation", "emt",
    "hypoxia", "tnf_nfkb",
]


def read_tsv(p: Path) -> list[dict]:
    with p.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def rankdata(xs: list[float]) -> list[float]:
    """平均秩，处理并列。"""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(x: list[float], y: list[float]) -> tuple[float, float]:
    rx, ry = rankdata(x), rankdata(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    if dx == 0 or dy == 0:
        return float("nan"), float("nan")
    rho = num / (dx * dy)
    # Fisher z 变换的近似 p 值（n 小，仅作排序用，不当作精确推断）
    z = 0.5 * math.log((1 + rho) / (1 - rho)) * math.sqrt(max(n - 3, 1))
    p = math.erfc(abs(z) / math.sqrt(2))
    return rho, p


def bh(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    prev = 1.0
    for rank, idx in enumerate(reversed(order), 1):
        k = m - rank + 1
        val = min(prev, pvals[idx] * m / k)
        adj[idx] = val
        prev = val
    return adj


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    qc = Path(args.outdir) / "qc"

    states: dict[str, dict] = {}
    for r in read_tsv(qc / "GSE308434_state_scores.tsv"):
        states[GSM_PREFIX.sub("", r["sample"])] = r

    tcr: dict[str, dict] = {}
    for r in read_tsv(qc / "GSE308435_tcr_clonality.tsv"):
        # 优先用脚本 08 已写好的 library 列；没有就现场剥 GSM 前缀
        key = r.get("library") or GSM_PREFIX.sub("", r["sample"])
        tcr[key] = r

    common = sorted(set(states) & set(tcr))
    if not common:
        raise SystemExit(
            f"文库代号没有交集（状态分 {len(states)} 个、TCR {len(tcr)} 个），配对键可能用错了"
        )

    paired = []
    dropped = []
    for s in common:
        cells = int(tcr[s]["total_cells"])
        if cells < MIN_TCR_CELLS:
            dropped.append((s, cells))
            continue
        paired.append(s)

    print(f"共同样本名 {len(common)} 个；TCR ≥{MIN_TCR_CELLS} 细胞后剩 {len(paired)} 个；剔除 {len(dropped)} 个")
    for s, c in sorted(dropped, key=lambda x: x[1]):
        print(f"  剔除 {s}: 仅 {c} 个 T 细胞")

    with (qc / "module2_pairing.tsv").open("w", encoding="utf-8", newline="") as fh:
        cols = ["library", "tcr_cells", "tcr_clones", "tcr_clonality", "n_nuclei_qc", "qc_kept_pct"] + \
               [f"state_{m}_sd" for m in STATE_METRICS]
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)  # csv.writer 没有 writeheader，必须自己写表头行
        for s in paired:
            st, tc = states[s], tcr[s]
            w.writerow([s, tc["total_cells"], tc["n_clonotypes"], tc["clonality"],
                        st["n_nuclei_qc"], st["qc_kept_pct"]]
                       + [st.get(f"{m}_sd", "") for m in STATE_METRICS])

    results = []
    for sm in STATE_METRICS:
        for cm in CLONALITY_METRICS:
            xs, ys = [], []
            for s in paired:
                xv, yv = states[s].get(f"{sm}_sd", ""), tcr[s].get(cm, "")
                if xv in ("", "nan") or yv in ("", "nan"):
                    continue
                xs.append(float(xv)); ys.append(float(yv))
            if len(xs) < 5:
                continue
            rho, p = spearman(xs, ys)
            # 坑：PRIMARY_PAIRS 里是元组，这里必须用元组查。
            # 用 [sm, cm]（列表）去比元组列表恒为 False，会把事先指定的三对主比较
            # 全部误标成探索性、还混进探索性的 BH-FDR 家族里一起校正。
            primary = (sm, cm) in PRIMARY_PAIRS
            results.append({
                "state": sm, "tcr_metric": cm, "n": len(xs),
                "spearman_rho": round(rho, 4), "p_approx": f"{p:.4g}",
                "is_primary": "yes" if primary else "exploratory",
            })

    # 只对探索性部分做 BH，探索性排在主比较之后不参与校正
    expl_idx = [i for i, r in enumerate(results) if r["is_primary"] == "exploratory"]
    if expl_idx:
        adj = bh([float(results[i]["p_approx"]) for i in expl_idx])
        for i, a in zip(expl_idx, adj):
            results[i]["p_fdr_bh"] = f"{a:.4g}"
    for r in results:
        r.setdefault("p_fdr_bh", "")

    order = {"yes": 0, "exploratory": 1}
    results.sort(key=lambda r: (order[r["is_primary"]], float(r["p_approx"])))

    with (qc / "module2_correlations.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["state", "tcr_metric", "n", "spearman_rho",
                                           "p_approx", "is_primary", "p_fdr_bh"], delimiter="\t")
        w.writeheader()
        w.writerows(results)

    print(f"\n配对明细已写入 {qc / 'module2_pairing.tsv'}")
    print(f"相关性结果已写入 {qc / 'module2_correlations.tsv'}\n")
    print("主比较（事先指定）：")
    for r in results:
        if r["is_primary"] == "yes":
            print(f"  {r['state']:<22}× {r['tcr_metric']:<18} n={r['n']:<3} "
                  f"rho={r['spearman_rho']:>7} p≈{r['p_approx']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())