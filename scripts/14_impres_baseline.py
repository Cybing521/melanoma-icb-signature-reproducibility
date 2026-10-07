#!/usr/bin/env python3
"""IMPRES 基线（Auslander et al., Nat Med 2018, PMID 30127394）。

15 个特征对的来源（不凭记忆写，全部回溯到作者仓库）
--------------------------------------------------
论文正文只给出构造方法，具体的 15 对基因在 **Supp. Table 2**，而该文非开放获取
（Europe PMC `supplementaryFiles` 返回 "not open access"；PMC 网页直下触发
JS cookie 门，Nature/Springer 静态资源 403）。因此改走**作者 GitHub 仓库**
`https://github.com/noamaus/IMPRES-codes`：

* `ADDITIONAL_CODES/Additional feature sets/FEATS.mat`  → 15 个特征索引（uint16，1-based）
  `[31, 61, 128, 169, ...]`
* `ADDITIONAL_CODES/Additional feature sets/CPall.mat`  → 28 个 IC 基因（字母序）
* 索引 → 基因对的映射由 `getIMPRESRAT.m` / `classifyImmuneCOMP.m` 的双层
  for 循环枚举顺序决定，见 `scripts/impres_pairs.py`

**映射正确性的独立验证**：论文写明候选对必须"至少含一个基因属于 6 个直接
ICB 靶点（CTLA4、CD28、CD80、CD86、PD-1、PD-L1）"。还原出的 15 对**全部满足**
这一条，若索引偏移一位即会破坏该性质。

评分公式（论文 Methods）
------------------------
    F_{i,j}(x) = 1  若 exp_i(x) < exp_j(x)，否则 0
    IMPRES(x) = Σ_{k=1..15} F(x)，取值 0–15，越高越预测应答

输入
----
expr/<cohort>_symbol_matrix.tsv.gz（行 = 基因符号）
输出
----
qc/IMPRES_scores_<cohort>.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from impres_pairs import FEATURE_IDX, IC_GENES, idx_to_pair  # noqa: E402

PAIRS = [idx_to_pair(i) for i in FEATURE_IDX]

# HGNC 官方改名后的别名。只处理已知的、被证实的改名，不是"找不到就随便换一个"：
#   C10orf54 → VSTM1（HGNC 已更名）
# 各队列注释体系不同，GSE91061 与 GSE244982 用旧名 GSE78220 / GSE294272 用新名。
SYMBOL_ALIASES = {"C10orf54": "VSTM1"}


def resolve(sym: str, have: set[str]) -> str | None:
    """返回该基因在当前矩阵中的实际符号；不存在返回 None。"""
    if sym in have:
        return sym
    alt = SYMBOL_ALIASES.get(sym)
    return alt if alt and alt in have else None


def read_matrix(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    op = gzip.open if path.suffix == ".gz" else open
    expr: dict[str, list[float]] = {}
    samples: list[str] = []
    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        cols = next(csv.reader(fh, delimiter="\t"))
        samples = cols[1:]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            expr[parts[0]] = [float(v) for v in parts[1:]]
    return samples, expr


def quantile_normalise(mat: list[list[float]]) -> list[list[float]]:
    """按行（基因）在样本维度做分位数归一化——论文用 quantile-normalized expression。"""
    n_samp = len(mat[0])
    order_by_row: list[list[int]] = []
    vals_by_row: list[list[float]] = []
    for row in mat:
        order = sorted(range(n_samp), key=lambda k: row[k])
        order_by_row.append(order)
        vals_by_row.append(sorted(row))

    # 目标分布：所有基因排序值的均值（逐列取均值）
    target = [sum(vals_by_row[g][k] for g in range(len(mat))) / len(mat)
              for k in range(n_samp)]

    out = [[0.0] * n_samp for _ in mat]
    for g, order in enumerate(order_by_row):
        for rank, s in enumerate(order):
            out[g][s] = target[rank]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--no-qnorm", action="store_true",
                    help="跳过分位数归一化（仅作敏感性分析，默认执行）")
    ap.add_argument("--orientation", choices=["g1_low", "g2_low"], default="g1_low",
                    help="F=1 的方向。默认 g1_low：F=1 当 exp_{g1} < exp_{g2}，"
                         "与 getIMPRESRAT.m / classifyImmuneCOMP.m 的 ratio=exp_g1/exp_g2 一致。"
                         "g2_low 为镜像方向，仅作敏感性分析（见文件末尾说明）")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    mtx = outdir / "expr" / f"{args.cohort}_symbol_matrix.tsv.gz"
    if not mtx.exists():
        print(f"缺表达矩阵 {mtx}", file=sys.stderr)
        return 1

    samples, expr = read_matrix(mtx)
    have = set(expr)

    resolved: dict[str, str] = {}
    missing: list[str] = []
    for p in PAIRS:
        for g in p:
            if g in resolved:
                continue
            r = resolve(g, have)
            if r is None:
                missing.append(g)
            else:
                if r != g:
                    print(f"  别名映射：{g} → {r}（HGNC 官方改名）")
                resolved[g] = r
    if missing:
        print(f"矩阵内缺失基因 {sorted(set(missing))}，无法计算 IMPRES", file=sys.stderr)
        return 1
    print(f"{args.cohort}: {len(samples)} 样本 × {len(expr)} 基因，"
          f"15 对涉及 {len(resolved)} 个基因全部命中")

    genes = sorted(set(resolved.values()))
    mat = []
    for g in genes:
        row = expr[g]
        # 该队列的量纲不同（raw 计数 / FPKM / 作者标准化值），统一取 log2 后再归一化
        mat.append([math.log2(v + 1.0) for v in row])
    if not args.no_qnorm:
        mat = quantile_normalise(mat)

    gpos = {g: i for i, g in enumerate(genes)}
    scores = []
    for s in range(len(samples)):
        tot = 0
        for g1, g2 in PAIRS:
            a, b = resolved[g1], resolved[g2]
            lo, hi = (a, b) if args.orientation == "g1_low" else (b, a)
            tot += 1 if mat[gpos[lo]][s] < mat[gpos[hi]][s] else 0
        scores.append(tot)

    tag = "" if args.orientation == "g1_low" else "_mirror"
    out = outdir / "qc" / f"IMPRES_scores_{args.cohort}{tag}.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "IMPRES", "orientation"])
        for s, sc in zip(samples, scores):
            w.writerow([s, sc, args.orientation])
    print(f"\n15 对：{PAIRS}")
    dist = {v: scores.count(v) for v in sorted(set(scores))}
    print(f"IMPRES 取值分布（0–15）：{dist}")
    print(f"已写入 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


"""
⚠️ 方向约定的诚实说明
--------------------
作者仓库里两个 .m 文件本身就不一致：
  * `getIMPRESRAT.m`        区块 1 用 GE(c(i))/GE(c(j))，区块 2 用 GE(c(j))/GE(c(i))
  * `classifyImmuneCOMP.m`  区块 1 用 GE(c(j))/GE(c(i))，区块 2 用 GE(c(i))/GE(c(j))
两者的比值方向在区块间是反的。论文正文只写 F_{i,j}(x)=1 当 exp_i < exp_j，
没说 Supp. Table 2 里每对的书写方向。

因此本脚本**默认 `--orientation g1_low`**（与 ratio=exp_{g1}/exp_{g2} 一致），
并提供 `--orientation g2_low`（镜像）作为敏感性分析。

**两个方向都报，不挑好看的那个当主结果。** 若某方向与本队列判别方向相反，
只能说明二者之一与作者实际用法不符，不能据此反推"哪个才对"——
那需要作者原始 Supp. Table 2，而该文件不可获取（见文件头说明）。
"""
