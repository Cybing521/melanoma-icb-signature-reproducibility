#!/usr/bin/env python3
"""检查 Hallmark 基因集在各队列已建矩阵中的覆盖率。

这是**看到任何相关性/分类结果之前**必须先确认的事：
若某队列的注释体系（Entrez / Ensembl / 芯片探针）与 Hallmark 符号对不上，
后面所有签名打分都是空的，而这件事只有在做特征可及性筛选时才会暴露。
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def read_gmt(p: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    with p.open(encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            out[parts[0]] = {g.strip().upper() for g in parts[2:] if g.strip()}
    return out


def read_symbols(matrix: Path) -> set[str]:
    op = gzip.open if matrix.suffix == ".gz" else open
    with op(matrix, "rt", encoding="utf-8", errors="replace") as fh:
        fh.readline()
        return {ln.split("\t", 1)[0].strip().upper() for ln in fh if ln.strip()}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--gmt", required=True)
    ap.add_argument("--label", default="")
    args = ap.parse_args()

    have = read_symbols(Path(args.matrix))
    sets = read_gmt(Path(args.gmt))
    label = args.label or Path(args.matrix).stem

    print(f"\n=== {label} ===")
    print(f"矩阵基因符号 {len(have)} 个")
    rows = []
    for name, genes in sets.items():
        hit = genes & have
        rows.append((name, len(hit), len(genes), len(hit) / len(genes) * 100))
    rows.sort(key=lambda r: r[3])
    low = [r for r in rows if r[3] < 50]
    print(f"基因集 {len(rows)} 个；覆盖率 <50% 的有 {len(low)} 个")
    if low:
        print("  最低 10 个：")
        for name, hit, tot, pct in rows[:10]:
            print(f"    {name:<46} {hit:>4}/{tot:<4} {pct:6.1f}%")
    else:
        best = rows[0]
        print(f"  最低覆盖率：{best[0]} {best[1]}/{best[2]} = {best[3]:.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
