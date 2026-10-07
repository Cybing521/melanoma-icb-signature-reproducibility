#!/usr/bin/env python3
"""模块 2 免疫侧：从 GSE308435 的克隆型表计算每份活检的克隆性指标。

每份活检对应一个 TCR 文库，与 GSE308434 的细胞核文库一一对应（同 pre/on/post 命名）。
本脚本只做汇总统计，不做任何跨样本比较或结论推断。

输入：data_ext/GSE308435/*_tcr_clonotypes.csv.gz
      列：clonotype_id, frequency, proportion, cdr3s_aa, cdr3s_nt, inkt_evidence, mait_evidence
输出：qc/GSE308435_tcr_clonality.tsv

指标定义（全部写在这里，便于复核）：
  n_clonotypes        该活检中检出的克隆型总数
  total_cells         frequency 求和，等于该文库中参与比对的 T 细胞总数
  top1_frac           最大克隆型频率占比，克隆扩增的直接指标
  top10_frac          前 10 个克隆型频率占比之和
  clonality           1 - Shannon 熵 / ln(n_clonotypes)，0 = 完全均匀（无克隆性），1 = 全部集中在一个克隆
  clones_for_50pct    累计频率达到 50% 所需的克隆型个数
  n_inkt              带 iNKT 证据的克隆型数
  n_mait              带 MAIT 证据的克隆型数

用法：
  python3 scripts/08_tcr_clonality.py --outdir /path/to/project --tcr-dir /path/to/project/data_ext/GSE308435
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from pathlib import Path


def metrics(path: Path) -> dict:
    freq: list[float] = []
    n_inkt = n_mait = 0
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            f = (row.get("frequency") or "").strip()
            if f:
                freq.append(float(f))
            if (row.get("inkt_evidence") or "").strip():
                n_inkt += 1
            if (row.get("mait_evidence") or "").strip():
                n_mait += 1

    freq.sort(reverse=True)
    total = sum(freq)
    n = len(freq)
    if n == 0 or total == 0:
        return {
            "n_clonotypes": 0, "total_cells": 0, "top1_frac": 0.0, "top10_frac": 0.0,
            "clonality": 0.0, "clones_for_50pct": 0, "n_inkt": n_inkt, "n_mait": n_mait,
        }

    p = [f / total for f in freq]
    entropy = -sum(x * math.log(x) for x in p if x > 0)
    clonality = 0.0 if n <= 1 else 1.0 - entropy / math.log(n)

    cum = 0.0
    need = 0
    for x in p:
        cum += x
        need += 1
        if cum >= 0.5:
            break

    return {
        "n_clonotypes": n,
        "total_cells": int(total),
        "top1_frac": round(p[0], 6),
        "top10_frac": round(sum(p[:10]), 6),
        "clonality": round(clonality, 6),
        "clones_for_50pct": need,
        "n_inkt": n_inkt,
        "n_mait": n_mait,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tcr-dir", required=True)
    args = ap.parse_args()

    files = sorted(Path(args.tcr_dir).glob("*_tcr_clonotypes.csv.gz"))
    if not files:
        raise SystemExit(f"{args.tcr_dir} 下没有 *_tcr_clonotypes.csv.gz")

    rows = []
    for f in files:
        # GSM9245440_F01_on_tcr_clonotypes.csv.gz
        #
        # 两个系列（snRNA GSE308434 / TCR GSE308435）的 GSM 号完全不同
        # （F01_on 在前者是 GSM9245398、在后者是 GSM9245440），
        # 唯一可靠的配对键是去掉 GSM 前缀后的**文库代号**。
        # 因此同时输出 sample（完整，带 GSM，便于溯源）与 library（配对用）。
        sample = f.name[: -len("_tcr_clonotypes.csv.gz")]
        library = re.sub(r"^GSM\d+_", "", sample)
        m = metrics(f)
        rows.append({"sample": sample, "library": library, **m})
        print(
            f"{sample} (library={library}): clones={m['n_clonotypes']} cells={m['total_cells']} "
            f"top1={m['top1_frac']:.4f} clonality={m['clonality']:.4f} "
            f"need50%={m['clones_for_50pct']}",
            flush=True,
        )

    qc = Path(args.outdir) / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    out = qc / "GSE308435_tcr_clonality.tsv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"\n已写入 {out}（{len(rows)} 个文库）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())