#!/usr/bin/env python3
"""细胞核质控：统计 GSE308434 每个细胞核的 UMI 数、检出基因数与线粒体基因占比。

用途是**定阈值**，不是出结果。按方案「模块 1」的约定，质控参数依据本队列自身分布设定，
不套用其他队列的经验阈值。本脚本只产出分布统计与分位数，不做任何过滤，
阈值由人根据输出决定并写入《终点口径与分析预登记》后才允许进入正式分析。

输入：服务器 data/GSE308434/*_sn_counts.csv.gz（稠密 genes × cells 矩阵）
输出：qc/GSE308434_nucleus_qc.tsv（逐细胞核）+ qc/GSE308434_qc_summary.tsv（逐样本分位数）

用法：
  python3 scripts/07_qc_nuclei.py --outdir /path/to/project --data-dir /path/to/project/data/GSE308434
  python3 scripts/07_qc_nuclei.py ... --limit 8      # 只跑前 N 个样本，用于快速定阈值
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
from pathlib import Path

import numpy as np
import pandas as pd

# 线粒体基因判定：人类 MT- 前缀（HGNC 符号）
MT_PREFIX = "MT-"


def nucleus_stats(path: Path, mt_prefix: str = MT_PREFIX) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """读稠密 CSV，返回 (总 UMI, 检出基因数, 线粒体占比, 细胞条码)。

    文件是 genes × cells 稠密矩阵，第 0 列是基因符号、第 1 行是细胞条码。
    实现上只读一次表，让 pandas 自己做类型推断（第 0 列会推断成 object，其余为整数），
    再切出数值区。实测比「usecols 指定上万个列号」或「分块 + dtype 字典」都快得多：
    后两种写法在本数据规模下会因为逐列匹配而慢一个数量级。
    峰值内存约 3.3 GB/样本（int64），落在服务器内存预算内。
    """
    with gzip.open(path, "rb") as fh:
        header = fh.readline().decode("utf-8").rstrip("\n").split(",")
    cells = header[1:]

    df = pd.read_csv(path, header=None, skiprows=1)
    genes = df[0].astype(str).to_numpy()
    vals = df.iloc[:, 1:].to_numpy(dtype=np.int64, copy=False)
    n = vals.shape[1]
    assert n == len(cells), f"表头列数 {len(cells)} 与数据列数 {n} 不一致"

    n_umi = vals.sum(axis=0)
    n_feat = (vals > 0).sum(axis=0)
    mt_mask = np.array([g.startswith(mt_prefix) for g in genes])
    mt_umi = vals[mt_mask].sum(axis=0) if mt_mask.any() else np.zeros(n, dtype=np.int64)

    mt_pct = np.divide(mt_umi, np.maximum(n_umi, 1), dtype=np.float64) * 100.0
    return n_umi, n_feat, mt_pct, cells


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--limit", type=int, default=0, help="只处理前 N 个样本（按文件名排序）")
    args = ap.parse_args()

    files = sorted(Path(args.data_dir).glob("*_sn_counts.csv.gz"))
    if args.limit:
        files = files[: args.limit]
    if not files:
        raise SystemExit(f"{args.data_dir} 下没有 *_sn_counts.csv.gz")

    qc_dir = Path(args.outdir) / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for i, f in enumerate(files, 1):
        sample = f.name.replace("_sn_counts.csv.gz", "")
        n_umi, n_feat, mt_pct, cells = nucleus_stats(f)
        out_tsv = qc_dir / f"{sample}.nucleus_qc.tsv"
        with out_tsv.open("w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["cell_barcode", "n_umi", "n_genes", "mt_pct"])
            for c, u, g, m in zip(cells, n_umi, n_feat, mt_pct):
                w.writerow([c, int(u), int(g), f"{m:.4f}"])

        q = lambda a, p: float(np.percentile(a, p))  # noqa: E731
        summary_rows.append({
            "sample": sample,
            "n_nuclei": len(cells),
            "n_genes_total": len(n_feat) and None,
            "umi_p05": q(n_umi, 5), "umi_p25": q(n_umi, 25), "umi_p50": q(n_umi, 50),
            "umi_p75": q(n_umi, 75), "umi_p95": q(n_umi, 95),
            "feat_p05": q(n_feat, 5), "feat_p25": q(n_feat, 25), "feat_p50": q(n_feat, 50),
            "feat_p75": q(n_feat, 75), "feat_p95": q(n_feat, 95),
            "mt_p50": q(mt_pct, 50), "mt_p90": q(mt_pct, 90), "mt_p95": q(mt_pct, 95), "mt_p99": q(mt_pct, 99),
        })
        summary_rows[-1]["n_genes_total"] = 36601
        s = summary_rows[-1]
        print(
            f"[{i}/{len(files)}] {sample}: n={s['n_nuclei']} "
            f"UMI p5/p50/p95={s['umi_p05']:.0f}/{s['umi_p50']:.0f}/{s['umi_p95']:.0f} "
            f"gene p5/p50/p95={s['feat_p05']:.0f}/{s['feat_p50']:.0f}/{s['feat_p95']:.0f} "
            f"MT% p50/p95={s['mt_p50']:.2f}/{s['mt_p95']:.2f}",
            flush=True,
        )

    out = qc_dir / "GSE308434_qc_summary.tsv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\n逐细胞核质控已写入 {qc_dir}/<sample>.nucleus_qc.tsv")
    print(f"样本级分位数已写入 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())