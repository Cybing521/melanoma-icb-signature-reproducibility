#!/usr/bin/env python3
"""逐队列打印 characteristics 键的取值分布，并检查列数不齐的样本行。

用于数据门禁第二步：确认响应状态、访视时点、瘤种等关键字段是否可直接使用，
还是必须回原文献补充表。

用法：
  python3 scripts/01_inspect_metadata.py --outdir /path/to/project [GSE ...]
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def load(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    # GEO 元数据取值含引号，禁用 csv 引号解析，避免表头被当作引号字段拆坏
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
        cols = list(reader.fieldnames or [])
        rows = [{(k or ""): (v if v is not None else "") for k, v in r.items()} for r in reader]
    return cols, rows


def describe(series: str, path: Path, max_values: int = 25) -> None:
    cols, rows = load(path)
    print(f"\n{'=' * 72}\n{series}  n={len(rows)}\n{'=' * 72}")

    key_cols = [c for c in cols if c.startswith("char::")]
    for col in key_cols:
        name = col[6:]
        values = Counter(r.get(col, "") for r in rows)
        blank = values.pop("", 0)
        total = sum(values.values()) + blank
        print(f"\n-- {name}  (非空 {total - blank}/{total}, 取值 {len(values)} 种)")
        for value, count in values.most_common(max_values):
            print(f"     {value[:70]:<70} {count}")

    # 访视时点与患者标识常藏在样本标题里，单独统计
    titles = [r.get("title", "") for r in rows]
    if titles:
        patients = {t.split("_")[0] for t in titles if t}
        print(f"\n-- 标题推断: 唯一患者/个体前缀 {len(patients)} 个; 例: {sorted(patients)[:8]}")
        print(f"-- 标题例: {titles[:4]}")

    ragged = [r.get("gsm", "?") for r in rows if not r.get("gsm") or not r.get("title")]
    if ragged:
        print(f"\n-- !! 列不齐或缺少 GSM/标题的行: {ragged}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("series", nargs="*", help="只检查指定 GSE，默认全部")
    args = ap.parse_args()

    meta = Path(args.outdir) / "meta"
    targets = args.series or sorted(p.stem.replace("_samples", "") for p in meta.glob("*_samples.tsv"))
    for gse in targets:
        path = meta / f"{gse}_samples.tsv"
        if not path.exists():
            print(f"[缺失] {gse}")
            continue
        describe(gse, path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())