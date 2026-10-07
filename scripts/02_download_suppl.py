#!/usr/bin/env python3
"""下载已通过标签核对的 bulk 队列补充文件（表达矩阵），并记录清单与校验值。

只处理 characteristics 中已确认应答/分组字段的队列：
  GSE91061  主开发队列（response + visit，109 样本）
  GSE78220  外部验证队列（anti-pd-1 response，28 样本）
  GSE294272 参考谱来源（treatment 三组，29 样本）
  GSE244982 耐药队列（groups，41 样本，全耐药）

GSE308434 未确认应答标签，按数据门禁结论暂不下载。

用法：
  python3 scripts/02_download_suppl.py --outdir /path/to/project            # dry-run，只列计划
  python3 scripts/02_download_suppl.py --outdir /path/to/project --execute  # 实际下载
"""

from __future__ import annotations

import argparse
import hashlib
import re
import time
import urllib.request
from pathlib import Path

TARGETS = {
    "GSE91061": "主开发队列：response + visit，109 样本",
    "GSE78220": "外部验证队列：anti-pd-1 response，28 样本",
    "GSE294272": "细胞状态参考谱：treatment 三组，29 样本",
    "GSE244982": "耐药队列：groups，41 样本，无应答者",
}
DEFER = {
    "GSE308434": "GEO 无应答字段，按数据门禁结论暂缓",
}

ROW_RE = re.compile(
    r'<a href="(?P<name>[^"]+)">[^<]+</a>\s+(?P<date>[\d-]+\s+[\d:]+)\s+(?P<size>[\d.]+\s*[KMG]?)'
)
SIZE_RE = re.compile(r"([\d.]+)\s*([KMG]?)", re.I)
MULTIPLIER = {"": 1, "K": 1024, "M": 1024**2, "G": 1024**3}


def supp_url(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/suppl/"


def to_bytes(text: str) -> int:
    match = SIZE_RE.search(text)
    if not match:
        return 0
    return int(float(match.group(1)) * MULTIPLIER[match.group(2).upper()])


def list_files(gse: str, retries: int = 3) -> list[tuple[str, int]]:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(supp_url(gse), headers={"User-Agent": "geo-suppl-fetch/1.0"})
            with urllib.request.urlopen(req, timeout=90) as resp:
                html = resp.read().decode("utf-8", errors="replace")
            files = [
                (m.group("name"), to_bytes(m.group("size")))
                for m in ROW_RE.finditer(html)
                if not m.group("name").endswith("/")
            ]
            if files:
                return files
        except Exception as exc:  # noqa: BLE001
            last = exc
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"无法列出 {gse} 补充文件: {last}")


def download(url: str, dest: Path, retries: int = 3) -> None:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "geo-suppl-fetch/1.0"})
            with urllib.request.urlopen(req, timeout=600) as resp, dest.open("wb") as fh:
                while chunk := resp.read(1 << 20):
                    fh.write(chunk)
            return
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"下载失败 {url}: {last}")


def md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--execute", action="store_true", help="实际下载；缺省只打印计划")
    ap.add_argument("--max-total-gb", type=float, default=8.0, help="总量上限，超过则中止")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)

    plan: list[tuple[str, str, int]] = []
    for gse, role in TARGETS.items():
        files = list_files(gse)
        total = sum(size for _, size in files)
        print(f"\n{gse}  {role}")
        for name, size in files:
            print(f"    {size / 1024**2:9.2f} MB  {name}")
            plan.append((gse, name, size))
        print(f"    {'':>9}       小计 {total / 1024**2:.2f} MB")

    for gse, reason in DEFER.items():
        print(f"\n{gse}  暂缓：{reason}")

    grand = sum(size for _, _, size in plan)
    print(f"\n合计 {grand / 1024**3:.3f} GB，共 {len(plan)} 个文件")
    if not args.execute:
        print("dry-run 结束，加 --execute 执行下载")
        return 0
    if grand / 1024**3 > args.max_total_gb:
        print(f"超过体积上限 {args.max_total_gb} GB，已中止")
        return 1

    manifest = ["gse\tfile\tbytes\tmd5"]
    for gse, name, size in plan:
        dest_dir = outdir / "data" / gse
        dest_dir.mkdir(parents=True, exist_ok=True)
        dest = dest_dir / name
        url = f"{supp_url(gse)}{name}"
        if dest.exists() and dest.stat().st_size == size:
            print(f"[跳过] {gse}/{name}")
        else:
            print(f"[下载] {gse}/{name} ({size / 1024**2:.2f} MB)", flush=True)
            download(url, dest)
        actual = dest.stat().st_size
        digest = md5(dest)
        manifest.append(f"{gse}\t{name}\t{actual}\t{digest}")
        print(f"[完成] {actual} bytes  md5={digest}")

    out = outdir / "meta" / "suppl_manifest.tsv"
    out.write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(f"\n清单已写入 {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())