#!/usr/bin/env python3
"""下载 Izar 组三姊妹系列（GSE308433 / GSE308434 / GSE308435）与 IMPRES 开发队列 GSE115821。

为什么不整包下载 GSE308434_RAW.tar：
  该系列按样本另存了 5 类文件，其中 `sn_matrix.mtx.gz`（稀疏矩阵，单样本 130–180 MB）
  与 `sn_raw_feature_bc_matrix.h5`（未过滤原始矩阵，50–75 MB）占了大头，整包 10 GB。
  本方案需要的是细胞状态打分所需的**基因 × 细胞计数**，GEO 另存了等价的稠密
  `sn_counts.csv.gz`（单样本约 34 MB）+ `sn_features.tsv.gz`（288 KB）。
  42 个样本按此取法约 1.5 GB，且避免了解压时的第二份拷贝（解压会再占 10 GB，
  而 /root/autodl-tmp 只剩约 26 GB）。

  GSE308435（TCR）与 GSE308433（lpWGS）整包很小，直接整包下。

用法：
  python3 scripts/04_download_izar.py --outdir /path/to/project              # dry-run
  python3 scripts/04_download_izar.py --outdir /path/to/project --execute
"""

from __future__ import annotations

import argparse
import hashlib
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

SAMPLE_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/samples/{stub}/{gsm}/suppl/"
SERIES_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/suppl/"

# GSE308434 只取这两类：稠密计数矩阵 + 基因特征表
KEEP_SUFFIX = ("_sn_counts.csv.gz", "_sn_features.tsv.gz")
# 整包下的小系列
FULL_TARS = {
    "GSE308433": ["GSE308433_RAW.tar", "filelist.txt"],
    "GSE308435": ["GSE308435_RAW.tar", "filelist.txt"],
    "GSE115821": ["GSE115821_MGH_counts.csv.gz"],
    # GSE294273（去卷积参考谱来源）只有稀疏 matrix.mtx.gz，按样本取合计约 900 MB，
    # 比整包 742 MB（tar 已压缩）还大，因此直接整包取。
    "GSE294273": ["GSE294273_RAW.tar", "filelist.txt"],
}


def gsm_stub(gsm: str) -> str:
    """GSM9245398 -> GSM9245nnn（GEO 按后三位分桶）"""
    return gsm[:-3] + "nnn"


def series_stub(gse: str) -> str:
    return gse[:-3] + "nnn"


def fetch_text(url: str, retries: int = 3) -> str:
    last: Exception | None = None
    for attempt in range(retries):
        try:
            proc = subprocess.run(
                ["curl", "-sL", "--noproxy", "*", "--max-time", "120", "-A", "geo-suppl-fetch/1.0", url],
                capture_output=True,
                text=True,
            )
            if proc.returncode == 0 and proc.stdout.strip():
                return proc.stdout
        except Exception as exc:  # noqa: BLE001
            last = exc
        time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"取不到 {url}: {last}")


def md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def integrity_ok(path: Path) -> tuple[bool, str]:
    import gzip

    name = path.name.lower()
    try:
        if name.endswith(".gz"):
            with gzip.open(path, "rb") as fh:
                while fh.read(1 << 22):
                    pass
            return True, "gzip 全流解压通过"
        if name.endswith(".tar"):
            proc = subprocess.run(["tar", "tf", str(path)], capture_output=True, text=True, timeout=600)
            n = len([x for x in proc.stdout.splitlines() if x.strip()])
            return proc.returncode == 0 and n > 0, f"tar 可列出 {n} 项"
        return path.stat().st_size > 0, "可读非空"
    except Exception as exc:  # noqa: BLE001
        return False, f"完整性检查失败：{exc}"


def aria2(url: str, dest: Path, conns: int = 8, retries: int = 5) -> int:
    sidecar = dest.with_name(dest.name + ".aria2")
    if dest.exists() and dest.stat().st_size > 0 and not sidecar.exists():
        print(f"[跳过] {dest.name}", flush=True)
        return 0
    if sidecar.exists():
        print(f"[续传] {dest.name}", flush=True)
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "aria2c", "-x", str(conns), "-s", str(conns), "-k", "1M",
        "--file-allocation=none", "--summary-interval=0",
        "--max-tries", str(retries), "--retry-wait=3",
        "-d", str(dest.parent), "-o", dest.name, url,
    ]
    start = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"aria2c 失败 {url}: {proc.stderr[-400:]}")
    if sidecar.exists():
        raise RuntimeError(f"aria2c 退出但控制文件残留，判定未完成: {dest.name}")
    return int(time.time() - start)


def build_plan(gse: str) -> list[tuple[str, str, int]]:
    """返回 [(本地相对路径, url, 期望字节数)]；期望值来自 GEO 的 filelist.txt。"""
    if gse == "GSE308434":
        listing = fetch_text(SERIES_BASE.format(stub=series_stub(gse), gse=gse) + "filelist.txt")
        plan = []
        # 行格式：File<TAB>文件名<TAB>MM/DD/YYYY HH:MM:SS<TAB>字节数<TAB>类型
        row = re.compile(
            r"^File\t(\S+)\t\d{2}/\d{2}/\d{4}\s+\d{2}:\d{2}:\d{2}\t(\d+)\t"
        )
        for line in listing.splitlines():
            m = row.match(line)
            if not m:
                continue
            fname, size = m.group(1), int(m.group(2))
            if not fname.endswith(KEEP_SUFFIX):
                continue
            gsm = fname.split("_", 1)[0]
            plan.append((f"{gse}/{fname}", SAMPLE_BASE.format(stub=gsm_stub(gsm), gsm=gsm) + fname, size))
        if not plan:
            raise RuntimeError(f"{gse} 未从 filelist.txt 解析出任何目标文件")
        return plan

    plan = []
    for name in FULL_TARS.get(gse, []):
        size = 0
        listing = fetch_text(SERIES_BASE.format(stub=series_stub(gse), gse=gse))
        m = re.search(rf'>{re.escape(name)}</a>.*?(\d+(?:\.\d+)?)\s*([KMG])?\s*</td>', listing, re.S)
        if m:
            mult = {"K": 1 << 10, "M": 1 << 20, "G": 1 << 30}.get((m.group(2) or "").upper(), 1)
            size = int(float(m.group(1)) * mult)
        plan.append((f"{gse}/{name}", SERIES_BASE.format(stub=series_stub(gse), gse=gse) + name, size))
    return plan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--workers", type=int, default=4, help="并发下载的文件数（总连接数 = workers × 8）")
    ap.add_argument("--only", default="", help="逗号分隔的队列名，只处理这些队列")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    cohorts = ["GSE308434", "GSE308435", "GSE308433", "GSE115821", "GSE294273"]
    if args.only:
        want = {c.strip() for c in args.only.split(",") if c.strip()}
        cohorts = [c for c in cohorts if c in want]
        if not cohorts:
            raise SystemExit(f"--only {args.only} 未匹配到任何已知队列")

    plan: list[tuple[str, str, int]] = []
    for gse in cohorts:
        print(f"\n=== 规划 {gse} ===")
        for rel, url, size in build_plan(gse):
            print(f"    {size / 1024**2:9.1f} MB  {rel}")
            plan.append((rel, url, size))

    total = sum(s for _, _, s in plan)
    print(f"\n合计 {len(plan)} 个文件 / {total / 1024**3:.2f} GB")
    if not args.execute:
        print("dry-run 结束，加 --execute 执行")
        return 0

    manifest = ["cohort\tfile\tlisted_bytes\tactual_bytes\tmd5\tintegrity\tstatus"]
    failures: list[str] = []

    # 单连接实测约 13 KB/s，aria2 单文件 8 连接约 75 KB/s。GSE308434 有 42 个样本文件，
    # 串行跑完要近 5 小时，因此按文件并发（每个文件内部仍用多连接）。
    # 总并发连接数 = workers × conns，对公开数据仓库保持克制。
    def get_one(item: tuple[str, str, int]) -> tuple[str, str, int, str, str, str, int]:
        rel, url, listed = item
        try:
            secs = aria2(url, outdir / "data" / rel)
            return rel, url, listed, "", "ok", "", secs
        except Exception as exc:  # noqa: BLE001
            return rel, url, listed, "", "failed", str(exc)[:200], 0

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for rel, url, listed, _, state, err, secs in pool.map(get_one, plan):
            dest = outdir / "data" / rel
            if state == "failed":
                print(f"[失败] {rel}: {err}", flush=True)
                failures.append(rel)
                manifest.append(f"{rel.split('/')[0]}\t{rel}\t{listed}\t0\t\t\tfailed")
                continue
            actual = dest.stat().st_size
            ok, note = integrity_ok(dest)
            digest = md5(dest)
            # filelist 里的字节数是权威值；整包 tar 未登记体积时退化为非空判定
            size_ok = actual >= listed * 0.98 if listed else actual > 0
            good = ok and size_ok
            if not good:
                failures.append(rel)
            print(
                f"[{'完成' if good else '异常'}] {rel} {actual} bytes"
                f"（登记 {listed}）{note} md5={digest} 用时 {secs}s",
                flush=True,
            )
            manifest.append(
                f"{rel.split('/')[0]}\t{rel}\t{listed}\t{actual}\t{digest}\t{note}\t{'ok' if good else 'failed'}"
            )

    out = outdir / "meta" / "izar_manifest.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(f"\n清单已写入 {out}")
    if failures:
        print(f"未通过：{failures}")
        return 1
    print("全部下载并校验通过")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())