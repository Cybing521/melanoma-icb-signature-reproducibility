#!/usr/bin/env python3
"""用 aria2c 多连接下载 GEO 补充文件，并校验体积、写入清单。

单连接从 NCBI FTP 拉取实测约 13 KB/s，多连接（-x16）可到约 96 KB/s，
因此本脚本统一下载到暂存目录，校验通过后再移入 data/。

用法：
  python3 scripts/03_download_aria2.py --outdir /path/to/project            # dry-run
  python3 scripts/03_download_aria2.py --outdir /path/to/project --execute
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import subprocess
import time
import zipfile
from pathlib import Path

# 复用 02 的清单逻辑，避免两处各写一份目标列表
import importlib.util

_util = importlib.util.spec_from_file_location(
    "dl02", Path(__file__).resolve().parent / "02_download_suppl.py"
)
assert _util and _util.loader
dl02 = importlib.util.module_from_spec(_util)
_util.loader.exec_module(dl02)


def fetch_with_aria2(url: str, dest: Path, conns: int = 16, retries: int = 5) -> int:
    """返回下载耗时（秒）。

    跳过条件：文件已存在、非零、**且没有 .aria2 控制文件**。
    带 .aria2 旁车文件说明上次下载未完成，此时必须调用 aria2c 让它断点续传，
    否则残片会被反复判为 integrity_failed 却永不重下。
    """
    sidecar = dest.with_name(dest.name + ".aria2")
    if dest.exists() and dest.stat().st_size > 0 and not sidecar.exists():
        print(f"[跳过] {dest.name} 已存在 {dest.stat().st_size} bytes", flush=True)
        return 0
    if sidecar.exists():
        print(f"[续传] {dest.name} 上次未完成（{sidecar.name} 存在），由 aria2c 续传", flush=True)
    cmd = [
        "aria2c", "-x", str(conns), "-s", str(conns), "-k", "1M",
        "--file-allocation=none", "--summary-interval=0",
        "--max-tries", str(retries), "--retry-wait=3",
        "-d", str(dest.parent), "-o", dest.name, url,
    ]
    start = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"aria2c 失败 {url}: {proc.stderr[-500:]}")
    # aria2c 正常完成会自行清掉控制文件；残留则视为未完成
    if sidecar.exists():
        raise RuntimeError(f"aria2c 退出但控制文件仍存在，判定为未完成: {dest.name}")
    return int(time.time() - start)


def md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def size_ok(actual: int, listed: int, tol: float = 0.02) -> bool:
    """GEO 目录列表的体积四舍五入到 MB/KM/G，不是精确字节数。

    实测：列表显示 16M 的文件实际 16,354,905 字节。故按容差判定，
    容差需覆盖列表舍入误差（最大约半个单位），并留出截断余量。
    """
    if listed <= 0:
        return actual > 0
    slack = max(listed * tol, 1 << 20)  # 至少容忍 1 MB 的舍入/截断差
    return actual >= listed - slack


def integrity_ok(path: Path) -> tuple[bool, str]:
    """按文件类型做完整性检查，比比对体积更能证明下载未截断。"""
    name = path.name.lower()
    try:
        if name.endswith(".gz"):
            with gzip.open(path, "rb") as fh:
                while fh.read(1 << 22):
                    pass
            return True, "gzip 解压通过"
        if name.endswith((".tar", ".tar.gz")):
            proc = subprocess.run(
                ["tar", "tf", str(path)], capture_output=True, text=True, timeout=600
            )
            n = len([x for x in proc.stdout.splitlines() if x.strip()])
            return proc.returncode == 0 and n > 0, f"tar 可列出 {n} 项"
        if name.endswith((".xlsx", ".zip")):
            with zipfile.ZipFile(path) as zf:
                bad = zf.testzip()
                return bad is None, f"zip 结构完整（{len(zf.namelist())} 项）"
        with path.open("rb") as fh:
            head = fh.read(4096)
        return bool(head), "可读非空"
    except Exception as exc:  # noqa: BLE001
        return False, f"完整性检查失败：{exc}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--conns", type=int, default=16)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    staging = outdir / "data_staging"
    staging.mkdir(parents=True, exist_ok=True)

    plan: list[tuple[str, str, int]] = []
    for gse, role in dl02.TARGETS.items():
        for name, size in dl02.list_files(gse):
            plan.append((gse, name, size))
    grand = sum(size for _, _, size in plan)
    print(f"共 {len(plan)} 个文件，合计 {grand / 1024**2:.1f} MB")
    for gse, name, size in plan:
        print(f"  {gse:<10} {size / 1024**2:8.2f} MB  {name}")
    if not args.execute:
        print("dry-run 结束，加 --execute 执行")
        return 0

    manifest = ["gse\tfile\tlisted_bytes\tactual_bytes\tmd5\tintegrity\tstatus\tlocation"]
    failures = []
    for gse, name, size in plan:
        # 已落库的文件不再重下，但仍重新核验并写进清单，
        # 保证清单描述的是磁盘真实状态，而不是上一轮的下载意图。
        final = outdir / "data" / gse / name
        if final.exists() and final.stat().st_size > 0:
            actual = final.stat().st_size
            ok_int, note = integrity_ok(final)
            ok = size_ok(actual, size) and ok_int
            digest = md5(final)
            if ok:
                status = "ok"
                print(f"[已落库] {gse}/{name} {actual} bytes（列表 {size}）{note}", flush=True)
            else:
                status = "size_mismatch" if not size_ok(actual, size) else "integrity_failed"
                print(f"[已落库但不合格] {gse}/{name} {actual} bytes（列表 {size}）{note}", flush=True)
                failures.append(f"{gse}/{name}")
            manifest.append(
                f"{gse}\t{name}\t{size}\t{actual}\t{digest}\t{note}\t{status}\tdata/{gse}"
            )
            continue

        sdir = staging / gse
        sdir.mkdir(parents=True, exist_ok=True)
        url = f"{dl02.supp_url(gse)}{name}"
        try:
            secs = fetch_with_aria2(url, sdir / name, args.conns)
        except Exception as exc:  # noqa: BLE001
            print(f"[失败] {gse}/{name}: {exc}", flush=True)
            failures.append(f"{gse}/{name}")
            manifest.append(f"{gse}\t{name}\t{size}\t0\t\t\tfailed\t-")
            continue

        src = sdir / name
        actual = src.stat().st_size
        ok_size = size_ok(actual, size)
        ok_int, note = integrity_ok(src)
        # 体积在容差内且完整性检查通过才算通过
        ok = ok_size and ok_int
        digest = md5(src)
        if ok:
            dest_dir = outdir / "data" / gse
            dest_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src), str(dest_dir / name))
        else:
            status = "size_mismatch" if not ok_size else "integrity_failed"
            failures.append(f"{gse}/{name}")
        print(
            f"[{'完成' if ok else '异常'}] {gse}/{name} {actual} bytes"
            f"（列表 {size}）{note} md5={digest} 用时 {secs}s",
            flush=True,
        )
        manifest.append(
            f"{gse}\t{name}\t{size}\t{actual}\t{digest}\t{note}\t{'ok' if ok else status}\t"
            f"{'data/' + gse if ok else 'data_staging/' + gse}"
        )

    out = outdir / "meta" / "suppl_manifest.tsv"
    out.write_text("\n".join(manifest) + "\n", encoding="utf-8")
    print(f"\n清单已写入 {out}")
    if failures:
        print(f"以下文件未通过校验：{failures}")
        return 1
    print("全部文件下载并校验通过")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())