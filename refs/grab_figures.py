#!/usr/bin/env python3
"""抓取参考文献的图件原图，用于提取其视觉风格（字体 / 配色 / 版式）。

路线：Europe PMC 解析 PMCID → PMC 文章页 → 抽出 cdn.ncbi.nlm.nih.gov 上的
blobs 图片直链 → 下载。订阅制且不在 PMC OA 子集的文献记为 SKIP。
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
IMG = HERE / "figures_style"
IMG.mkdir(exist_ok=True)

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36")

# 自然语言期刊（可直接下 PDF）与 Cell 系（只能抓 PMC 原图）分开处理
PDF_TARGETS = {
    "Vathiotis2022_npjPrecisOncol": "10.1038/s41698-022-00330-3",
    "Lauss2024_NatCommun":         "10.1038/s41467-024-47425-y",
    "DiPietro2026_NatCommun":      "10.1038/s41467-026-74076-y",
}

PMC_TARGETS = {
    "Hugo2016_Cell":      "PMC4808437",
    "Charoentong2017_CellRep": None,   # PMCID 由 PMID 解析
    "Riaz2017_Cell":      None,
    "Auslander2018_NatMed": None,
}

PMC_BY_PMID = {
    "Charoentong2017_CellRep": "28052254",
    "Riaz2017_Cell":          "29033130",
    "Auslander2018_NatMed":   "30127394",
}


def curl(url: str, out: Path | None = None) -> bytes:
    cmd = ["curl", "-sL", "--noproxy", "*", "-A", UA, "--max-time", "90", url]
    if out:
        cmd += ["-o", str(out)]
        subprocess.run(cmd, check=False)
        return b""
    return subprocess.run(cmd, capture_output=True, check=False).stdout


def pmcid_of(pmid: str) -> str | None:
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
           f"?query=EXT_ID:{pmid}&format=json&resultType=core")
    try:
        rec = json.loads(curl(url))["resultList"]["result"][0]
        return rec.get("pmcid")
    except Exception as exc:  # noqa: BLE001
        print(f"    PMCID lookup failed PMID {pmid}: {exc}")
        return None


def grab_figures(name: str, pmcid: str, max_n: int = 8) -> int:
    html = curl(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/")
    text = html.decode("utf-8", "ignore")
    urls = re.findall(r'src="(https://cdn\.ncbi\.nlm\.nih\.gov/pmc/blobs/[^"]+\.(?:jpg|png|gif))"', text)
    seen, keep = set(), []
    for u in urls:
        if u not in seen:
            seen.add(u)
            keep.append(u)
    keep = keep[:max_n]
    if not keep:
        print(f"[SKIP] {name}: PMC 页内无 blobs 直链")
        return 0
    d = IMG / name
    d.mkdir(exist_ok=True)
    n = 0
    for i, u in enumerate(keep, 1):
        dest = d / f"fig{i}.jpg"
        curl(u, dest)
        if dest.exists() and dest.stat().st_size > 8000 and dest.read_bytes()[:2] == b"\xff\xd8":
            n += 1
        else:
            dest.unlink(missing_ok=True)
    print(f"[ OK ] {name}: {n} 张原图 -> {d}")
    return n


def main() -> int:
    ok = 0
    for name, doi in PDF_TARGETS.items():
        dest = HERE / "papers" / f"{name}.pdf"
        url = f"https://www.nature.com/articles/{doi.split('/')[-1]}.pdf"
        curl(url, dest)
        if dest.exists() and dest.read_bytes()[:4] == b"%PDF":
            print(f"[ OK ] {name}.pdf  {dest.stat().st_size/1024:.0f} KB")
            ok += 1
        else:
            print(f"[FAIL] {name}: {url}")
            dest.unlink(missing_ok=True)

    for name, pmcid in PMC_TARGETS.items():
        mc = pmcid or pmcid_of(PMC_BY_PMID[name])
        if not mc:
            print(f"[SKIP] {name}: 未在 PMC 检索到（可能非 OA 子集）")
            continue
        ok += grab_figures(name, mc)

    print(f"\n成功 {ok} 项")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
