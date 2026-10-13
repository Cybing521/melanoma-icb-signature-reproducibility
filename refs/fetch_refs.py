#!/usr/bin/env python3
"""下载参考文献原文 PDF，用于提取其图件视觉风格。

只走 Europe PMC 的 OA 路由；订阅制文献记录为 SKIP，不做任何绕过。
NCBI / EBI 域名在本网段必须 `curl --noproxy '*'`，Python urllib 走代理会 SSL 失败。
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPERS = HERE / "papers"
PAPERS.mkdir(exist_ok=True)

TARGETS = {
    "Auslander2018_IMPRES":      ("30127394", "Nat Med - IMPRES"),
    "Charoentong2017_IPS":       ("28052254", "Cell Rep - Immunophenoscore"),
    "Riaz2017_GSE91061":         ("29033130", "Cell - GSE91061 source"),
    "Hugo2016_GSE78220":         ("26997480", "Cell - GSE78220 source"),
    "Vathiotis2022_GSE215868":   ("36522538", "NPJ Precis Oncol - GSE215868 source"),
    "Lauss2024_GSE244982":       (None,       "Nat Commun - GSE244982 source"),
    "DiPietro_GSE294272":        (None,       "Nat Commun - GSE294272 source"),
}


def curl(url: str, out: Path | None = None) -> bytes:
    cmd = ["curl", "-sL", "--noproxy", "*", "--max-time", "90", url]
    if out:
        cmd += ["-o", str(out)]
        subprocess.run(cmd, check=False)
        return b""
    return subprocess.run(cmd, capture_output=True, check=False).stdout


def lookup(pmid: str) -> dict:
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
           f"?query=EXT_ID:{pmid}&format=json&resultType=core")
    try:
        data = json.loads(curl(url))
        return data["resultList"]["result"][0]
    except Exception as exc:  # noqa: BLE001
        print(f"    lookup failed: {exc}")
        return {}


def best_pdf(rec: dict) -> str | None:
    urls = (rec.get("fullTextUrlList") or {}).get("fullTextUrl", [])
    cands = [u["url"] for u in urls
             if u.get("availabilityCode") == "OA" and u.get("documentStyle") == "pdf"]
    # Europe PMC 的 pdf=render 最稳，优先
    cands.sort(key=lambda u: (0 if "europepmc.org" in u and "pdf=render" in u else 1))
    return cands[0] if cands else None


def resolve_missing() -> dict[str, str]:
    """两条 Nat Commun 的 PMID 未硬编码，按刊+卷+页硬定位（沿用 docs/12 的做法）。"""
    found: dict[str, str] = {}
    for key, query in (
        ("Lauss2024_GSE244982", 'JOURNAL:"Nat Commun" AND VOLUME:15 AND PAGE:3075'),
        ("DiPietro_GSE294272", 'TITLE:"Tumor-resident T cells and dendritic cells"'),
    ):
        url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
               f"?query={query.replace(' ', '%20')}&format=json&resultType=core")
        try:
            res = json.loads(curl(url)).get("resultList", {}).get("result", [])
            if res:
                found[key] = res[0]["pmid"]
                print(f"    resolved {key} -> PMID {res[0]['pmid']}")
        except Exception as exc:  # noqa: BLE001
            print(f"    resolve failed {key}: {exc}")
    return found


def main() -> int:
    resolved = resolve_missing()
    targets = {k: (resolved.get(k) or p, note) for k, (p, note) in TARGETS.items()}

    ok, skipped = [], []
    for name, (pmid, note) in targets.items():
        if not pmid:
            print(f"[SKIP] {name}: PMID 未能定位 ({note})")
            skipped.append(name)
            continue
        rec = lookup(pmid)
        if not rec:
            print(f"[SKIP] {name}: Europe PMC 无记录 (PMID {pmid})")
            skipped.append(name)
            continue
        pdf = best_pdf(rec)
        title = (rec.get("title") or "")[:70]
        if not pdf:
            print(f"[SKIP] {name}: 无 OA PDF (PMID {pmid}) — {title}")
            skipped.append(name)
            continue
        dest = PAPERS / f"{name}.pdf"
        curl(pdf, dest)
        head = dest.read_bytes()[:5] if dest.exists() else b""
        if head == b"%PDF-":
            kb = dest.stat().st_size / 1024
            print(f"[ OK ] {name}.pdf  {kb:6.0f} KB  PMID {pmid}  {title}")
            ok.append(name)
        else:
            print(f"[FAIL] {name}: 下载内容不是 PDF ({pdf})")
            skipped.append(name)

    print(f"\n成功 {len(ok)} / {len(TARGETS)}；跳过 {len(skipped)}: {skipped}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
