#!/usr/bin/env python3
"""回 PubMed 核对参考文献，生成 MDPI 体例的参考文献条目。

为什么自己再核一遍
------------------
本项目 `docs/12` 记过一次教训：稿件参考文献表曾凭记忆写成，10 条里 5 条的
题名、刊卷页或作者与数据库记录不符。因此任何进入稿件的条目都必须回到
PubMed efetch 的 XML 取字段，不采信任何转述。

本脚本对给定 PMID 列表逐条 efetch，输出：
  - 作者（全部，与稿件既有体例一致）
  - 题名 / 刊名缩写 / 年 / 卷 / 页
  - 作者总数（供人工核对「是否该用 et al.」）
  - PMID 与 DOI（DOI 以 Europe PMC 交叉核对，PubMed 的 DOI 字段有过错配）

用法：
    python refs/verify_refs.py 30127393 28650338 ...
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

HERE = Path(__file__).resolve().parent
CACHE = HERE / "pmid_records.json"


def curl(url: str) -> str:
    return subprocess.run(
        ["curl", "-s", "--noproxy", "*", "--max-time", "60", url],
        capture_output=True, check=False).stdout.decode("utf-8", "ignore")


def efetch(pmids: list[str]) -> dict[str, dict]:
    """批量 efetch，返回 {pmid: 记录}。"""
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=pubmed&id={','.join(pmids)}&retmode=xml")
    root = ET.fromstring(curl(url))
    out = {}
    for art in root.findall(".//PubmedArticle"):
        pmid = art.findtext(".//PMID")
        j = art.find(".//Journal")
        issue = j.findtext("Issue") or ""
        page = art.findtext(".//MedlinePgn") or ""
        # 文章号（如 BMJ 的 g7594、BMC 的 91）只在没有 MedlinePgn 时才用 ELocationID。
        # 注意 JCI 记录里 pii 形如 "91190"（取自 JCI91190），不是页码，误用会得到错页。
        if not page:
            for eid in art.findall(".//ELocationID"):
                if eid.get("EIdType") == "pii" and eid.text:
                    page = eid.text
                    break
        authors = []
        for a in art.findall(".//AuthorList/Author"):
            ln, ini = a.findtext("LastName"), a.findtext("Initials")
            if ln:
                authors.append(f"{ln} {ini}".strip())
            elif a.findtext("CollectiveName"):
                authors.append(a.findtext("CollectiveName"))
        year = j.findtext(".//Year") or (j.findtext(".//MedlineDate") or "")[:4]
        out[pmid] = {
            "authors": authors,
            "n_authors": len(authors),
            "title": (art.findtext(".//ArticleTitle") or "").strip(),
            "journal": j.findtext("ISOAbbreviation") or j.findtext("Title") or "",
            "year": year,
            "volume": j.findtext("JournalIssue/Volume") or "",
            "issue": issue,
            "pages": page,
        }
    return out


def epmc_doi(pmid: str) -> str:
    url = ("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
           f"?query=EXT_ID:{pmid}&format=json&resultType=core")
    try:
        rec = json.loads(curl(url))["resultList"]["result"][0]
        return rec.get("doi") or ""
    except Exception:  # noqa: BLE001
        return ""


def main(pmids: list[str]) -> int:
    recs = efetch(pmids)
    CACHE.write_text(json.dumps(recs, ensure_ascii=False, indent=1))
    print(f"{'PMID':<10}{'n':>3}  {'刊':<26}{'年':<6}{'卷':<6}{'页':<12}DOI")
    for p in pmids:
        r = recs.get(p)
        if not r:
            print(f"{p:<10}  !! PubMed 无记录")
            continue
        doi = epmc_doi(p)
        print(f"{p:<10}{r['n_authors']:>3}  {r['journal']:<26}{r['year']:<6}"
              f"{r['volume']:<6}{r['pages']:<12}{doi}")
        print(f"           {r['title']}")
        time.sleep(0.35)
    missing = [p for p in pmids if p not in recs]
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
