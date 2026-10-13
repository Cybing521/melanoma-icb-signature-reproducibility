#!/usr/bin/env python3
"""把新参考文献的排版统一到稿件既有 13 条的体例。

PubMed efetch 给的字段与稿件既有体例有三处不一致，统一后才不会在稿件里
出现两种写法：

| 项 | PubMed efetch | 稿件既有体例 | 本脚本处理 |
|---|---|---|---|
| 作者 | `Jiang P; Gu S` | `Auslander N.; Zhang G.` | 名字首字母后补句点 |
| 刊名 | `Nat Med` | `Nat. Med.` | 长度 ≥2 的词补句点（`U S A` 等单字母保留） |
| 页码 | `640-8`（省略式） | `1545–1549`（全写 + 破折号） | 展开省略式，破折号统一为 en dash |

用法：
    python refs/normalize_refs.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CONTENT = ROOT / "manuscript" / "content.json"
RECORDS = ROOT / "refs" / "pmid_records.json"


# 刊名缩写里本身已是缩写、不能加句点的词
JOURNAL_KEEP = {"BMC", "BMJ", "JAMA", "PLoS", "PLOS", "eLife"}


def fix_authors(s: str) -> str:
    """`Jiang P; Gu S; Traugh N.` → `Jiang P.; Gu S.; Traugh N.`

    每个作者各自补句点后直接拼接即可——末尾不能再补一次，否则末位作者
    会变成 `Simon R..`。
    """
    s = s.rstrip(".")          # expand_refs 已在整串末尾补过一个句点，先去掉
    parts = []
    for a in s.split("; "):
        a = a.rstrip(".")
        if " " in a and len(a.split()[-1]) == 1:      # 末段是单字母首字母
            a = f"{a}."
        parts.append(a)
    return "; ".join(parts)


def fix_journal(s: str) -> str:
    """`Proc Natl Acad Sci U S A` → `Proc. Natl. Acad. Sci. U S A`

    三个例外：单词刊名（`Epidemiology`、`Radiology`）本身就是全称，不加句点；
    已是缩写的词（`BMC`、`BMJ`）加了反而变成 `BMC.`；刊名**末词**后不加句点，
    因为它后面紧跟年份，加了会写成 `Nat. Methods. 2015`。
    """
    toks = s.split()
    if len(toks) == 1:
        return s
    out = []
    for i, tok in enumerate(toks):
        last = i == len(toks) - 1
        if last or tok in JOURNAL_KEEP or len(tok) == 1 or tok.endswith("."):
            out.append(tok)
        else:
            out.append(f"{tok}.")
    return " ".join(out)


def expand_pages(pages: str) -> str:
    """`640-8` → `640–648`；`895-900` → `895–900`（en dash）。"""
    m = re.fullmatch(r"(\d+)-(\d+)", pages.strip())
    if not m:
        return pages.replace("-", "–")
    lo, hi = m.group(1), m.group(2)
    if len(hi) < len(lo):
        # 省略式：按末位补齐被省略的高位数字
        hi = lo[:len(lo) - len(hi)] + hi
    return f"{lo}–{hi}"


def main() -> int:
    doc = json.loads(CONTENT.read_text())
    blocks = doc["blocks"]
    head = next(i for i, b in enumerate(blocks)
                if b["type"] == "heading" and b["text"] == "References")

    # 编号重排后新旧条目已交错（旧的落在 1,2,8-14,16,18,20,23），不能按位置切分，
    # 只能按题名认出这 12 条新条目。
    titles = [r["title"] for r in json.loads(RECORDS.read_text()).values()]
    refs = blocks[head + 1:]
    changed = 0
    for n, b in enumerate(refs):
        t = b["text"]
        if not any(ti in t for ti in titles):
            continue
        # 结构固定为 "<序号>. <作者> <题名> <刊名> <年>, <卷>, <页>. <注解>"
        body = re.sub(r"^\d+\.\s*", "", t)
        # 作者段 = 题名之前、以 ". " 收尾的一段（作者内部以 "; " 分隔，无 ". "）
        m = re.match(r"^(.*?\.)\s+(.*)$", body, re.S)
        authors, rest = m.group(1), m.group(2)
        new = fix_authors(authors)
        # 刊名 / 年 / 卷 / 页 是尾部形如 "Jrnl 2017, 127, 2930-2940." 的一段
        m2 = re.match(r"^(.*?)\s(\d{4}),\s([^,]+),\s([0-9g]+[-–]?[0-9]*)\.\s(.*)$",
                      rest, re.S)
        if not m2:
            raise SystemExit(f"参考文献尾部解析失败: {t[:120]}")
        head_txt, year, vol, pages, note = m2.groups()
        # head_txt = "<题名> <刊名>"，题名以句号收尾，刊名是其后的一串词
        title, journal = head_txt.rsplit(". ", 1)
        title += "."
        new += (f" {title} {fix_journal(journal)} {year}, {vol}, "
                f"{expand_pages(pages)}. {note}")
        b["text"] = f"{n + 1}. {new}"
        changed += 1

    CONTENT.write_text(json.dumps(doc, ensure_ascii=False, indent=1))
    print(f"已统一 {changed} 条新参考文献的排版体例")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
