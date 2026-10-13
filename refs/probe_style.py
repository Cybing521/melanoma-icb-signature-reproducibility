#!/usr/bin/env python3
"""从参考文献 PDF 中提取图件风格证据：嵌入字体族 + 主色板。

字体：PyMuPDF 只能给出子集名（如 `CPPJKO+AdvOT511ca30e.B`），因此把嵌入字体
文件取出，用 fontTools 读 name 表，取 family（id 1/16）真实族名。
配色：把每页栅格化后统计饱和度 >45 的彩色像素，按 16 级量化取众数簇。
只做测量，不做判断；结论写入 `docs/13_图件风格仿制依据.md`。
"""
from __future__ import annotations

import collections
import io
from pathlib import Path

import fitz
import numpy as np
from PIL import Image
from fontTools.ttLib import TTFont

HERE = Path(__file__).resolve().parent
PAPERS = HERE / "papers"


def family_of(data: bytes, ext: str) -> str:
    """从嵌入字体字节读真实族名；失败则退回子集名。"""
    try:
        if ext in ("cff", "cff2"):           # OTF/CFF 用 TTFont 读不了，走 sfnt 兜底
            raise ValueError("cff")
        f = TTFont(io.BytesIO(data), fontNumber=0, lazy=True)
        for nid in (16, 1):
            for rec in f["name"].names:
                if rec.nameID == nid:
                    try:
                        return rec.toUnicode().strip()
                    except Exception:  # noqa: BLE001
                        pass
        return "?"
    except Exception:  # noqa: BLE001
        return "(CFF subset)"


def pdf_fonts(pdf: Path) -> list[str]:
    doc = fitz.open(pdf)
    fams: collections.Counter = collections.Counter()
    sizes: collections.Counter = collections.Counter()
    for page in doc:
        seen_xref = set()
        for f in page.get_fonts():
            xref, ext = f[0], f[1]
            if xref in seen_xref:
                continue
            seen_xref.add(xref)
            try:
                name, _n, _t, buf = doc.extract_font(xref)
                fams[family_of(buf, ext)] += 1
            except Exception:  # noqa: BLE001
                fams[f[3].split("+")[-1]] += 1
        for block in page.get_text("dict")["blocks"]:
            for line in block.get("lines", []):
                for span in line.get("spans", []):
                    n = len(span["text"].strip())
                    if n:
                        sizes[round(span["size"], 1)] += n
    n_pages = doc.page_count
    doc.close()
    print(f"\n=== {pdf.name}  ({n_pages} 页) ===")
    print("  字体族:", ", ".join(f"{k}" for k, _ in fams.most_common(5)))
    print("  字号pt(字符数):", ", ".join(f"{s}pt×{c}" for s, c in
                                          sorted(sizes.items(), key=lambda kv: -kv[1])[:6]))
    return [k for k, _ in fams.most_common(5)]


def pdf_palette(pdf: Path, max_pages: int = 16) -> list[tuple[str, float]]:
    doc = fitz.open(pdf)
    acc: collections.Counter = collections.Counter()
    for i in range(min(max_pages, doc.page_count)):
        pix = doc[i].get_pixmap(dpi=100)
        img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        a = np.asarray(img).reshape(-1, 3).astype(np.int16)
        mx, mn = a.max(1), a.min(1)
        m = (mx - mn) > 45
        if m.any():
            for px in a[m]:
                acc[(int(px[0]) // 16 * 16, int(px[1]) // 16 * 16,
                     int(px[2]) // 16 * 16)] += 1
    doc.close()
    if not acc:
        print("  色板: 无显著彩色像素")
        return []
    top = acc.most_common(12)
    denom = sum(c for _, c in top)
    print("  色板:")
    out = []
    for (r, g, b), c in top:
        pct = 100 * c / denom
        print(f"    #{r:02X}{g:02X}{b:02X}  {pct:5.1f}%")
        out.append((f"#{r:02X}{g:02X}{b:02X}", pct))
    return out


def main() -> None:
    for pdf in sorted(PAPERS.glob("*.pdf")):
        try:
            pdf_fonts(pdf)
            pdf_palette(pdf)
        except Exception as exc:  # noqa: BLE001
            print(f"[FAIL] {pdf.name}: {exc}")


if __name__ == "__main__":
    main()
