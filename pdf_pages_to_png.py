#!/usr/bin/env python3
"""
pdf_pages_to_png.py — Render user-specified PDF pages to high-resolution PNGs.

Purpose
-------
Export chosen pages of a PDF (e.g. pages with complex math formulas) as
lossless PNG images at a user-controlled DPI, so that vector math at any
scale — tiny sub/superscripts to large display equations — is sampled with
enough pixels to stay crisp.

Why DPI is the lever
--------------------
A PDF page is 72 pt/inch. Rendering at R DPI samples a glyph of height h pt
into  h * R / 72  pixels. Legibility of the smallest element needs roughly
>= 20-25 px, i.e. R >= ~1440/h. So:
    300 DPI  -> 6 pt sub/superscript = 25 px   (practical floor for math)
    600 DPI  -> 4 pt prime/index      = 33 px   (dense/tiny elements)
    1200 DPI -> for archival / deep zoom
Big formulas simply get more pixels than they need; small ones get the
pixels they need. The whole page is rasterized at the same DPI, so the
"tiny to big" scale range is preserved automatically.

Caveat
------
If the source PDF's math is a SCANNED/EMBEDDED RASTER (not vector), no DPI
setting recovers detail that isn't in the source. Vector (LaTeX/typeset)
PDFs are the good case and benefit directly from higher DPI.

Usage
-----
    python pdf_pages_to_png.py INPUT.pdf --pages 3 --dpi 300
    python pdf_pages_to_png.py INPUT.pdf --pages 1,4,7-9 --dpi 600 -o out/
    python pdf_pages_to_png.py INPUT.pdf --pages all --dpi 300

    --pages   comma list and/or ranges, 1-based: "3", "1-5", "1,3,7-9", "all"
    --dpi     rasterization resolution (default 300; 600 for dense math)
    -o, --out output directory (default: <input_stem>_pages)
    --png     (default) lossless PNG
    --jpeg    lossy JPEG (smaller files, quality via --jpeg-quality)
    --jpeg-quality  1-100 (default 95)
    --prefix  output filename prefix (default: <input_stem>)

Output files:  <prefix>_p003.png, <prefix>_p007.png, ...  (zero-padded, 3 digits)
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

try:
    import pymupdf  # PyMuPDF >= 1.24 (new module name); falls back below
except ImportError:  # pragma: no cover - older installs
    try:
        import fitz as pymupdf
    except ImportError:
        sys.exit(
            "PyMuPDF is required. Install with:\n"
            "    python -m pip install pymupdf"
        )


# --------------------------------------------------------------------------- #
# Page-spec parsing
# --------------------------------------------------------------------------- #
def parse_pages(spec: str, total: int) -> list[int]:
    """Parse '1,3,7-9' / 'all' into a sorted, deduped list of 1-based pages.

    Raises ValueError on bad tokens or out-of-range pages.
    """
    spec = spec.strip()
    if spec.lower() == "all":
        return list(range(1, total + 1))

    pages: set[int] = set()
    for token in spec.split(","):
        token = token.strip()
        if not token:
            continue
        if "-" in token:
            a, _, b = token.partition("-")
            a, b = a.strip(), b.strip()
            if not (a.isdigit() and b.isdigit()):
                raise ValueError(f"bad range: {token!r}")
            lo, hi = int(a), int(b)
            if lo > hi:
                lo, hi = hi, lo
            pages.update(range(lo, hi + 1))
        elif token.isdigit():
            pages.add(int(token))
        else:
            raise ValueError(f"bad page token: {token!r}")

    bad = [p for p in pages if p < 1 or p > total]
    if bad:
        raise ValueError(
            f"page(s) {sorted(bad)} out of range (document has {total} pages)"
        )
    return sorted(pages)


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #
@dataclass
class Result:
    path: Path
    page: int
    width_px: int
    height_px: int


def render_pages(
    pdf_path: Path,
    pages: list[int],
    dpi: int,
    out_dir: Path,
    prefix: str,
    fmt: str = "png",
    jpeg_quality: int = 95,
    progress_cb=None,
) -> list[Result]:
    """Rasterize the given pages to images. Returns per-page results.

    progress_cb: optional callable (done, total, page_number) invoked after
    each page is written — used by the GUI to drive a live progress bar.
    """
    doc = pymupdf.open(pdf_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    zoom = dpi / 72.0
    mat = pymupdf.Matrix(zoom, zoom)
    results: list[Result] = []
    total = len(pages)

    try:
        for i, p in enumerate(pages, start=1):
            page = doc[p - 1]  # 0-based index
            pix = page.get_pixmap(matrix=mat, alpha=False)
            pad = 3 if len(pages) > 99 else 2
            name = f"{prefix}_p{p:0{pad}d}.{fmt}"
            out = out_dir / name
            if fmt == "png":
                pix.save(out)
            else:  # jpeg
                pix.save(out, jpg_quality=jpeg_quality)
            results.append(Result(out, p, pix.width, pix.height))
            if progress_cb is not None:
                progress_cb(i, total, p)
    finally:
        doc.close()
    return results


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Render user-specified PDF pages to high-resolution PNGs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Usage", 1)[-1],
    )
    ap.add_argument("pdf", type=Path, help="input PDF file")
    ap.add_argument("--pages", required=True,
                    help="page list: '3', '1-5', '1,3,7-9', or 'all' (1-based)")
    ap.add_argument("--dpi", type=int, default=300,
                    help="rasterization DPI (default 300; use 600 for dense math)")
    ap.add_argument("-o", "--out", type=Path, default=None,
                    help="output directory (default: <input_stem>_pages)")
    ap.add_argument("--prefix", default=None,
                    help="output filename prefix (default: <input_stem>)")
    fmt = ap.add_mutually_exclusive_group()
    fmt.add_argument("--png", dest="fmt", action="store_const", const="png",
                     default="png", help="lossless PNG (default)")
    fmt.add_argument("--jpeg", dest="fmt", action="store_const", const="jpg",
                     help="lossy JPEG")
    ap.add_argument("--jpeg-quality", type=int, default=95,
                    help="JPEG quality 1-100 (default 95)")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    pdf_path = args.pdf.expanduser().resolve()
    if not pdf_path.is_file():
        print(f"error: no such file: {pdf_path}", file=sys.stderr)
        return 2

    # Open once to get page count for validation.
    probe = pymupdf.open(pdf_path)
    total = probe.page_count
    probe.close()

    try:
        pages = parse_pages(args.pages, total)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if not pages:
        print("error: page list resolved to no pages", file=sys.stderr)
        return 2

    out_dir = args.out or pdf_path.with_name(pdf_path.stem + "_pages")
    prefix = args.prefix or pdf_path.stem

    print(f"PDF      : {pdf_path}  ({total} pages)")
    print(f"Pages    : {pages}")
    print(f"DPI      : {args.dpi}   (zoom = {args.dpi / 72.0:.3f}x)")
    print(f"Output   : {out_dir}  [{args.fmt}]")
    print("-" * 60)

    results = render_pages(
        pdf_path, pages, args.dpi, out_dir, prefix,
        fmt=args.fmt, jpeg_quality=args.jpeg_quality,
    )

    for r in results:
        print(f"  p{r.page:<4} {r.width_px}x{r.height_px}px  ->  {r.path.name}")
    print("-" * 60)
    print(f"Done: {len(results)} image(s) in {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
