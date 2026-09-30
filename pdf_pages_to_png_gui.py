#!/usr/bin/env python3
r"""
pdf_pages_to_png_gui.py — Tiny tkinter GUI for pdf_pages_to_png.py.

Renders user-specified PDF pages to high-resolution PNG/JPEG at a chosen DPI,
reusing the exact same render core as the CLI (single source of truth).

Run:
    python pdf_pages_to_png_gui.py            # from the folder, or
    .venv\Scripts\python pdf_pages_to_png_gui.py

Why DPI: a glyph of height h pt rasterizes into h*DPI/72 px. For the
smallest math element to stay legible you need ~20-25 px, i.e. DPI >= ~1440/h.
    300 DPI -> practical floor for typical math (6 pt sub/superscript = 25 px)
    600 DPI -> dense/tiny elements (4 pt prime/index = 33 px)
    1200 DPI -> archival / deep zoom
Vector (LaTeX/typeset) PDFs benefit directly; scanned-raster math does not.
"""
from __future__ import annotations

import os
import queue
import subprocess
import sys
import threading
import traceback
from dataclasses import dataclass
from pathlib import Path

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

# Reuse the tested render core from the CLI module (same directory).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pdf_pages_to_png import parse_pages, render_pages  # noqa: E402

try:
    import pymupdf
except ImportError:  # pragma: no cover
    pymupdf = None

DPI_PRESETS = {
    "300  (math floor — typical sub/superscripts)": 300,
    "600  (dense math — tiny primes/indexes)": 600,
    "1200 (archival / deep zoom)": 1200,
    "Custom…": 0,
}


@dataclass
class Job:
    pdf: Path
    pages: list[int]
    dpi: int
    out_dir: Path
    prefix: str
    fmt: str
    jpeg_quality: int


class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("PDF Pages → PNG (high-res)")
        root.geometry("680x520")
        root.minsize(600, 460)

        self.q: queue.Queue = queue.Queue()
        self.job: Job | None = None
        self._thread: threading.Thread | None = None

        self._build()
        self.root.after(80, self._poll)

    # ------------------------------------------------------------------ #
    # UI construction
    # ------------------------------------------------------------------ #
    def _build(self):
        pad = {"padx": 10, "pady": 5}
        outer = ttk.Frame(self.root)
        outer.pack(fill="both", expand=True)

        # --- PDF row ---------------------------------------------------- #
        row = ttk.Frame(outer)
        row.pack(fill="x", **pad)
        ttk.Label(row, text="PDF file:").pack(side="left")
        self.var_pdf = tk.StringVar()
        ent = ttk.Entry(row, textvariable=self.var_pdf)
        ent.pack(side="left", fill="x", expand=True, padx=(8, 6))
        ttk.Button(row, text="Browse…", command=self._pick_pdf).pack(side="left")

        # --- Pages row -------------------------------------------------- #
        row = ttk.Frame(outer)
        row.pack(fill="x", **pad)
        ttk.Label(row, text="Pages:").pack(side="left")
        self.var_pages = tk.StringVar(value="all")
        ttk.Entry(row, textvariable=self.var_pages).pack(
            side="left", fill="x", expand=True, padx=(8, 6)
        )
        self.lbl_total = ttk.Label(row, text="", foreground="#666")
        self.lbl_total.pack(side="left")

        # --- DPI row ---------------------------------------------------- #
        row = ttk.Frame(outer)
        row.pack(fill="x", **pad)
        ttk.Label(row, text="DPI:").pack(side="left")
        self.var_dpi = tk.StringVar(value=list(DPI_PRESETS)[0])
        ttk.Combobox(
            row, textvariable=self.var_dpi, state="readonly",
            values=list(DPI_PRESETS), width=42,
        ).pack(side="left", padx=(8, 6))
        self.var_dpi_custom = tk.StringVar(value="600")
        self.ent_dpi_custom = ttk.Entry(row, textvariable=self.var_dpi_custom, width=8)
        self.ent_dpi_custom.pack(side="left")
        self._sync_dpi_entry()
        self.var_dpi.trace_add("write", self._sync_dpi_entry)

        # --- Format row ------------------------------------------------- #
        row = ttk.Frame(outer)
        row.pack(fill="x", **pad)
        ttk.Label(row, text="Format:").pack(side="left")
        self.var_fmt = tk.StringVar(value="png")
        ttk.Radiobutton(row, text="PNG (lossless)", value="png",
                        variable=self.var_fmt).pack(side="left", padx=(8, 4))
        ttk.Radiobutton(row, text="JPEG", value="jpg",
                        variable=self.var_fmt).pack(side="left")
        ttk.Label(row, text="quality").pack(side="left", padx=(8, 2))
        self.var_jq = tk.StringVar(value="95")
        ttk.Spinbox(row, from_=1, to=100, textvariable=self.var_jq,
                    width=5).pack(side="left")

        # --- Output row ------------------------------------------------- #
        row = ttk.Frame(outer)
        row.pack(fill="x", **pad)
        ttk.Label(row, text="Output:").pack(side="left")
        self.var_out = tk.StringVar()
        ttk.Entry(row, textvariable=self.var_out).pack(
            side="left", fill="x", expand=True, padx=(8, 6)
        )
        ttk.Button(row, text="Browse…", command=self._pick_out).pack(side="left")

        # --- Render button + progress ----------------------------------- #
        row = ttk.Frame(outer)
        row.pack(fill="x", **pad)
        self.btn = ttk.Button(row, text="Render", command=self._start)
        self.btn.pack(side="left")
        self.progress = ttk.Progressbar(row, mode="determinate", length=220)
        self.progress.pack(side="left", fill="x", expand=True, padx=(10, 0))
        self.lbl_prog = ttk.Label(row, text="")
        self.lbl_prog.pack(side="left", padx=(8, 0))

        # --- Log -------------------------------------------------------- #
        logf = ttk.LabelFrame(outer, text="Log")
        logf.pack(fill="both", expand=True, **pad)
        self.txt = tk.Text(logf, height=10, state="disabled", wrap="word",
                           font=("Consolas", 9))
        sb = ttk.Scrollbar(logf, command=self.txt.yview)
        self.txt.configure(yscrollcommand=sb.set)
        self.txt.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")

        # --- Hint ------------------------------------------------------- #
        ttk.Label(
            outer,
            text="Vector math: 300 DPI = floor, 600 = dense.  "
                 "Scanned-raster math cannot be improved by DPI.",
            foreground="#888",
        ).pack(anchor="w", padx=12, pady=(0, 6))

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    def log(self, msg: str):
        self.txt.configure(state="normal")
        self.txt.insert("end", msg + "\n")
        self.txt.see("end")
        self.txt.configure(state="disabled")

    def _sync_dpi_entry(self, *_):
        val = DPI_PRESETS.get(self.var_dpi.get(), 0)
        self.ent_dpi_custom.configure(state="normal" if val == 0 else "disabled")
        if val != 0:
            self.var_dpi_custom.set(str(val))

    def _pick_pdf(self):
        p = filedialog.askopenfilename(
            title="Choose PDF", filetypes=[("PDF", "*.pdf"), ("All", "*.*")]
        )
        if p:
            self.var_pdf.set(p)
            self._refresh_total()

    def _pick_out(self):
        d = filedialog.askdirectory(title="Choose output directory")
        if d:
            self.var_out.set(d)

    def _refresh_total(self):
        pdf = self.var_pdf.get().strip()
        if not pdf or pymupdf is None:
            self.lbl_total.configure(text="")
            return
        try:
            doc = pymupdf.open(pdf)
            n = doc.page_count
            doc.close()
            self.lbl_total.configure(text=f"({n} pages)")
        except Exception:
            self.lbl_total.configure(text="(unreadable)")

    def _resolved_dpi(self) -> int:
        preset = DPI_PRESETS.get(self.var_dpi.get(), 0)
        if preset != 0:
            return preset
        try:
            v = int(self.var_dpi_custom.get())
            if 1 <= v <= 4800:
                return v
        except ValueError:
            pass
        raise ValueError("Enter a valid custom DPI (1–4800)")

    # ------------------------------------------------------------------ #
    # Job lifecycle
    # ------------------------------------------------------------------ #
    def _start(self):
        if self._thread is not None and self._thread.is_alive():
            messagebox.showinfo("Busy", "A render is already in progress.")
            return
        if pymupdf is None:
            messagebox.showerror(
                "Missing dependency",
                "PyMuPDF is not installed.\nRun:  python -m pip install pymupdf",
            )
            return

        pdf = Path(self.var_pdf.get().strip()).expanduser()
        if not pdf.is_file():
            messagebox.showerror("No file", "Choose an existing PDF first.")
            return
        try:
            dpi = self._resolved_dpi()
            pages = parse_pages(self.var_pages.get(),
                                pymupdf.open(pdf).page_count)
            fmt = self.var_fmt.get()
            jq = int(self.var_jq.get())
        except ValueError as e:
            messagebox.showerror("Invalid input", str(e))
            return

        out_dir = Path(self.var_out.get().strip()).expanduser() if \
            self.var_out.get().strip() else pdf.with_name(pdf.stem + "_pages")
        prefix = pdf.stem

        self.job = Job(pdf, pages, dpi, out_dir, prefix, fmt, jq)
        self.btn.configure(state="disabled")
        self.progress.configure(value=0, maximum=len(pages))
        self.log(f"\n=== Render: {pdf.name} | pages {pages} | {dpi} DPI | {fmt} ===")
        self._thread = threading.Thread(target=self._worker, daemon=True)
        self._thread.start()

    def _worker(self):
        j = self.job
        assert j is not None
        try:
            def cb(done, total, page):
                self.q.put(("progress", done, total, page))
            results = render_pages(
                j.pdf, j.pages, j.dpi, j.out_dir, j.prefix,
                fmt=j.fmt, jpeg_quality=j.jpeg_quality, progress_cb=cb,
            )
            self.q.put(("done", results))
        except Exception:
            self.q.put(("error", traceback.format_exc()))

    def _poll(self):
        try:
            while True:
                item = self.q.get_nowait()
                kind = item[0]
                if kind == "progress":
                    _, done, total, page = item
                    self.progress.configure(value=done)
                    self.lbl_prog.configure(text=f"{done}/{total}")
                    self.log(f"  p{page:<4} done  ({done}/{total})")
                elif kind == "done":
                    results = item[1]
                    self._on_done(results)
                elif kind == "error":
                    self._on_error(item[1])
        except queue.Empty:
            pass
        self.root.after(80, self._poll)

    def _on_done(self, results):
        self.btn.configure(state="normal")
        for r in results:
            self.log(f"  p{r.page:<4} {r.width_px}x{r.height_px}px  ->  {r.path.name}")
        self.log(f"Done: {len(results)} image(s) in {self.job.out_dir if self.job else '?'}")
        if messagebox.askyesno(
            "Finished",
            f"Rendered {len(results)} page(s).\n\nOpen the output folder?",
        ):
            self._open_dir(self.job.out_dir if self.job else None)

    def _on_error(self, tb: str):
        self.btn.configure(state="normal")
        self.log("ERROR:\n" + tb)
        messagebox.showerror("Render failed", tb.strip().splitlines()[-1])

    @staticmethod
    def _open_dir(path: Path | None):
        if not path or not path.exists():
            return
        try:
            if sys.platform.startswith("win"):
                os.startfile(path)  # type: ignore[attr-defined]
            elif sys.platform == "darwin":
                subprocess.Popen(["open", str(path)])
            else:
                subprocess.Popen(["xdg-open", str(path)])
        except Exception:
            pass


def main():
    root = tk.Tk()
    try:
        style = ttk.Style(root)
        if "vista" in style.theme_names():
            style.theme_use("vista")
        elif "clam" in style.theme_names():
            style.theme_use("clam")
    except Exception:
        pass
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
