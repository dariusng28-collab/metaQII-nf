#!/usr/bin/env python3
"""Build a self-contained HTML index for a metaQII-nf run.

The report lists the interactive QIIME 2 visualisations (``.qzv``) staged
alongside it — which open at https://view.qiime2.org — and embeds the static
PDF figures produced by the plotting scripts. It uses only the standard
library so it runs inside the QIIME 2 environment without extra dependencies.
"""

import argparse
import datetime as _dt
import html
from pathlib import Path


QZV_DESCRIPTIONS = {
    "demux.qzv": "Per-base sequence quality of the imported reads (use to pick DADA2 truncation lengths).",
    "trimmed-demux.qzv": "Per-base sequence quality after ITSxpress primer/region trimming.",
    "denoising-stats.qzv": "DADA2 denoising statistics (reads in/out per sample).",
    "table.qzv": "Feature-table summary: per-sample frequencies and feature detail.",
    "rep-seqs.qzv": "Representative sequences (ASVs) with lengths and BLAST links.",
    "taxonomy.qzv": "Taxonomic classification assigned to each feature.",
    "taxa-bar-plots.qzv": "Interactive stacked taxonomic bar plots by metadata group.",
    "alpha-rarefaction.qzv": "Alpha-rarefaction curves (richness vs. sampling depth).",
}


def describe(name: str) -> str:
    return QZV_DESCRIPTIONS.get(name, "QIIME 2 visualisation.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", required=True, type=Path, help="Directory to scan and write index.html into.")
    parser.add_argument("--title", default="metaQII-nf report")
    args = parser.parse_args()

    root = args.outdir
    qzvs = sorted(p.name for p in root.glob("*.qzv"))
    pdfs = sorted(root.glob("plots/**/*.pdf"))

    generated = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    title = html.escape(args.title)

    parts = [
        "<!doctype html>",
        "<html lang='en'><head><meta charset='utf-8'>",
        f"<title>{title}</title>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<style>",
        "body{font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;"
        "max-width:960px;margin:2rem auto;padding:0 1rem;color:#1f2937;line-height:1.5}",
        "h1{margin-bottom:.2rem} .meta{color:#6b7280;font-size:.9rem;margin-bottom:1.5rem}",
        "h2{border-bottom:2px solid #e5e7eb;padding-bottom:.3rem;margin-top:2rem}",
        "ul{list-style:none;padding-left:0}",
        "li{margin:.5rem 0;padding:.6rem .8rem;background:#f9fafb;border:1px solid #e5e7eb;border-radius:8px}",
        "a{color:#2563eb;text-decoration:none;font-weight:600} a:hover{text-decoration:underline}",
        ".desc{display:block;color:#4b5563;font-size:.88rem;font-weight:400;margin-top:.15rem}",
        ".note{background:#eff6ff;border:1px solid #bfdbfe;border-radius:8px;padding:.8rem 1rem;font-size:.9rem}",
        "img{max-width:100%;height:auto;border:1px solid #e5e7eb;border-radius:8px;margin:.5rem 0}",
        "figure{margin:1.2rem 0} figcaption{color:#6b7280;font-size:.85rem}",
        "</style></head><body>",
        f"<h1>{title}</h1>",
        f"<div class='meta'>Generated {generated}</div>",
        "<div class='note'>Interactive <code>.qzv</code> files below are QIIME 2 visualisations. "
        "Download one and drag it onto <a href='https://view.qiime2.org' target='_blank' rel='noopener'>"
        "view.qiime2.org</a> — nothing is uploaded, rendering happens in your browser.</div>",
    ]

    parts.append("<h2>Interactive visualisations</h2>")
    if qzvs:
        parts.append("<ul>")
        for name in qzvs:
            parts.append(
                f"<li><a href='{html.escape(name)}' download>{html.escape(name)}</a>"
                f"<span class='desc'>{html.escape(describe(name))}</span></li>"
            )
        parts.append("</ul>")
    else:
        parts.append("<p>No <code>.qzv</code> visualisations were collected.</p>")

    parts.append("<h2>Static figures</h2>")
    if pdfs:
        parts.append("<ul>")
        for pdf in pdfs:
            rel = pdf.relative_to(root).as_posix()
            parts.append(
                f"<li><a href='{html.escape(rel)}' target='_blank' rel='noopener'>"
                f"{html.escape(pdf.stem)}</a>"
                f"<span class='desc'>{html.escape(rel)}</span></li>"
            )
        parts.append("</ul>")
    else:
        parts.append("<p>No PDF figures were produced.</p>")

    parts.append(
        "<h2>Other outputs</h2>"
        "<p>Diversity significance tests (<code>05_diversity/</code>), differential abundance "
        "(<code>07_differential_abundance/</code>), phylogeny (<code>04_phylogeny/</code>) and "
        "FastQC (<code>00_qc/</code>) are written to their numbered folders under the run output "
        "directory.</p>"
    )
    parts.append("</body></html>")

    (root / "index.html").write_text("\n".join(parts), encoding="utf-8")


if __name__ == "__main__":
    main()
