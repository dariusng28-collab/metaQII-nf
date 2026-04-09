#!/usr/bin/env python3

import argparse
import io
from pathlib import Path

import matplotlib
from biom import load_table
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt


LEVELS = [
    ("Phylum", 1),
    ("Class", 2),
    ("Order", 3),
    ("Family", 4),
    ("Genus", 5),
    ("Species", 6),
]


def read_qiime_metadata(path: Path) -> pd.DataFrame:
    lines = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#q2:"):
                continue
            lines.append(line)

    frame = pd.read_csv(io.StringIO("".join(lines)), sep="\t", dtype=str)
    first = frame.columns[0]
    frame = frame.rename(columns={first: first.lstrip("#")})
    sample_col = frame.columns[0]
    return frame.set_index(sample_col, drop=False)


def clean_taxon(raw_taxon: str, level_index: int) -> str:
    if pd.isna(raw_taxon):
        return "Unassigned"
    parts = [part.strip() for part in str(raw_taxon).split(";")]
    if len(parts) <= level_index:
        return "Unassigned"
    token = parts[level_index]
    token = token.split("__", 1)[-1].strip()
    return token if token else "Unassigned"


def main() -> None:
    parser = argparse.ArgumentParser(description="Create stacked taxonomic bar plots from exported QIIME2 data.")
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--metadata-column", required=True)
    parser.add_argument("--biom", required=True, type=Path)
    parser.add_argument("--taxonomy", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--top-n", type=int, default=15)
    args = parser.parse_args()

    metadata = read_qiime_metadata(args.metadata)
    if args.metadata_column not in metadata.columns:
        raise SystemExit(f"Metadata column '{args.metadata_column}' was not found in {args.metadata}")

    table = load_table(str(args.biom)).to_dataframe(dense=True)
    table = table.loc[:, [sample for sample in table.columns if sample in metadata.index]]
    if table.empty:
        return

    relative = table.div(table.sum(axis=0), axis=1).fillna(0)
    taxonomy = pd.read_csv(args.taxonomy, sep="\t")
    taxonomy = taxonomy.rename(columns={"Feature ID": "FeatureID"})
    taxonomy = taxonomy.set_index("FeatureID")

    plots_dir = args.outdir / "taxa_barplots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    color_map = plt.get_cmap("tab20")

    for level_name, level_index in LEVELS:
        labels = taxonomy.reindex(relative.index)["Taxon"].map(lambda tax: clean_taxon(tax, level_index)).fillna("Unassigned")
        collapsed = relative.groupby(labels).sum()

        grouped = collapsed.T.join(metadata[[args.metadata_column]], how="inner")
        grouped = grouped.groupby(args.metadata_column).mean()
        if grouped.empty:
            continue

        mean_abundance = grouped.mean(axis=0).sort_values(ascending=False)
        top_taxa = mean_abundance.head(args.top_n).index.tolist()
        other = grouped.drop(columns=top_taxa, errors="ignore").sum(axis=1)
        grouped = grouped[top_taxa]
        grouped["Other"] = other

        fig, ax = plt.subplots(figsize=(10, 6))
        grouped.plot(
            kind="bar",
            stacked=True,
            ax=ax,
            color=[color_map(i % color_map.N) for i in range(grouped.shape[1])],
            width=0.8,
        )
        ax.set_xlabel(args.metadata_column)
        ax.set_ylabel("Mean relative abundance")
        ax.set_title(f"{level_name} composition by {args.metadata_column}")
        ax.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            title=level_name,
            frameon=False,
        )
        fig.tight_layout()
        fig.savefig(plots_dir / f"{level_name.lower()}_taxonomy.pdf")
        plt.close(fig)


if __name__ == "__main__":
    main()
