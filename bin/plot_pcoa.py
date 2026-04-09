#!/usr/bin/env python3

import argparse
import io
from pathlib import Path

import matplotlib
import pandas as pd
import seaborn as sns
from skbio import OrdinationResults

matplotlib.use("Agg")
import matplotlib.pyplot as plt


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


def main() -> None:
    parser = argparse.ArgumentParser(description="Create PCoA scatter plots from exported QIIME2 ordinations.")
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--metadata-column", required=True)
    parser.add_argument("--pcoa-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    metadata = read_qiime_metadata(args.metadata)
    if args.metadata_column not in metadata.columns:
        raise SystemExit(f"Metadata column '{args.metadata_column}' was not found in {args.metadata}")

    output_dir = args.outdir / "beta_pcoa"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not args.pcoa_dir.exists():
        return

    ordination_dirs = sorted(p for p in args.pcoa_dir.iterdir() if p.is_dir())
    if not ordination_dirs:
        return

    sns.set_theme(style="ticks")

    for ord_dir in ordination_dirs:
        ordination_file = ord_dir / "ordination.txt"
        if not ordination_file.exists():
            continue

        ordination = OrdinationResults.read(str(ordination_file))
        coords = ordination.samples.iloc[:, :2].copy()
        coords.index.name = "SampleID"
        coords = coords.reset_index()

        merged = coords.merge(
            metadata[[args.metadata_column]],
            left_on="SampleID",
            right_index=True,
            how="inner",
        ).dropna()

        if merged.empty:
            continue

        var = ordination.proportion_explained.iloc[:2] * 100
        axis_x = coords.columns[1]
        axis_y = coords.columns[2]

        fig, ax = plt.subplots(figsize=(7, 6))
        sns.scatterplot(
            data=merged,
            x=axis_x,
            y=axis_y,
            hue=args.metadata_column,
            s=70,
            alpha=0.85,
            ax=ax,
        )
        ax.set_xlabel(f"PC1 ({var.iloc[0]:.2f}%)")
        ax.set_ylabel(f"PC2 ({var.iloc[1]:.2f}%)")
        ax.set_title(ord_dir.name.replace("_", " "))
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
        fig.tight_layout()
        fig.savefig(output_dir / f"{ord_dir.name}.pdf")
        plt.close(fig)


if __name__ == "__main__":
    main()
