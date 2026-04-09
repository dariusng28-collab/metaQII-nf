#!/usr/bin/env python3

import argparse
import io
from pathlib import Path

import matplotlib
import pandas as pd
import seaborn as sns

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


def read_alpha_vector(vector_path: Path) -> pd.DataFrame:
    frame = pd.read_csv(vector_path, sep="\t", comment="#")
    first = frame.columns[0]
    second = frame.columns[1]
    frame = frame.rename(columns={first: "SampleID", second: "value"})
    return frame[["SampleID", "value"]]


def metric_label(metric_dir_name: str) -> str:
    return metric_dir_name.replace("_vector", "").replace("_", " ").title()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot alpha diversity boxplots from exported QIIME2 vectors.")
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--metadata-column", required=True)
    parser.add_argument("--alpha-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    metadata = read_qiime_metadata(args.metadata)
    if args.metadata_column not in metadata.columns:
        raise SystemExit(f"Metadata column '{args.metadata_column}' was not found in {args.metadata}")

    output_dir = args.outdir / "alpha_boxplots"
    output_dir.mkdir(parents=True, exist_ok=True)

    if not args.alpha_dir.exists():
        return

    metric_dirs = sorted(p for p in args.alpha_dir.iterdir() if p.is_dir())
    if not metric_dirs:
        return

    sns.set_theme(style="whitegrid")

    for metric_dir in metric_dirs:
        vector_file = metric_dir / "alpha-diversity.tsv"
        if not vector_file.exists():
            continue

        alpha_frame = read_alpha_vector(vector_file)
        merged = alpha_frame.merge(
            metadata[[args.metadata_column]],
            left_on="SampleID",
            right_index=True,
            how="inner",
        ).dropna()

        if merged.empty:
            continue

        metric_name = metric_label(metric_dir.name)
        fig, ax = plt.subplots(figsize=(8, 5))
        sns.boxplot(
            data=merged,
            x=args.metadata_column,
            y="value",
            ax=ax,
            color="#99c1f1",
            fliersize=0,
        )
        sns.stripplot(
            data=merged,
            x=args.metadata_column,
            y="value",
            ax=ax,
            color="#1f2937",
            alpha=0.7,
            size=4,
        )

        ax.set_xlabel(args.metadata_column)
        ax.set_ylabel(metric_name)
        ax.set_title(f"{metric_name} by {args.metadata_column}")
        ax.tick_params(axis="x", rotation=30)
        fig.tight_layout()
        fig.savefig(output_dir / f"{metric_dir.name}.pdf")
        plt.close(fig)


if __name__ == "__main__":
    main()
