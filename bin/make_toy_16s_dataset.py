#!/usr/bin/env python3

from pathlib import Path
import csv


SAMPLES = {
    "H1": ("Healthy", [("ASV1", 30), ("ASV2", 10)]),
    "H2": ("Healthy", [("ASV1", 28), ("ASV2", 12)]),
    "D1": ("Disease", [("ASV1", 10), ("ASV2", 30)]),
    "D2": ("Disease", [("ASV1", 12), ("ASV2", 28)]),
}

SEQUENCES = {
    "ASV1": "ACGTAGGGTGCGAGCGTTGTCCGGAATTACTGGGCGTAAAGCGTGCGTAGGCGGTTTTTTAAGTCTGATGTGAAAGCCCTCGGCTCAACCGAGGAATTGCATCGGAAACTGGAAA",
    "ASV2": "ACGTAGGGTGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGCGCGCGCAGGCGGTCTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAATTGCATTTGAAACTGGCGAG",
}


def write_fastq(path: Path, reads: list[str]) -> None:
    with path.open("w", newline="\n") as handle:
        for i, seq in enumerate(reads, start=1):
            handle.write(f"@{path.stem}_{i}\n")
            handle.write(f"{seq}\n")
            handle.write("+\n")
            handle.write(f"{'I' * len(seq)}\n")


def main() -> None:
    root = Path(__file__).resolve().parents[1] / "assets" / "smoke_test"
    fastq_dir = root / "fastq"
    fastq_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = root / "manifest.tsv"
    metadata_path = root / "metadata.tsv"

    with manifest_path.open("w", newline="") as manifest_handle:
        writer = csv.writer(manifest_handle, delimiter="\t")
        writer.writerow(["sample-id", "absolute-filepath", "direction"])

        for sample_id, (_, composition) in SAMPLES.items():
            reads = []
            for asv_name, count in composition:
                reads.extend([SEQUENCES[asv_name]] * count)

            fastq_path = (fastq_dir / f"{sample_id}.fastq").resolve()
            write_fastq(fastq_path, reads)
            writer.writerow([sample_id, str(fastq_path), "forward"])

    with metadata_path.open("w", newline="\n") as metadata_handle:
        metadata_handle.write("#SampleID\tstate_of_disease\tcollection_site\n")
        for sample_id, (state, _) in SAMPLES.items():
            metadata_handle.write(f"{sample_id}\t{state}\tSmokeTest\n")


if __name__ == "__main__":
    main()
