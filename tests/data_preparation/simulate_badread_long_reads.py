"""Simulate long reads for tests using badread.

This script builds a closely related genome from a reference by introducing SNPs
and selecting methylated cytosines at user-specified densities. The mutated
reference is then used as input for ``badread simulate`` to generate long-read
fastq files.

Example:
    python simulate_badread_long_reads.py \\
        --reference tests/MN-03.fa \\
        --output-dir tests/data_preparation/output \\
        --prefix mn03_sim \\
        --quantity 20x \\
        --snp-density 0.001 \\
        --methylation-density 0.0005 \\
        --seed 42

Requirements:
    * badread (https://github.com/rrwick/Badread)
"""

from __future__ import annotations

import argparse
import random
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

BASES = ("A", "C", "G", "T")


def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    """Read a FASTA file into a list of (header, sequence) tuples."""

    records: List[Tuple[str, str]] = []
    header: str | None = None
    seq_parts: List[str] = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts)))

    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    return records


def mutate_sequence(
    chrom: str,
    sequence: str,
    rng: random.Random,
    snp_density: float,
    methylation_density: float,
    methylation_label: str,
) -> Tuple[str, List[str], List[str]]:
    """Mutate bases and collect SNP and methylation annotations.

    Returns mutated sequence, VCF records (strings), and methylation BED lines.
    """

    seq_list = list(sequence.upper())
    vcf_records: List[str] = []
    methylation_records: List[str] = []

    for idx, base in enumerate(seq_list):
        if base not in BASES:
            continue

        # SNP mutations
        if rng.random() < snp_density:
            alternatives = [b for b in BASES if b != base]
            alt = rng.choice(alternatives)
            seq_list[idx] = alt
            pos = idx + 1  # VCF is 1-based
            vcf_records.append(
                f"{chrom}\t{pos}\t.\t{base}\t{alt}\t.\tPASS\tSNP=simulated"
            )
        else:
            alt = base

        # Methylation annotations on cytosines of the mutated sequence
        if alt == "C" and rng.random() < methylation_density:
            start = idx
            end = idx + 1
            methylation_records.append(
                f"{chrom}\t{start}\t{end}\t{methylation_label}"
            )

    return "".join(seq_list), vcf_records, methylation_records


def write_fasta(records: Iterable[Tuple[str, str]], path: Path) -> None:
    with path.open("w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def write_vcf(
    vcf_records: Dict[str, List[str]],
    path: Path,
    reference_path: Path,
) -> None:
    with path.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(f"##reference={reference_path}\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom in sorted(vcf_records):
            for record in vcf_records[chrom]:
                handle.write(record + "\n")


def write_bed(methylation_records: Dict[str, List[str]], path: Path) -> None:
    with path.open("w") as handle:
        for chrom in sorted(methylation_records):
            for record in methylation_records[chrom]:
                handle.write(record + "\n")


def simulate_reads_with_badread(
    reference: Path,
    output_fastq: Path,
    quantity: str,
    length: str,
    identity: str,
    error_model: str,
    qscore_model: str,
    seed: int,
) -> None:
    cmd = [
        "badread",
        "simulate",
        "--reference",
        str(reference),
        "--quantity",
        quantity,
        "--length",
        length,
        "--identity",
        identity,
        "--error_model",
        error_model,
        "--qscore_model",
        qscore_model,
        "--seed",
        str(seed),
    ]
    subprocess.run(cmd, check=True, stdout=output_fastq.open("wb"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Simulate long reads that reflect SNPs and methylated positions "
            "in a genome closely related to the input reference."
        )
    )
    parser.add_argument(
        "--reference",
        required=True,
        type=Path,
        help="Reference FASTA file used as the template genome.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory where outputs will be written.",
    )
    parser.add_argument(
        "--prefix",
        required=True,
        help="Prefix for all generated files.",
    )
    parser.add_argument(
        "--quantity",
        default="20x",
        help="Total quantity passed to badread (e.g. '20x' or '1G').",
    )
    parser.add_argument(
        "--length",
        default="10000,4000",
        help="Read length distribution for badread (mean,sd).",
    )
    parser.add_argument(
        "--identity",
        default="95,5",
        help="Identity distribution for badread (mean,sd).",
    )
    parser.add_argument(
        "--error-model",
        default="nanopore2023",
        help="Badread error model to use (e.g. 'nanopore2023').",
    )
    parser.add_argument(
        "--qscore-model",
        default="nanopore2023",
        help="Badread qscore model to use (e.g. 'nanopore2023').",
    )
    parser.add_argument(
        "--snp-density",
        type=float,
        default=0.001,
        help="Probability of a SNP at each base of the reference genome.",
    )
    parser.add_argument(
        "--methylation-density",
        type=float,
        default=0.0005,
        help="Probability of methylation at each cytosine of the mutated genome.",
    )
    parser.add_argument(
        "--methylation-label",
        default="5mC",
        help="Label written to the methylation BED file.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed for SNP placement and badread.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(args.seed)

    records = parse_fasta(args.reference)

    mutated_records: List[Tuple[str, str]] = []
    vcf_records: Dict[str, List[str]] = {}
    methylation_records: Dict[str, List[str]] = {}

    for header, sequence in records:
        mutated_seq, chrom_vcf, chrom_methylation = mutate_sequence(
            header,
            sequence,
            rng,
            args.snp_density,
            args.methylation_density,
            args.methylation_label,
        )
        mutated_records.append((header, mutated_seq))
        if chrom_vcf:
            vcf_records.setdefault(header, []).extend(chrom_vcf)
        if chrom_methylation:
            methylation_records.setdefault(header, []).extend(chrom_methylation)

    mutated_fasta = args.output_dir / f"{args.prefix}.mutated.fa"
    vcf_path = args.output_dir / f"{args.prefix}.vcf"
    methylation_bed = args.output_dir / f"{args.prefix}.methylation.bed"
    fastq_output = args.output_dir / f"{args.prefix}.fastq"

    write_fasta(mutated_records, mutated_fasta)
    write_vcf(vcf_records, vcf_path, args.reference)
    write_bed(methylation_records, methylation_bed)

    simulate_reads_with_badread(
        mutated_fasta,
        fastq_output,
        args.quantity,
        args.length,
        args.identity,
        args.error_model,
        args.qscore_model,
        args.seed,
    )


if __name__ == "__main__":
    main()
