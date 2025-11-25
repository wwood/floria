from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

SCRIPT_PATH = Path(__file__).with_name("simulate_badread_long_reads.py")


@pytest.mark.integration
def test_simulation_outputs(tmp_path: Path) -> None:
    reference = tmp_path / "ref.fa"
    reference.write_text(
        ">chr1\n" "ACGTACGTACGT\n",
        encoding="utf-8",
    )

    output_dir = tmp_path / "output"
    command = [
        "python",
        str(SCRIPT_PATH),
        "--reference",
        str(reference),
        "--output-dir",
        str(output_dir),
        "--prefix",
        "sample",
        "--quantity",
        "1k",
        "--length",
        "500,50",
        "--identity",
        "95,5",
        "--snp-density",
        "0.3",
        "--methylation-density",
        "0.5",
        "--seed",
        "42",
    ]

    subprocess.run(command, check=True)

    mutated_fasta = output_dir / "sample.mutated.fa"
    vcf_path = output_dir / "sample.vcf"
    methylation_bed = output_dir / "sample.methylation.bed"
    fastq_path = output_dir / "sample.fastq"

    assert mutated_fasta.exists(), "Mutated FASTA was not created"
    assert vcf_path.exists(), "VCF was not created"
    assert methylation_bed.exists(), "Methylation BED was not created"
    assert fastq_path.exists(), "FASTQ was not created"

    mutated_sequence = "".join(
        line.strip()
        for line in mutated_fasta.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )
    assert mutated_sequence == "AGATACGATAGT"

    vcf_records = [
        line
        for line in vcf_path.read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ]
    assert vcf_records == [
        "chr1\t2\t.\tC\tG\t.\tPASS\tSNP=simulated",
        "chr1\t3\t.\tG\tA\t.\tPASS\tSNP=simulated",
        "chr1\t8\t.\tT\tA\t.\tPASS\tSNP=simulated",
        "chr1\t9\t.\tA\tT\t.\tPASS\tSNP=simulated",
        "chr1\t10\t.\tC\tA\t.\tPASS\tSNP=simulated",
    ]

    methylation_records = methylation_bed.read_text(encoding="utf-8").splitlines()
    assert methylation_records == ["chr1\t5\t6\t5mC"]

    with fastq_path.open("r", encoding="utf-8") as handle:
        first_line = handle.readline()
        assert first_line.startswith("@"), "FASTQ output missing read header"
