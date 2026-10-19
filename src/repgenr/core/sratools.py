"""The SRA toolkit specification: ``fastq-dump`` for streaming runs.

``repgenr reads-sketch`` streams the reads of a run from NCBI with
``fastq-dump -Z`` and pipes them into sourmash, so no FASTQ file is written.
``fasterq-dump`` is not used: it has no spot cap (``-X``), which the command
needs to read only the first spots of a run.
"""

from __future__ import annotations

from .binaries import BinarySpec
from .plugins import ToolCapabilities

SRATOOLS_IMAGE = "quay.io/biocontainers/sra-tools:3.4.1--h4304569_1"
SRATOOLS_CONDA: tuple[str, ...] = ("bioconda::sra-tools",)
FASTQ_DUMP_BINARY = BinarySpec("fastq-dump", version_args=("--version",), min_version="3.0")

SRATOOLS_TOOL = ToolCapabilities(
    name="sra-tools",
    required_binaries=(FASTQ_DUMP_BINARY,),
    container=SRATOOLS_IMAGE,
    conda=SRATOOLS_CONDA,
)


def fastq_dump_command(accession: str, max_spots: int | None) -> list[str]:
    """The ``fastq-dump`` call that streams a run's biological reads to stdout.

    ``--split-spot`` writes both mates of a pair as separate records (mate
    order does not matter to a k-mer sketch) and ``--skip-technical`` drops
    barcode and adapter reads. ``-X`` caps the stream at the first
    ``max_spots`` spots; a paired-end spot gives two reads.
    """
    cmd = ["fastq-dump", "-Z", "--split-spot", "--skip-technical"]
    if max_spots is not None:
        cmd += ["-X", str(max_spots)]
    cmd.append(accession)
    return cmd
