"""SibeliaZ aligner.

SibeliaZ builds locally collinear blocks and emits a MAF for moderately
divergent genomes (up to ~0.09 substitutions/site). The MAF is projected onto
reference coordinates by :mod:`repgenr.converters.maf_to_fasta`.
"""

from __future__ import annotations

import logging
import sys
from collections.abc import Sequence
from pathlib import Path

from ..converters.maf_to_fasta import maf_to_fasta
from ..core.binaries import BinarySpec
from ..core.bindirs import host_which
from ..core.containers import get_config, run_tool
from ..core.contracts import MSA_FASTA, glob_visible, record_name
from ..core.errors import WorkdirError
from ..core.plugins import ToolCapabilities
from ..core.process import is_gzip, warn_argv_bytes
from .base import Aligner, AlignParams, AlignResult

# SibeliaZ's bash wrapper uses GNU/Linux-only constructs in its alignment step
# (`free`, GNU `find -printf`, `stat -c`, `mktemp --suffix`). On macOS/BSD these
# fail silently, so blocks are found but no MAF is written. We run a patched copy
# of the wrapper on those platforms.
_BSD_PATCHES: tuple[tuple[str, str], ...] = (
    (
        "free -g -w | head -2 | tail -1 | awk '{print $2}'",
        "echo $(( $(sysctl -n hw.memsize) / 1073741824 ))",
    ),
    (
        "free -k -w | head -2 | tail -1 | awk '{print $2}'",
        "echo $(( $(sysctl -n hw.memsize) / 1024 ))",
    ),
    # macOS writes a binary "._" AppleDouble sibling for every file on a non-HFS
    # volume (exFAT, SMB); both finds must skip them, or the binary data lands
    # in the block list and in alignment.maf.
    ('find $outdir -name "*.tmp" -printf "%p\\n"', 'find $outdir -name "*.tmp" ! -name "._*"'),
    ('find $outdir -name "*.msa"  -print0', 'find $outdir -name "*.msa" ! -name "._*" -print0'),
    ('stat -c "%s" $i', "stat -f%z $i"),
    # BSD mktemp has no --suffix, and spoa needs a known FASTA extension: rename
    # the reserved file to <name>.fa, so the bare name is not left behind for
    # every block. Patch the whole $(...) so the result still ends in .fa.
    (
        "$(mktemp --suffix=.fa $outdir/block.XXXXX)",
        '$(t=$(mktemp $outdir/block.XXXXX) && mv "$t" "$t.fa" && echo "$t.fa")',
    ),
    ("ulimit $memory_min", ":"),
)


class SibeliazAligner(Aligner):
    capabilities = ToolCapabilities(
        name="sibeliaz",
        # No pinned BioContainer: the sibeliaz images are BusyBox-based, and
        # the wrapper's `mktemp --suffix` (GNU only) fails there, so the
        # alignment step writes an empty MAF. A Wave image minted from the
        # conda spec carries GNU coreutils and works.
        conda=("bioconda::sibeliaz",),
        required_binaries=(BinarySpec("sibeliaz", version_args=("-v",)),),
        recommended_max_genomes=2000,
        accepted_extras=frozenset({"filtermemory", "kmer", "abundance", "bubble"}),
    )

    def align(
        self,
        genomes: Sequence[Path],
        reference: Path | None,
        out_dir: Path,
        params: AlignParams,
        logger: logging.Logger,
    ) -> AlignResult:
        genomes = list(genomes)
        if reference is None:
            reference = genomes[0]
        out_dir.mkdir(parents=True, exist_ok=True)

        # Default mode runs the block alignment step, producing alignment.maf.
        # (The -n flag would emit only block coordinates, no MAF.)
        sibeliaz_cmd = _sibeliaz_invocation(out_dir, logger)
        # Pass -f (twopaco bloom-filter memory, GB) explicitly so the wrapper
        # skips its system-memory probe, which uses Linux-only `free`/`stat`.
        filtermemory = params.extra.get("filtermemory", 8)
        # Divergence/sensitivity tuning: a smaller -k (k-mer) finds more anchors
        # across divergent genomes; -a (vertex abundance) and -b (max bubble size)
        # trade sensitivity for speed/memory. Unset keys keep SibeliaZ's defaults.
        tuning: list[str] = []
        for key, flag in (("kmer", "-k"), ("abundance", "-a"), ("bubble", "-b")):
            if key in params.extra:
                tuning += [flag, str(params.extra[key])]
        sib_argv: list[str | Path] = [
            *sibeliaz_cmd,
            *tuning,
            "-f",
            str(filtermemory),
            "-t",
            str(params.threads),
            "-o",
            out_dir,
            *[str(g.resolve()) for g in genomes],
        ]
        warn_argv_bytes("sibeliaz", sib_argv, logger)
        run_tool(
            self.capabilities,
            sib_argv,
            logger=logger,
            log_prefix="sibeliaz",
        )
        maf = out_dir / "alignment.maf"
        if not maf.exists():
            # SibeliaZ writes blocks_coords.gff + alignment.maf; locate the MAF
            candidates = glob_visible(out_dir, "*.maf", recursive=True)
            if not candidates:
                raise WorkdirError("SibeliaZ did not produce a MAF file")
            maf = candidates[0]

        # SibeliaZ may find blocks (blocks_coords.gff) yet write a MAF with no
        # alignment rows: its wrapper runs spoa per block with stderr suppressed,
        # so an OOM-killed spoa (global alignment of a very large collinear block
        # in a memory-limited container) is silent and yields an empty MAF.
        if not any(line.startswith("s") for line in maf.read_text().splitlines()):
            gzipped = [g for g in genomes if is_gzip(g)]
            if gzipped:
                # TwoPaCo aborts with "Can't read the input file" on gzip and
                # the wrapper still exits 0, so the cause is not memory.
                raise WorkdirError(
                    f"SibeliaZ wrote an empty MAF at {maf}: it cannot read gzipped "
                    f"FASTA, and {len(gzipped)} input(s) are gzipped (first: "
                    f"{gzipped[0]}). The phylo stage gives it decompressed copies; "
                    "when calling the aligner directly, pass plain FASTA."
                )
            raise WorkdirError(
                f"SibeliaZ wrote an empty MAF (no alignment rows) at {maf}. Its "
                "blocks were found but the per-block spoa alignment produced "
                "nothing -- most often spoa was OOM-killed on a large collinear "
                "block (its global alignment is O(n^2) memory). Give the container "
                "engine more memory (e.g. raise Docker Desktop's VM RAM) or run "
                "this aligner natively."
            )

        # SibeliaZ's MAF uses sequence/contig IDs (FASTA header first token), not
        # genome filenames; build the seqid -> record-name map for the converter.
        name_map = _build_seqid_map(genomes)
        msa = out_dir / MSA_FASTA
        maf_to_fasta(maf, record_name(reference), msa, name_map=name_map)
        return AlignResult(msa_fasta=msa, native_format=maf)


def _sibeliaz_invocation(out_dir: Path, logger: logging.Logger) -> list[str]:
    """Return the command prefix to run SibeliaZ.

    On non-macOS, or when running in a (Linux) container, just ``["sibeliaz"]``.
    On native macOS, write a BSD-compatible copy of the wrapper (its alignment
    step otherwise fails silently) and run it via bash.
    """
    if sys.platform != "darwin" or get_config().active:
        return ["sibeliaz"]

    real = host_which("sibeliaz", "sibeliaz")
    if real is None:
        return ["sibeliaz"]
    script = Path(real).read_text(encoding="utf-8")
    if "-printf" not in script:  # already BSD-friendly / unexpected layout
        return ["sibeliaz"]

    patched = script
    for old, new in _BSD_PATCHES:
        if old in patched:
            patched = patched.replace(old, new)
    patched_path = out_dir / "sibeliaz_bsd.sh"
    patched_path.write_text(patched)
    logger.info("Using BSD-compatible SibeliaZ wrapper on macOS: %s", patched_path)
    return ["bash", str(patched_path)]


def _build_seqid_map(genomes) -> dict[str, str]:
    """Map each FASTA sequence ID (header first token) to its genome's record name."""
    name_map: dict[str, str] = {}
    for genome in genomes:
        stem = record_name(genome)
        with open(genome, encoding="utf-8") as fo:
            for line in fo:
                if line.startswith(">"):
                    seqid = line[1:].split()[0]
                    name_map[seqid] = stem
    return name_map
