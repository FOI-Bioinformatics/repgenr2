"""Discrete genome-download step for data-channel orchestration.

The shared-workdir ``genome`` stage reads selected accessions from the SQLite
manifest and downloads them into the working directory. For the data-channel
pipeline the same download is exposed as a stateless step that reads a portable
``selection.tsv`` (produced by the metadata stage) and writes genome FASTAs into
an output directory -- no manifest, no shared workdir.

Output layout under ``out_dir``::

    genomes/                downloaded selected genomes (named by selection.tsv)
    outgroup/               the outgroup genome, if the selection includes one
    missing_accessions.txt  selected accessions NCBI served no genome for
    substituted_accessions.tsv  RefSeq accessions whose genome was taken from
                            the paired GenBank assembly (selected, downloaded)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

from ..core.contracts import (
    MISSING_ACCESSIONS_TXT,
    SUBSTITUTED_ACCESSIONS_TSV,
    read_selection,
)
from ..core.errors import WorkdirError
from ..core.plugins import preflight
from ..core.versions import write_versions_fragment
from .genome import DATASETS_CAPS, download_with_fallback, write_substituted


@dataclass
class GenomeFetchParams:
    selection_tsv: Path
    out_dir: Path
    keep_files: bool = False
    versions_out: Path | None = None
    # Take the paired GenBank assembly when NCBI no longer serves a selected
    # RefSeq accession (as the genome stage does).
    genbank_fallback: bool = True


def genome_fetch(params: GenomeFetchParams, logger: logging.Logger) -> int:
    """Download the genomes listed in ``selection.tsv`` into ``out_dir``.

    Returns the number of selected (non-outgroup) genomes requested.
    """
    if not params.selection_tsv.exists():
        raise WorkdirError(f"genome-fetch: selection file not found: {params.selection_tsv}")
    # Always check the downloader, so a missing datasets binary is reported as
    # such (exit 4) rather than surfacing from the first download.
    versions = preflight(DATASETS_CAPS)
    if params.versions_out is not None:
        write_versions_fragment(params.versions_out, versions)
    rows = read_selection(params.selection_tsv)
    if not rows:
        raise WorkdirError(f"genome-fetch: selection file is empty: {params.selection_tsv}")

    selected = [r for r in rows if not r.is_outgroup]
    outgroup = [r for r in rows if r.is_outgroup]
    scratch = params.out_dir / "scratch"

    missing: list[str] = []
    substituted: dict[str, str] = {}
    if selected:
        filenames = {r.accession: r.filename for r in selected}
        missing, substituted = download_with_fallback(
            list(filenames),
            filenames,
            params.out_dir / "genomes",
            scratch,
            logger,
            params.keep_files,
            genbank_fallback=params.genbank_fallback,
        )
    # The same record the genome stage keeps, so the report and the
    # completeness guard of the data-channel route can excuse these accessions.
    params.out_dir.mkdir(parents=True, exist_ok=True)
    (params.out_dir / MISSING_ACCESSIONS_TXT).write_text(
        "".join(f"{acc}\n" for acc in sorted(missing)), encoding="utf-8"
    )
    write_substituted(params.out_dir / SUBSTITUTED_ACCESSIONS_TSV, substituted)
    if outgroup:
        og_filenames = {r.accession: r.filename for r in outgroup}
        og_missing, og_substituted = download_with_fallback(
            list(og_filenames),
            og_filenames,
            params.out_dir / "outgroup",
            scratch,
            logger,
            params.keep_files,
            genbank_fallback=params.genbank_fallback,
        )
        substituted.update(og_substituted)
        write_substituted(params.out_dir / SUBSTITUTED_ACCESSIONS_TSV, substituted)
        if og_missing:
            # A selection that names an outgroup must not continue without it:
            # the tree would be built and rooted differently with no notice.
            raise WorkdirError(
                f"NCBI serves no assembly for the outgroup {', '.join(og_missing)}; "
                "choose another with 'repgenr metadata --outgroup-accession'."
            )

    logger.info(
        "genome-fetch: requested %d genomes + %d outgroup into %s (%d not served by NCBI, "
        "%d taken from the paired GenBank assembly)",
        len(selected),
        len(outgroup),
        params.out_dir,
        len(missing),
        len(substituted),
    )
    return len(selected)
