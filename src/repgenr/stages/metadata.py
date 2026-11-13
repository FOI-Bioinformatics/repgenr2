"""metadata stage: select a taxon's genomes from a selection source.

The source (``--source``) finds the target's genomes and an outgroup; see
:mod:`repgenr.sources`. Two ship with RepGenR2, both reading GTDB:

* ``tsv`` (default): the full GTDB metadata table of a release.
* ``api``: the GTDB API, for the target taxon only.

Either way the selection is recorded in the SQLite manifest plus ``repgenr.yaml``
provenance.
"""

from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import (
    SELECTION_TSV,
    SelectionRow,
    genome_filename,
    write_selection,
)
from ..core.errors import UserInputError
from ..core.integrity import refuse_foreign_rows
from ..core.manifest import GenomeRecord
from ..core.outgroup_choice import (
    OUTGROUP_MODES,
    read_choice,
    write_choice,
)
from ..core.selection_summary import (
    write_metadata_selected,
    write_selection_summary,
)
from ..sources.base import (
    SOURCE_OPTIONS,
    GenomeSource,
    MetadataParams,
    registry,
    target_species_list,
)

__all__ = ["MetadataParams", "gtdb_provenance", "run", "target_species_list"]


def run(ctx: WorkdirContext, params: MetadataParams) -> int:
    logger = ctx.logger
    refuse_foreign_rows(ctx, "metadata", drop_foreign=params.drop_foreign, logger=logger)
    source = _validate(params)

    # Written again by the source when an outgroup is chosen automatically.
    write_choice(ctx.workdir, None)
    # Candidate genomes before a --limit cut, for selection_stats.json.
    counts: dict[str, int] = {}
    # What another source left for census does not describe this selection.
    for name in registry.names():
        if name != source.name and not registry.is_broken(name):
            for filename in registry.get(name).workdir_files:
                (ctx.workdir / filename).unlink(missing_ok=True)
    provenance = source.prepare(ctx, params, logger)
    selected, outgroup, details = source.select(ctx, params, logger, counts)

    logger.info(
        "Selected %d genomes; outgroup: %s",
        len(selected),
        outgroup.accession if outgroup is not None else "none (--no-outgroup)",
    )
    _populate_manifest(ctx, selected, outgroup)
    _write_outgroup_file(ctx, outgroup.accession if outgroup is not None else None)
    rows = _write_selection(ctx, selected, outgroup)
    write_metadata_selected(ctx.workdir, rows, details)
    write_selection_summary(
        ctx.workdir,
        rows,
        logger,
        stats={
            "source": source.stats_source,
            "candidates": counts.get("candidates", len(selected)),
            "limit": params.limit,
        },
    )

    ctx.config.record_stage(
        "metadata",
        # A web service or a metadata table; no external binary.
        tool=source.record_tool,
        params={
            "source": params.source,
            "release": params.release,
            "api_query_date": provenance.pop("api_query_date", None),
            "version": params.version,
            "dataset": params.dataset,
            "level": params.level,
            "target_family": params.target_family,
            "target_genus": params.target_genus,
            "target_species": params.target_species,
            "selected_count": len(selected),
            "outgroup": outgroup.accession if outgroup is not None else None,
            "no_outgroup": params.no_outgroup,
            "outgroup_mode": params.outgroup_mode,
            # The rule that chose the outgroup: sister, quality (also when
            # sister fell back), first, or named (--outgroup-accession).
            "outgroup_rule": _applied_rule(ctx, params, outgroup),
            "gtdb_tree_path": params.gtdb_tree_path,
            **provenance,
        },
        tool_versions=source.tool_versions or None,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    return len(selected)


def _applied_rule(ctx: WorkdirContext, params: MetadataParams, outgroup) -> str | None:
    if outgroup is None:
        return None
    if params.outgroup_accession:
        return "named"
    rows = read_choice(ctx.workdir)
    return rows[0]["mode"] if rows else "first"


def gtdb_provenance(params: dict) -> dict[str, str]:
    """The GTDB reference behind a metadata record, for ``versions`` and ``status``.

    The table path names its release (``gtdb_release``); the API path has no
    release number and names the UTC date of its query (``gtdb_api_query_date``).
    A record from before the query date was recorded yields nothing for the API.
    """
    if params.get("source") == "ncbi":
        queried = params.get("query_date")
        return {"ncbi_query_date": str(queried)} if queried else {}
    if params.get("source") == "api":
        queried = params.get("api_query_date")
        return {"gtdb_api_query_date": str(queried)} if queried else {}
    release = params.get("release")
    return {"gtdb_release": str(release)} if release else {}


def _validate(params: MetadataParams) -> GenomeSource:
    """Check the parameters; return the selection source they name."""
    names = registry.names()
    if params.source not in names:
        raise UserInputError(f"--source must be one of: {', '.join(names)}.")
    source = registry.create(params.source)
    unread = [
        flag
        for field, (flag, default) in SOURCE_OPTIONS.items()
        if field not in source.options and getattr(params, field) != default
    ]
    if unread:
        raise UserInputError(
            f"{', '.join(unread)} cannot be used with --source {params.source}; "
            "the option belongs to another selection source."
        )
    source.validate(params)
    if not (params.target_genus or params.target_family):
        raise UserInputError("Supply --target-genus or --target-family.")
    if params.level == "species" and not params.target_species:
        raise UserInputError("Level 'species' needs --target-species.")
    species = target_species_list(params.target_species)
    if len(species) > 1 and params.level != "species":
        raise UserInputError(
            f"{len(species)} species were given with -ts; select them with --level species "
            "(at --level genus the whole genus is selected)."
        )
    if len(species) > 1 and not params.target_genus:
        raise UserInputError("Several species (-ts) need their genus (-tg).")
    if params.outgroup_mode not in OUTGROUP_MODES:
        raise UserInputError(
            f"--outgroup-mode must be one of {', '.join(OUTGROUP_MODES)}; "
            f"got {params.outgroup_mode!r}."
        )
    if params.gtdb_tree_path and not Path(params.gtdb_tree_path).expanduser().is_file():
        raise UserInputError(f"--gtdb-tree-path not found: {params.gtdb_tree_path}")
    if params.no_outgroup and params.outgroup_accession:
        raise UserInputError("Give either --outgroup-accession or --no-outgroup, not both.")
    if params.level == "genus" and not params.target_genus:
        raise UserInputError("Level 'genus' needs --target-genus.")
    if params.level == "family" and not params.target_family:
        raise UserInputError("Level 'family' needs --target-family.")
    return source


def _populate_manifest(ctx, selected: list[GenomeRecord], outgroup: GenomeRecord | None) -> None:
    # Replace, not upsert: a re-selection must remove de-selected rows.
    ctx.manifest.replace_genomes([*selected, *([outgroup] if outgroup is not None else [])])


def _write_outgroup_file(ctx, outgroup_acc: str | None) -> None:
    """Name the outgroup; without one (--no-outgroup) remove the file, as
    ingest does, so no earlier outgroup is named."""
    path = ctx.workdir / "outgroup_accession.txt"
    if outgroup_acc is None:
        path.unlink(missing_ok=True)
    else:
        path.write_text(outgroup_acc + "\n")


def _write_selection(
    ctx, selected: list[GenomeRecord], outgroup: GenomeRecord | None
) -> list[SelectionRow]:
    """Publish the portable selection.tsv (the metadata -> genome data-channel hand-off)."""
    rows = []
    for r in (*selected, *([outgroup] if outgroup is not None else [])):
        family, genus, species = r.family or "", r.genus or "", r.species or ""
        rows.append(
            SelectionRow(
                accession=r.accession,
                family=family,
                genus=genus,
                species=species,
                is_outgroup=r.is_outgroup,
                filename=genome_filename(family, genus, species, r.accession),
                completeness=r.completeness,
                contamination=r.contamination,
                gtdb_representative=r.gtdb_representative,
            )
        )
    write_selection(ctx.workdir / SELECTION_TSV, rows)
    return rows
