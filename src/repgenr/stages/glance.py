"""glance stage: quick all-vs-all ANI overview of the genomes.

Ports ``glance.py``: compare all genomes with a dereplicator that implements
``compare`` (dRep compare, Mash primary clustering only, or sourmash sketches
with the ANI estimate used by ``dereplicate --tool sourmash``, c^(1/k) from
the average containment c), copy out the
clustering dendrogram, and plot a box plot + histogram of the pairwise
similarities. ``--tool auto`` prefers dRep and falls back to sourmash.
"""

from __future__ import annotations

import logging
import shutil
from dataclasses import asdict, dataclass, replace
from datetime import UTC, datetime
from pathlib import Path

from ..core.context import WorkdirContext
from ..core.contracts import list_fasta
from ..core.errors import MissingBinaryError, UserInputError, WorkdirError
from ..core.plugins import AUTO
from ..core.process import remove_tree
from ..core.sketches import adapter_sketches, workdir_provider

# File names predate the sourmash backend and stay fixed; the plots' axes name
# the measure the chosen tool reports.
GLANCE_OUTPUTS = (
    "glance_clustering_dendrogram.pdf",
    "glance_MASH_ANI_similarity_boxplot.png",
    "glance_MASH_ANI_similarity_histogram.png",
)

# --tool auto: dRep first (the original glance tool), then sourmash; any other
# compare-capable adapter ranks after these, alphabetically.
_AUTO_PREFERENCE = ("drep", "sourmash")


@dataclass
class GlanceParams:
    threads: int = 16  # same default as the CLI's -t
    tool: str = AUTO
    plot_max: float = 1.0
    plot_min: float = 0.0
    keep_files: bool = False


def resolve_auto_tool(reg=None) -> str | None:
    """The compare-capable dereplicator ``--tool auto`` picks, or None.

    A tool qualifies when it can run here: its binaries are on the PATH, or a
    container backend is active and the adapter declares an image (the same
    availability test ``dereplicate --tool auto`` uses).
    """
    from ..core.plugins import tool_available
    from ..dereplicators.base import compare_supporters
    from ..dereplicators.base import registry as derep_registry

    reg = derep_registry if reg is None else reg

    def rank(name: str) -> tuple[int, str]:
        pref = _AUTO_PREFERENCE.index(name) if name in _AUTO_PREFERENCE else len(_AUTO_PREFERENCE)
        return pref, name

    for name in sorted(compare_supporters(reg), key=rank):
        if tool_available(reg.get(name).capabilities):
            return name
    return None


def no_compare_tool_message(reg=None) -> str:
    from ..core.containers import get_config
    from ..dereplicators.base import compare_supporters

    names = ", ".join(compare_supporters(reg)) or "(none registered)"
    if get_config().active:
        # Under a backend availability means a declared image, not the PATH.
        return (
            "glance --tool auto found no comparison tool to run: the container "
            f"backend is active, and none of {names} declares a container image. "
            "Name one with --tool, or run without --container to use a tool on the PATH."
        )
    return (
        f"glance --tool auto found no comparison tool to run: none of {names} is on "
        "the PATH. Install one of them, use a container backend (--container docker "
        "or singularity), or name one with --tool."
    )


def run(ctx: WorkdirContext, params: GlanceParams) -> Path:
    logger = ctx.logger
    from ..dereplicators.base import compare_supporters, registry

    supporters = compare_supporters()
    if params.tool != AUTO and params.tool not in supporters:
        registry.get(params.tool)  # an unknown or broken plugin reports itself
        raise UserInputError(
            f"Dereplicator '{params.tool}' does not support glance comparisons. "
            f"Tools with compare support: {', '.join(supporters) or 'none'}."
        )
    genomes = list_fasta(ctx.genomes_dir)
    if not genomes:
        raise WorkdirError(f"No genomes under {ctx.genomes_dir}")
    if len(genomes) < 2:
        # dRep compare fails on an empty distance matrix with one genome.
        raise WorkdirError(
            f"glance needs at least two genomes; found {len(genomes)} under {ctx.genomes_dir}"
        )
    if params.tool == AUTO:
        # The CLI resolves auto before the resume fingerprint is taken; this
        # covers direct callers. The record names the concrete tool.
        resolved = resolve_auto_tool()
        if resolved is None:
            raise MissingBinaryError(no_compare_tool_message())
        log_auto_choice(logger, resolved)
        params = replace(params, tool=resolved)
    adapter = registry.create(params.tool)
    versions = adapter.preflight()

    glance_wd = ctx.workdir / "glance_wd"
    if glance_wd.exists():
        remove_tree(glance_wd)

    # An adapter that compares sourmash sketches reads those of the workdir
    # (written first when missing or stale) at its default parameters.
    sketches = adapter_sketches(
        adapter,
        adapter.capabilities.default_params,
        genomes,
        workdir_provider(ctx, logger, params.threads, "glance"),
        logger,
        "glance",
    )
    if sketches is None:
        result = adapter.compare(genomes, glance_wd, params.threads, logger)
    else:
        result = adapter.compare(genomes, glance_wd, params.threads, logger, sketches=sketches)

    # The comparison succeeded: drop the previous outputs before writing new
    # ones, so a plot that is not drawn this time (no similarity in range)
    # does not survive from an earlier run.
    for name in GLANCE_OUTPUTS:
        (ctx.workdir / name).unlink(missing_ok=True)
    out_pdf = ctx.workdir / "glance_clustering_dendrogram.pdf"
    if result.dendrogram is not None:
        shutil.copy2(result.dendrogram, out_pdf)
    else:
        logger.warning(
            "The comparison returned no dendrogram; %s not written, and the next run "
            "repeats the comparison",
            out_pdf.name,
        )

    if result.similarity_csv is not None:
        _plot(result.similarity_csv, ctx.workdir, params, logger, result.measure)
    else:
        logger.warning("No similarity table produced; skipping plots")

    if not params.keep_files and glance_wd.exists():
        remove_tree(glance_wd)
    ctx.config.record_stage(
        "glance",
        tool=params.tool,
        params=asdict(params),
        tool_versions=versions,
        completed=datetime.now(UTC).isoformat(),
    )
    ctx.save_config()
    logger.info("Glance outputs written to %s", ctx.workdir)
    return out_pdf


def log_auto_choice(logger: logging.Logger, tool: str) -> None:
    logger.info("glance --tool auto selected '%s' (dRep preferred, then sourmash)", tool)


def _pair_similarities(mdb: Path, low: float, high: float) -> list[float]:
    """Similarity of each unordered genome pair in the table within [low, high].

    dRep's ``Mdb.csv`` lists every pair in both orders and each genome against itself; a
    pair is counted once (its first readable row) and self-comparisons are
    skipped.
    """
    import csv

    index: dict[str, int] = {}
    seen: set[int] = set()
    values: list[float] = []
    with open(mdb, encoding="utf-8", newline="") as fo:
        for row in csv.DictReader(fo):
            g1, g2 = row.get("genome1"), row.get("genome2")
            if g1 is None or g2 is None or g1 == g2:
                continue
            i = index.setdefault(g1, len(index))
            j = index.setdefault(g2, len(index))
            key = (min(i, j) << 32) | max(i, j)
            if key in seen:
                continue
            try:
                sim = float(row["similarity"])
            except (KeyError, TypeError, ValueError):
                continue
            # Marked only once parsed, so an unreadable row does not hide the
            # valid row of the same pair in the other order.
            seen.add(key)
            if low <= sim <= high:
                values.append(sim)
    return values


def _plot(
    mdb: Path, workdir: Path, params: GlanceParams, logger, measure: str = "MASH ANI"
) -> None:
    from matplotlib import pyplot as plt

    values = _pair_similarities(mdb, params.plot_min, params.plot_max)
    if not values:
        logger.warning("No similarity values in range; skipping plots")
        return

    title = f"{measure}, all-vs-all ({len(values)} genome pairs)"
    fig, ax = plt.subplots()
    ax.boxplot(values)
    ax.set_xticklabels([""])
    ax.set_ylabel(measure)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(workdir / "glance_MASH_ANI_similarity_boxplot.png")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.hist(values, bins=100)
    ax.set_xlabel(measure)
    ax.set_ylabel("Genome pairs")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(workdir / "glance_MASH_ANI_similarity_histogram.png")
    plt.close(fig)
