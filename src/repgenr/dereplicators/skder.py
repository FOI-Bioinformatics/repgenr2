"""skDER dereplication adapter (skani-based).

skDER scales to tens of thousands of genomes natively (no chunking needed). It
estimates pairwise ANI/AF with skani and picks representatives in one round.

CLI used::

    skder -g <staged_dir> -o <out_dir> -i <ANI%> -f <AF%> -c <threads> \
          -d dynamic|greedy|low_mem_greedy [--skani-triangle-parameters="-s <screen%>"]

skani's screen (``-s``) skips the alignment of pairs whose k-mer sketches put
them below the given identity. skDER's own default is the ANI cutoff minus 10,
so at a cutoff of 99.95 every pair above about 90 percent is aligned, which in
a genus-level set is most same-species pairs. The adapter sets the screen to
the cutoff minus 1 when the cutoff is at or above 99 (``_SCREEN_MARGIN_PCT``,
``_SCREEN_MIN_CUTOFF_PCT``), and leaves skDER's default below that, where the
sketch estimate is less exact. ``--tool-arg screen=<percent>`` overrides it
and ``screen=skder`` restores skDER's own default; a screen above the cutoff
is refused, because pairs at the cutoff would be skipped. The dereplicate
stage checks the override against the final cutoff before any work; a call
at a lower cutoff (a ``--pre-*`` chunk, a ``--target-reps`` step) falls back
to the default rule for that call. The screen is a ``skani triangle`` option,
so it applies to the greedy and dynamic modes; ``low_mem_greedy`` runs
``skani search`` instead and does not read it.

Three robustness notes:

* skDER takes the genome set as one directory, not one path per genome: the
  adapter stages a directory of symbolic links (one per genome, same basename)
  and passes that single path to ``-g``. Thousands of paths on the command line
  can exceed the operating system's argument-size limit (ARG_MAX), and skDER
  has no file-of-filenames option. skDER lists the directory and keeps every
  entry with a FASTA suffix, so the directory holds the links and nothing else.
  It lives in the local scratch directory below, not in the working
  directory: exFAT cannot hold symbolic links, and macOS adds ``._*`` files
  there, which skDER would read as genomes.

* skDER is run in a local scratch directory and the results are copied into the
  working directory afterwards. skDER's ``determineN50`` does
  ``os.listdir(N50/)`` and parses every file as a result row, which breaks on
  exFAT/NTFS volumes where macOS leaves ``._*`` AppleDouble companions. Keeping
  skDER's output on the local (APFS/ext4) temp filesystem avoids that.
* Plain skDER reports representatives only (membership is a CiDDER feature), so
  non-representatives are assigned to their closest representative using skDER's
  ``Skani_Triangle_Edge_Output.txt`` pairwise ANI/AF table.
"""

from __future__ import annotations

import csv
import logging
import os
import shutil
import tempfile
from collections.abc import Iterator, Sequence
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.containers import run_tool
from ..core.contracts import FASTA_SUFFIXES, list_fasta
from ..core.errors import ToolExecutionError, UserInputError, WorkdirError
from ..core.logging import capped_names, console_extra
from ..core.plugins import ToolCapabilities
from ..core.process import link_or_copy, remove_tree
from .base import (
    STATUS_CONTAINED,
    STATUS_REPRESENTATIVE,
    Dereplicator,
    DerepParams,
    DerepResult,
)

# Placement warnings beyond this count are folded into one line.
_MAX_PLACED_LINES = 5
# skDER asks for confirmation below this ANI (percent); skani is unreliable there.
_MIN_ANI_PCT = 80.0
# Default skani screen: cutoff minus this margin, once the cutoff is at or above
# the minimum. skani's screen is a k-mer estimate of identity made before any
# alignment; its own help calls it approximate. The margin of one point at
# cutoffs of 99 and above is a judgment, not a measured bound: at that identity
# bacterial genomes share nearly all sketch k-mers, so the estimate sits close
# to the aligned value, while at lower identities it is looser and skDER's own
# default (cutoff minus 10) is kept.
_SCREEN_MARGIN_PCT = 1.0
_SCREEN_MIN_CUTOFF_PCT = 99.0
# Override-above-cutoff fallbacks already reported, so the log says it once.
_screen_fallbacks_logged: set[tuple[float, float]] = set()
# Screens already reported, so the log names the value and its source once per run.
_screens_logged: set[tuple[float | None, str]] = set()


class SkderDereplicator(Dereplicator):
    capabilities = ToolCapabilities(
        name="skder",
        required_binaries=(
            BinarySpec("skder", version_args=("-v",), min_version="1.0"),
            BinarySpec("skani", version_args=("-V",), min_version="0.2"),
        ),
        default_params={"mode": "greedy", "screen": None},
        recommended_max_genomes=None,  # scales natively
        supports_native_scaling=True,
        # No pinned BioContainer: the skder images (1.3.2 through 1.3.8) are
        # BusyBox-based and their `sort` rejects the `--parallel` flag skder's
        # greedy mode runs, so that mode silently returns no representatives
        # there (dynamic mode, which uses skDERcore, is unaffected). A Wave
        # image minted from the conda spec carries GNU coreutils and works.
        conda=("bioconda::skder",),
        accepted_extras=frozenset({"mode", "screen"}),
        # skDER clusters at one ANI cutoff; there is no primary pre-clustering.
        ignored_params=frozenset({"primary_ani"}),
    )

    def dereplicate(
        self,
        genomes: Sequence[Path],
        out_dir: Path,
        params: DerepParams,
        logger: logging.Logger,
    ) -> DerepResult:
        check_secondary_ani(params.secondary_ani)
        ani_pct = _as_percent(params.secondary_ani)
        out_dir.mkdir(parents=True, exist_ok=True)

        # Run skDER on a local temp filesystem to avoid exFAT/NTFS ._* breakage,
        # then copy its results back next to the working directory.
        local_tmp = Path(tempfile.mkdtemp(prefix="repgenr_skder_"))
        result_dir = local_tmp / "skder_out"
        mode = params.extra.get("mode", self.capabilities.default_params["mode"])
        af_pct = _as_percent(params.aligned_fraction)
        requested = params.extra.get("screen")
        screen = resolve_screen(requested, float(ani_pct), logger)
        source = "--tool-arg screen" if requested is not None else "default rule"
        if (screen, source) not in _screens_logged:
            _screens_logged.add((screen, source))
            logger.info(
                "skani screen %s (%s)",
                f"{screen:g} percent" if screen is not None else "left to skDER",
                source,
            )
        try:
            genome_dir = stage_genome_dir(genomes, local_tmp / "genomes")
        except BaseException:
            shutil.rmtree(local_tmp, ignore_errors=True)
            raise
        # The links point at the genome files' own directories; bind those so
        # a container backend resolves them.
        genome_dirs = sorted({os.path.dirname(os.path.abspath(g)) for g in genomes})
        cmd = [
            "skder",
            "-g",
            genome_dir,
            "-o",
            result_dir,
            "-i",
            ani_pct,
            "-f",
            af_pct,
            "-c",
            params.threads,
            "-d",
            mode,
        ]
        if screen is not None:
            # skDER takes the skani triangle options as one string. Its parser
            # reads a separate "-s ..." value as its own -s flag, so the value
            # is attached with "=", which argparse resolves before that.
            cmd.append(f"--skani-triangle-parameters=-s {screen:g}")
        try:
            run_tool(
                self.capabilities,
                cmd,
                logger=logger,
                log_prefix="skder",
                extra_mounts=genome_dirs,
            )
            staged = out_dir / "skder_out"
            if staged.exists():
                remove_tree(staged)
            # Hardlink the result files when possible (no extra disk; the temp dir
            # is removed below, but hardlinks keep the data alive). Falls back to a
            # copy across filesystems (temp vs workdir) or on exFAT.
            shutil.copytree(result_dir, staged, copy_function=link_or_copy)
        finally:
            shutil.rmtree(local_tmp, ignore_errors=True)

        return _parse_skder_output(staged, genomes, float(ani_pct), float(af_pct), logger)


def reset_screen_log() -> None:
    """Forget the screens and fallbacks reported, so the next run logs them again."""
    _screen_fallbacks_logged.clear()
    _screens_logged.clear()


def default_screen(ani_pct: float) -> float | None:
    """Cutoff minus the margin at or above the minimum cutoff, else None (skDER's own)."""
    if ani_pct >= _SCREEN_MIN_CUTOFF_PCT:
        return round(ani_pct - _SCREEN_MARGIN_PCT, 4)
    return None


def resolve_screen(
    requested: str | float | None,
    ani_pct: float,
    logger: logging.Logger | None = None,
) -> float | None:
    """The skani screen (percent) to pass, or None for skDER's own default.

    ``requested`` is the ``--tool-arg screen=`` value when given. Without it the
    screen is ``default_screen``. A requested screen above this call's cutoff
    is refused when no logger is given (the stage precheck, against the final
    cutoff); with a logger, the call keeps the default rule instead and says so
    once, because a chunk at ``--pre-secondary-ani`` or a ``--target-reps`` step
    runs below the final cutoff by design.
    """
    if requested is None:
        return default_screen(ani_pct)
    if isinstance(requested, str) and requested.strip().lower() == "skder":
        return None
    try:
        screen = float(requested)
    except (TypeError, ValueError):
        raise UserInputError(
            f"--tool-arg screen= must be an ANI percentage or 'skder', got {requested!r}."
        ) from None
    if not 0 <= screen < 100:
        raise UserInputError(f"--tool-arg screen= must be between 0 and 100, got {screen:g}.")
    if screen > ani_pct:
        if logger is None:
            raise UserInputError(
                f"--tool-arg screen={screen:g} is above the ANI cutoff {ani_pct:g} percent; "
                "pairs at the cutoff would be skipped before alignment."
            )
        fallback = default_screen(ani_pct)
        if (screen, ani_pct) not in _screen_fallbacks_logged:
            _screen_fallbacks_logged.add((screen, ani_pct))
            logger.warning(
                "--tool-arg screen=%g is above this call's ANI cutoff %g percent; "
                "using %s for calls at this cutoff.",
                screen,
                ani_pct,
                f"a screen of {fallback:g}" if fallback is not None else "skDER's own screen",
            )
        return fallback
    return screen


def stage_genome_dir(genomes: Sequence[Path], dest: Path) -> Path:
    """Stage ``genomes`` as a directory of entries named by their basenames.

    Each entry is a symbolic link to the absolute (not symlink-resolved) path
    of the genome, so a container backend that binds the genome directories
    at their host paths resolves it. Where the filesystem refuses a symbolic
    link the genome is hardlinked, or copied as a last resort
    (:func:`~repgenr.core.process.link_or_copy`). skDER reports genomes by
    basename, so two genomes with the same basename cannot be told apart;
    they are refused before anything is staged.
    """
    seen: dict[str, Path] = {}
    for g in genomes:
        other = seen.setdefault(Path(g).name, Path(g))
        if other != Path(g):
            raise UserInputError(
                f"Two genomes share the basename {Path(g).name!r} ({other} and {g}); "
                "skDER reports genomes by basename, so the results could not be "
                "assigned. Rename one of them."
            )
    dest.mkdir(parents=True, exist_ok=False)
    for name, src in seen.items():
        target = os.path.abspath(src)
        link = dest / name
        try:
            os.symlink(target, link)
        except OSError:
            link_or_copy(target, link)
    return dest


def check_secondary_ani(value: float, label: str = "--secondary-ani") -> None:
    """Refuse an ANI cutoff below skDER's floor (80 percent).

    skDER stops at an interactive question below it (skani's estimates are
    unreliable there). The dereplicate stage calls this from its precheck, so
    a refused rerun leaves the finished record intact; the adapter calls it
    again as a backstop for callers that bypass the stage.
    """
    if float(_as_percent(value)) < _MIN_ANI_PCT:
        raise UserInputError(
            f"skDER does not cluster below {_MIN_ANI_PCT:g} percent ANI "
            f"({label} {value:g}); skani's estimates are unreliable there. "
            "Use a higher value or --tool sourmash."
        )


def _as_percent(value: float) -> str:
    """skDER expects identity/fraction as a percentage (e.g. 99.0, 50.0)."""
    pct = value * 100 if value <= 1.0 else value
    return f"{pct:g}"


def _parse_skder_output(
    result_dir: Path,
    genomes: Sequence[Path],
    ani_cutoff: float,
    af_cutoff: float,
    logger: logging.Logger,
) -> DerepResult:
    rep_dir = _find_representatives_dir(result_dir)
    if rep_dir is None:
        raise WorkdirError(
            f"Could not locate skDER representative genomes under {result_dir}. "
            "Confirm skDER ran successfully."
        )
    representatives = list_fasta(rep_dir)
    if not representatives:
        raise ToolExecutionError(["skder"], 1, "skDER produced no representative genomes")

    rep_names = {p.name for p in representatives}
    clusters: dict[str, list[str]] = {p.name: [] for p in representatives}
    status: dict[str, str] = {p.name: STATUS_REPRESENTATIVE for p in representatives}

    # Assign each non-representative to its closest representative via the skani
    # edge table (best ANI above the cutoffs). Streamed, so the full edge list is
    # never materialized -- the table can be large for big inputs. As in skDER's
    # own selection, the aligned-fraction cutoff applies to the member's aligned
    # fraction: a partial genome is covered by a complete representative even
    # though only a small fraction of the representative aligns to it.
    best: dict[str, tuple[str, float]] = {}  # member -> (rep, ani)
    for a, b, ani, af_a, af_b in _iter_edges(result_dir):
        if ani < ani_cutoff:
            continue
        for member, rep, member_af in ((a, b, af_a), (b, a, af_b)):
            if member_af < af_cutoff:
                continue
            if rep in rep_names and member not in rep_names:
                if member not in best or ani > best[member][1]:
                    best[member] = (rep, ani)
    # skDER covered every genome it did not select, but its aligned-fraction
    # rule may differ from the one above at the margins (rounding at the
    # cutoff, a later skDER rule). Place such a genome under its highest-ANI
    # representative at or above the ANI cutoff instead of failing the whole
    # stage, and name it.
    unassigned = {g.name for g in genomes} - rep_names - best.keys()
    if unassigned:
        nearest: dict[str, tuple[str, float]] = {}
        for a, b, ani, _af_a, _af_b in _iter_edges(result_dir):
            if ani < ani_cutoff:
                continue
            for member, rep in ((a, b), (b, a)):
                if member in unassigned and rep in rep_names:
                    if member not in nearest or ani > nearest[member][1]:
                        nearest[member] = (rep, ani)
        _warn_placed(nearest, logger)
        best.update(nearest)

    for member, (rep, _ani) in best.items():
        clusters[rep].append(member)
        status[member] = STATUS_CONTAINED

    # any remaining input genome that is neither a representative nor assigned;
    # the stage's completeness check refuses the result and names it
    for g in genomes:
        status.setdefault(g.name, STATUS_CONTAINED)

    return DerepResult(
        representatives=representatives,
        clusters=clusters,
        genome_status=status,
    )


def _warn_placed(nearest: dict[str, tuple[str, float]], logger: logging.Logger) -> None:
    """Name the genomes placed under their closest representative.

    One line per genome up to ``_MAX_PLACED_LINES``; beyond that one line, of
    which the console shows the first names and a count and the run log keeps
    every placement.
    """
    placed = sorted(nearest.items())
    if len(placed) <= _MAX_PLACED_LINES:
        for member, (rep, ani) in placed:
            logger.warning(
                "skDER did not select %s, and no edge to a representative meets "
                "the aligned-fraction cutoff; placing it under its closest "
                "representative %s (ANI %.2f).",
                member,
                rep,
                ani,
            )
        return
    head = (
        f"skDER did not select {len(placed)} genomes, and no edge to a representative "
        "meets the aligned-fraction cutoff; placing each under its closest representative"
    )
    full = "; ".join(f"{member} under {rep} (ANI {ani:.2f})" for member, (rep, ani) in placed)
    logger.warning(
        "%s: %s.",
        head,
        full,
        extra=console_extra(f"{head}: {capped_names([m for m, _ in placed])}."),
    )


def _find_representatives_dir(result_dir: Path) -> Path | None:
    candidates = [
        result_dir / "Dereplicated_Representative_Genomes",
        result_dir / "skDER_Dereplicated_Representative_Genomes",
    ]
    for cand in candidates:
        if cand.is_dir():
            return cand
    for sub in result_dir.rglob("*"):
        if sub.is_dir() and ("epresentative" in sub.name or "Dereplicated" in sub.name):
            if any(c.name.endswith(FASTA_SUFFIXES) for c in sub.iterdir()):
                return sub
    return None


def _iter_edges(result_dir: Path) -> Iterator[tuple[str, str, float, float, float]]:
    """Stream skDER's skani edge table as (a, b, ANI, AF_a, AF_b) tuples.

    The table is skani's ``Ref_file Query_file ANI Align_fraction_ref
    Align_fraction_query ...``; ``AF_a`` is the aligned fraction of genome
    ``a`` (the reference), ``AF_b`` that of ``b``. Names are basenames.
    """
    edge_file = result_dir / "Skani_Triangle_Edge_Output.txt"
    if not edge_file.exists():
        matches = sorted(result_dir.rglob("*Edge_Output*.txt"))
        if not matches:
            return
        edge_file = matches[0]

    with open(edge_file, encoding="utf-8", newline="") as fo:
        for ln, fields in enumerate(csv.reader(fo, delimiter="\t")):
            if ln == 0:
                continue  # header
            if len(fields) < 5:
                continue
            try:
                ani = float(fields[2])
                af_a, af_b = float(fields[3]), float(fields[4])
            except ValueError:
                continue
            yield (Path(fields[0]).name, Path(fields[1]).name, ani, af_a, af_b)
