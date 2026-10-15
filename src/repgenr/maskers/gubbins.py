"""Gubbins recombination masking of a whole-genome alignment."""

from __future__ import annotations

import logging
import re
import shlex
import shutil
from pathlib import Path

from ..core.binaries import BinarySpec
from ..core.bindirs import host_which, tool_path
from ..core.containers import available_cpus, run_tool, runs_on_host
from ..core.contracts import parse_genome_filename
from ..core.errors import ToolExecutionError, UserInputError, WorkdirError
from ..core.plugins import ToolCapabilities
from .base import Masker, MaskParams

_ALLOWED = set("ACGTNacgtn-")


def sanitise_alignment(source: Path, dest: Path, logger: logging.Logger) -> Path:
    """Copy ``source`` with every character outside ACGTN- replaced by N.

    Gubbins accepts only those symbols; SNP typers can leave IUPAC ambiguity
    codes (e.g. a heterozygous call rendered by ``bcftools consensus``) in a
    whole-genome alignment. Returns ``dest`` (or ``source`` when unchanged).
    """
    replaced = 0
    with open(source, encoding="utf-8") as fi, open(dest, "w", encoding="utf-8") as fo:
        for line in fi:
            if line.startswith(">"):
                fo.write(line)
                continue
            body = line.rstrip("\n")
            fixed = "".join(c if c in _ALLOWED else "N" for c in body)
            replaced += sum(1 for a, b in zip(body, fixed, strict=True) if a != b)
            fo.write(fixed + "\n")
    if replaced:
        logger.warning(
            "Replaced %d ambiguous base(s) with N before Gubbins (only ACGTN- are accepted).",
            replaced,
        )
        return dest
    dest.unlink(missing_ok=True)
    return source


def read_fasta(path: Path) -> dict[str, str]:
    """Record id -> sequence (single-line or wrapped), in file order."""
    records: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        records[name] = "".join(parts)
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n{seq}\n")


def read_recombination_gff(path: Path) -> dict[str, list[tuple[int, int]]]:
    """Recombinant regions per taxon from Gubbins' GFF3 (1-based, inclusive).

    Each feature line carries ``taxa="  a  b"`` in its attributes; a region
    applies to every taxon listed there.
    """
    regions: dict[str, list[tuple[int, int]]] = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            start, end = int(cols[3]), int(cols[4])
            attrs = cols[8]
            key = 'taxa="'
            i = attrs.find(key)
            if i < 0:
                continue
            taxa = attrs[i + len(key) : attrs.find('"', i + len(key))].split()
            for taxon in taxa:
                regions.setdefault(taxon, []).append((start, end))
    return regions


def apply_masks(
    records: dict[str, str], regions: dict[str, list[tuple[int, int]]]
) -> dict[str, str]:
    """Set the recombinant regions to N in the taxa they were predicted for."""
    out: dict[str, str] = {}
    for name, seq in records.items():
        spans = regions.get(name)
        if not spans:
            out[name] = seq
            continue
        chars = list(seq)
        for start, end in spans:
            lo, hi = max(start - 1, 0), min(end, len(chars))
            chars[lo:hi] = ["N"] * (hi - lo)
        out[name] = "".join(chars)
    return out


def polymorphic_sites(records: dict[str, str]) -> dict[str, str]:
    """Columns with two or more distinct bases among ACGT (N and gaps ignored)."""
    names = list(records)
    seqs = [records[n].upper() for n in names]
    length = min(len(s) for s in seqs) if seqs else 0
    keep = [i for i in range(length) if len({s[i] for s in seqs if s[i] in "ACGT"}) >= 2]
    return {n: "".join(s[i] for i in keep) for n, s in zip(names, seqs, strict=True)}


# Gubbins asks for one of these whenever it runs with more than one thread;
# some conda builds (osx-arm64 among them) ship only the single-threaded
# ``raxmlHPC``, and Gubbins then exits without building a tree.
MULTITHREADED_RAXML = (
    "raxmlHPC-PTHREADS-AVX2",
    "raxmlHPC-PTHREADS-AVX",
    "raxmlHPC-PTHREADS-SSE3",
    "raxmlHPC-PTHREADS",
)
IQTREE_BINARIES = ("iqtree2", "iqtree")


def multithreaded_raxml_available(path: str | None = None) -> bool:
    return any(shutil.which(name, path=path) for name in MULTITHREADED_RAXML)


def resolve_tree_builder(
    requested: str | None,
    threads: int,
    logger: logging.Logger,
    *,
    on_host: bool,
    path: str | None = None,
) -> tuple[str | None, int]:
    """Pick the Gubbins tree builder and thread count that can actually run.

    A requested builder is passed through unchanged. Otherwise Gubbins' own
    default (RAxML) stands unless the run is native, multi-threaded and the
    host has no multi-threaded RAxML build: then IQ-TREE is used when present,
    else Gubbins runs single-threaded. Both fallbacks are logged. ``path`` is
    the PATH Gubbins runs with (its ``--bin-dir`` first), None for PATH.
    """
    if requested:
        return requested, threads
    if threads <= 1 or not on_host or multithreaded_raxml_available(path):
        return None, threads
    if any(shutil.which(name, path=path) for name in IQTREE_BINARIES):
        logger.warning(
            "No multi-threaded RAxML build (%s) on PATH; Gubbins would exit with "
            "--threads %d. Using --tree-builder iqtree instead (set "
            "--tool-arg gubbins_tree_builder=... to choose).",
            "/".join(MULTITHREADED_RAXML),
            threads,
        )
        return "iqtree", threads
    logger.warning(
        "No multi-threaded RAxML build (%s) or IQ-TREE on PATH; running Gubbins "
        "with a single thread.",
        "/".join(MULTITHREADED_RAXML),
    )
    return None, 1


# Above this fraction of variable columns a set is far outside what a
# recombination scanner is built for: Gubbins expects isolates of one species.
DIVERGENCE_WARN_FRACTION = 0.10
# Above this fraction masking is refused unless gubbins_force=true: on the
# Francisella subsets Gubbins finished at 8.3% variable columns and died
# with a bus error at 13.2% (docs/verification.md).
DIVERGENCE_REFUSE_FRACTION = 0.12
# Below this many alignment columns the variable fraction says nothing about
# divergence (a handful of columns is often all variable), so neither the
# warning nor the refusal applies.
MIN_COLUMNS_FOR_DIVERGENCE = 1000
# With IQ-TREE as Gubbins' tree builder, tree building took 76 of 87 minutes
# on a genus-level set (docs/verification.md); a rapidnj first tree makes the
# first iteration cheap, and the resolved builder still builds the later ones.
DEFAULT_FIRST_TREE_BUILDER = ("rapidnj", "JC")


# An NCBI assembly accession (GCA_/GCF_ with a version), a RefSeq/GenBank
# nucleotide accession (NC_000001.1, CP012345.1, MN908947.3) or an ENA/SRA
# run (SRR..., ERR..., DRR...): the forms repgenr itself names genomes by.
_ACCESSION_LIKE = re.compile(r"^(GC[AF]_\d+\.\d+|[A-Z]{1,2}_?\d{5,}(\.\d+)?|[SED]RR\d+)$")


def species_spread(names: list[str]) -> int:
    """Distinct species among record names of the canonical
    ``Family_Genus_species_Accession`` form, counted only where the last
    token is an accession, so a user's own four-token names (such as
    ``Francisella_tularensis_LVS_1``) are not read as taxonomy; 0 when no
    name qualifies."""
    species: set[tuple[str, str]] = set()
    for name in names:
        _family, genus, sp, acc = parse_genome_filename(name)
        if genus and acc and _ACCESSION_LIKE.match(acc):
            species.add((genus.lower(), sp.lower()))
    return len(species)


def check_species_spread(names: list[str], extra: dict, logger: logging.Logger) -> None:
    """Refuse a set that spans more than one species by genome name unless
    ``gubbins_force`` is set, in which case a warning is logged."""
    n_species = species_spread(names)
    if n_species <= 1:
        return
    if _forced(extra):
        logger.warning(
            "Gubbins runs on a set spanning %d species by genome name (gubbins_force=true); "
            "it is built for isolates of one species.",
            n_species,
        )
        return
    raise UserInputError(
        f"The alignment spans {n_species} species by genome name; Gubbins is built "
        "for isolates of one species. Use --mask none, a species-level target, or "
        "--tool-arg gubbins_force=true to run it anyway."
    )


def _forced(extra: dict) -> bool:
    raw = str(extra.get("gubbins_force", "")).strip().lower()
    if raw in ("", "0", "false", "no"):
        return False
    if raw in ("1", "true", "yes"):
        return True
    raise UserInputError(f"gubbins_force must be true or false, got {extra.get('gubbins_force')!r}")


def variable_fraction(records: dict[str, str], sample: int = 200_000) -> tuple[float, int]:
    """Fraction of columns with more than one base among ACGT, and the length.

    Estimated on up to ``sample`` evenly spaced columns so a genus-scale
    alignment (tens of genomes, millions of columns) costs a second or two.
    """
    seqs = [s.upper() for s in records.values()]
    length = min((len(s) for s in seqs), default=0)
    if length == 0 or len(seqs) < 2:
        return 0.0, length
    step = max(1, length // sample)
    columns = range(0, length, step)
    varying = sum(1 for i in columns if len({s[i] for s in seqs if s[i] in "ACGT"}) >= 2)
    return varying / len(list(columns)), length


class GubbinsMasker(Masker):
    capabilities = ToolCapabilities(
        name="gubbins",
        # 3.0: --tree-builder / --first-tree-builder, which the tool-args map onto.
        required_binaries=(
            BinarySpec("run_gubbins.py", version_args=("--version",), min_version="3.0"),
        ),
        container="quay.io/biocontainers/gubbins:3.4.3--py310hfc0ef84_1",
        conda=("bioconda::gubbins",),
        # gubbins_tree_builder / gubbins_first_tree_builder name Gubbins'
        # --tree-builder / --first-tree-builder; gubbins_args is a quoted
        # string of further run_gubbins.py arguments.
        accepted_extras=frozenset(
            {
                "gubbins_tree_builder",
                "gubbins_first_tree_builder",
                "gubbins_args",
                # gubbins_divergence_limit: the variable-column fraction above
                # which masking is refused; gubbins_force=true runs anyway, also
                # on a set spanning more than one species.
                "gubbins_divergence_limit",
                "gubbins_force",
            }
        ),
    )

    _tree: Path | None = None

    def _rapidnj_available(self) -> bool:
        """rapidnj on the host, or assumed in the container that carries Gubbins."""
        if not runs_on_host(self.capabilities):
            return True
        return host_which("rapidnj", self.capabilities.name) is not None

    def precheck(self, names: list[str], extra: dict, logger: logging.Logger) -> None:
        """What can be refused before SNP typing: a set spanning more than one
        species by genome name, and an unusable gubbins_divergence_limit."""
        _divergence_limit(extra)
        _forced(extra)
        check_species_spread(names, extra, logger)

    def tree_output(self) -> Path | None:
        """Gubbins' own tree from the last ``mask`` call (``final_tree.tre``):
        its recombination-aware maximum-likelihood tree, kept for comparison
        with the tree built from the filtered sites."""
        return self._tree if self._tree is not None and self._tree.exists() else None

    def mask(
        self,
        full_alignment: Path,
        out_dir: Path,
        params: MaskParams,
        logger: logging.Logger,
    ) -> Path:
        """Run Gubbins on the whole-genome alignment; return the
        recombination-filtered polymorphic-sites FASTA."""
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = out_dir / "gubbins"
        cleaned = sanitise_alignment(full_alignment, out_dir / "input_alignment.fasta", logger)
        records = read_fasta(cleaned)
        excluded = {name for name in records if name in params.exclude}
        unmatched = sorted(set(params.exclude) - excluded)
        if unmatched:
            # A typer that names records otherwise than core.contracts.record_name
            # (Path.stem gives 'x.fasta' for x.fasta.gz) would put the outgroup
            # into the recombination scan without notice.
            logger.warning(
                "Gubbins: no alignment record is named %s, so it is not left out of "
                "the recombination scan. The SNP typer should name records by the "
                "genome file name without its FASTA suffix and .gz.",
                ", ".join(unmatched),
            )
        if excluded:
            # Gubbins sees the ingroup only; its predictions are applied to
            # every record below, so the excluded outgroup keeps its place.
            subset = out_dir / "ingroup_alignment.fasta"
            write_fasta(subset, {n: s for n, s in records.items() if n not in excluded})
            logger.info(
                "Gubbins runs on %d ingroup genome(s); %s stay(s) out of the scan.",
                len(records) - len(excluded),
                ", ".join(sorted(excluded)),
            )
            gubbins_input = subset
        else:
            gubbins_input = cleaned
        requested = params.extra.get("gubbins_tree_builder")
        # Gubbins hands its thread count to the tree builder, and IQ-TREE
        # refuses more threads than CPU cores; the default is 16.
        threads = params.threads
        cpus = available_cpus(self.capabilities)
        if threads > cpus:
            logger.info(
                "Gubbins: %d threads requested, %d CPU(s) available; using %d.",
                threads,
                cpus,
                cpus,
            )
            threads = cpus
        tree_builder, threads = resolve_tree_builder(
            str(requested) if requested else None,
            threads,
            logger,
            on_host=runs_on_host(self.capabilities),
            path=tool_path(self.capabilities.name),
        )
        scan_records = read_fasta(gubbins_input) if excluded else records
        fraction, length = variable_fraction(scan_records)
        logger.info(
            "Gubbins input: %d genome(s), %d columns, about %.1f%% variable.",
            len(scan_records),
            length,
            100 * fraction,
        )
        force = _forced(params.extra)
        divergence_limit = _divergence_limit(params.extra)
        diverse = length >= MIN_COLUMNS_FOR_DIVERGENCE
        if diverse and fraction > divergence_limit and force:
            logger.warning(
                "About %.0f%% of the alignment is variable, above the %.0f%% at which Gubbins "
                "is otherwise refused; it runs because gubbins_force=true.",
                100 * fraction,
                100 * divergence_limit,
            )
        if diverse and fraction > divergence_limit and not force:
            raise UserInputError(
                f"About {100 * fraction:.0f}% of the alignment is variable, above the "
                f"{100 * divergence_limit:.0f}% at which Gubbins is refused: it is built for "
                "isolates of one species and has died on sets near this divergence. Use "
                "--mask none, "
                "a species-level target, --tool-arg gubbins_divergence_limit=FRACTION, or "
                "--tool-arg gubbins_force=true to run it anyway."
            )
        check_species_spread(list(scan_records), params.extra, logger)
        if diverse and fraction > DIVERGENCE_WARN_FRACTION:
            logger.warning(
                "About %.0f%% of the alignment is variable. Gubbins is built for "
                "isolates of one species, where a few percent is usual; on a set "
                "this diverse it may exhaust its stack or return regions that mean "
                "little. Consider a species-level target or --mask none.",
                100 * fraction,
            )
        argv: list[str | Path] = ["run_gubbins.py", "--threads", str(threads)]
        if tree_builder:
            argv += ["--tree-builder", tree_builder]
        first = params.extra.get("gubbins_first_tree_builder")
        user_args = shlex.split(str(params.extra.get("gubbins_args", "")))
        if first:
            argv += ["--first-tree-builder", str(first)]
        elif (
            "--first-tree-builder" not in user_args
            and not any(
                a.startswith("--first-model") or a == "--custom-first-model" for a in user_args
            )
            and self._rapidnj_available()
        ):
            builder_name, model = DEFAULT_FIRST_TREE_BUILDER
            argv += ["--first-tree-builder", builder_name, "--first-model", model]
            logger.info(
                "Gubbins first tree: %s with the %s model (gubbins_first_tree_builder=... "
                "to change); later iterations use %s.",
                builder_name,
                model,
                tree_builder or "Gubbins' default",
            )
        limit = _filter_percentage(user_args)
        if limit is None:
            # Gubbins leaves taxa with more than 25% gaps or N out of its
            # analysis by default, and still writes them to its outputs, so
            # they are never scanned for recombination. A SNP typer that masks
            # unaligned positions (the simple typer) yields such taxa
            # routinely; keep every taxon in the analysis.
            argv += ["--filter-percentage", "100"]
        else:
            _refuse_filtered_taxa(scan_records, limit)
        argv += user_args
        argv += ["--prefix", prefix, gubbins_input]
        try:
            run_tool(
                self.capabilities,
                argv,
                logger=logger,
                cwd=out_dir,
                log_prefix="gubbins",
            )
        except ToolExecutionError as exc:
            # Keep the tool error, and with it the return code the Nextflow
            # retry rule keys on; add the diagnosis to its output tail.
            diagnosis = (
                f"Gubbins failed on an alignment of {len(scan_records)} genome(s) "
                f"whose columns are about {100 * fraction:.0f}% variable. Its "
                "recombination scan allocates per-SNP arrays on the thread stack "
                "and dies on very diverse input; a within-species set is what it "
                "expects. Run the stage with --mask none, or narrow the target."
            )
            raise ToolExecutionError(
                exc.command, exc.returncode, f"{exc.output or ''}\n{diagnosis}".strip()
            ) from exc
        self._tree = Path(str(prefix) + ".final_tree.tre")
        if not excluded:
            filtered = Path(str(prefix) + ".filtered_polymorphic_sites.fasta")
            if not filtered.exists():
                raise WorkdirError("Gubbins did not produce a filtered polymorphic sites FASTA")
            return filtered
        gff = Path(str(prefix) + ".recombination_predictions.gff")
        if not gff.exists():
            raise WorkdirError("Gubbins did not produce a recombination predictions GFF")
        regions = read_recombination_gff(gff)
        masked = apply_masks(records, regions)
        n_sites = sum(len(r) for r in regions.values())
        logger.info(
            "Masked %d recombinant region(s) across %d genome(s); extracting variable sites.",
            n_sites,
            len(regions),
        )
        out = Path(str(prefix) + ".masked_polymorphic_sites.fasta")
        write_fasta(out, polymorphic_sites(masked))
        return out


def _divergence_limit(extra: dict) -> float:
    """The variable-column fraction above which masking is refused."""
    raw = extra.get("gubbins_divergence_limit")
    if raw is None or str(raw).strip() == "":
        return DIVERGENCE_REFUSE_FRACTION
    try:
        value = float(str(raw))
    except ValueError as exc:
        raise UserInputError(
            f"gubbins_divergence_limit needs a fraction between 0 and 1, got {raw!r}"
        ) from exc
    if not 0 < value <= 1:
        raise UserInputError(
            f"gubbins_divergence_limit needs a fraction between 0 and 1, got {raw!r}"
        )
    return value


def _filter_percentage(args: list[str]) -> float | None:
    """The ``--filter-percentage``/``-f`` value in user arguments, or None."""
    for i, arg in enumerate(args):
        name, eq, value = arg.partition("=")
        if name in ("--filter-percentage", "-f"):
            raw = value if eq else (args[i + 1] if i + 1 < len(args) else "")
            try:
                return float(raw)
            except ValueError:
                raise UserInputError(
                    f"--filter-percentage in gubbins_args needs a number, got {raw!r}"
                ) from None
    return None


def _refuse_filtered_taxa(records: dict[str, str], limit: float) -> None:
    """Refuse taxa that Gubbins' ``--filter-percentage`` would leave out of its analysis.

    Gubbins counts N, n and '-' as missing data, excludes a taxon above the
    limit from the recombination scan, and still writes it to its outputs, so
    the omission is not visible in the result. Checked here, before the run.
    """
    over = sorted(
        name
        for name, seq in records.items()
        if seq and 100 * (seq.count("N") + seq.count("n") + seq.count("-")) / len(seq) > limit
    )
    if over:
        raise WorkdirError(
            f"{len(over)} taxon/taxa have more than {limit:g}% gaps or N, so Gubbins "
            f"(--filter-percentage {limit:g}) would leave them out of the recombination "
            f"analysis while keeping them in its output: {', '.join(over[:10])}"
            f"{' ...' if len(over) > 10 else ''}. Remove --filter-percentage from "
            "--tool-arg gubbins_args (repgenr then passes 100), leave those genomes out, "
            "or choose a reference closer to them."
        )
