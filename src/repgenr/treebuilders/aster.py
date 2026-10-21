"""ASTER tree builder: a species tree summarised from per-gene trees.

The core-gene route concatenates per-gene alignments and records where each
gene lies (``partitions.json``). With ``--partitioned`` this builder cuts the
alignment back into its genes, infers one tree per gene (IQ-TREE or FastTree)
and summarises the gene trees into one species tree with ASTER
(github.com/chaoszhang/ASTER), the maintained successor of ASTRAL. Summary
methods account for gene-tree discordance from incomplete lineage sorting,
which concatenation does not; they do not apply to SNP alignments, which have
no genes.

Programs (``program=``):

* ``wastral`` (default): weighted ASTRAL, which weights each gene-tree quartet
  by its branch support and length. IQ-TREE gene trees carry aBayes supports
  (``-B``), FastTree gene trees SH-like local supports (``-L``).
* ``astral4``: unweighted ASTRAL-IV; branch lengths in substitutions per site
  (estimated by its integrated CASTLES-2).
* ``astral-pro3``: ASTRAL-Pro3, for multi-copy gene families; here every
  record is one genome, so no gene-to-genome mapping is needed.

Node labels of the species tree are local posterior probabilities in [0, 1].
wastral writes internal branch lengths in coalescent units and no terminal
branch lengths, so its tree cannot be rooted at the midpoint.
"""

from __future__ import annotations

import logging
import re
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from ..core.binaries import BinarySpec, check_binaries
from ..core.containers import get_config, run_chain, run_tool, runs_on_host
from ..core.contracts import atomic_path
from ..core.errors import UserInputError, WorkdirError
from ..core.executors import parallel_map
from ..core.plugins import ToolCapabilities, parse_extra_int, preflight, tool_available
from .base import InputKind, TreeBuilder, TreeParams, as_msa_path, partition_names
from .fasttree import FasttreeBuilder
from .iqtree import IqtreeBuilder

PROGRAMS = ("wastral", "astral4", "astral-pro3")
GENE_TREE_BUILDERS = ("iqtree", "fasttree")
GENE_TREES = "gene_trees.nwk"
# ASTER needs four taxa for a quartet; a gene with fewer informs nothing.
MIN_RECORDS = 4

_IQTREE = IqtreeBuilder.capabilities
_FASTTREE = FasttreeBuilder.capabilities
# Characters that count as a called site in a gene slice.
_CALLED_DNA = frozenset("ACGTUacgtu")
_UNCALLED_PROTEIN = frozenset("Xx-?.*")
# IQ-TREE 3 writes an aBayes-only support as "/0.95" after the closing
# parenthesis, which ASTER cannot parse; the last field is the support.
_SLASH_LABEL = re.compile(r"\)(?:[^:,();]*/)+([^:,();/]*)")


@dataclass(frozen=True)
class _Gene:
    name: str  # file-safe partition name
    fasta: Path
    records: int


class AsterBuilder(TreeBuilder):
    capabilities = ToolCapabilities(
        name="aster",
        container="quay.io/biocontainers/aster:1.25--hc52dbad_1",
        conda=("bioconda::aster",),
        # program=: wastral, astral4 or astral-pro3. gene_tree_builder=:
        # iqtree (-fast, aBayes supports) or fasttree. model=: the IQ-TREE
        # model of every gene tree (GTR+G, or LG+G for protein, by default).
        # rounds=: ASTER's rounds of placement and subsampling (-r, -s).
        # min_gene_length=: genes shorter than this many columns are skipped.
        accepted_extras=frozenset(
            {"program", "gene_tree_builder", "model", "rounds", "min_gene_length"}
        ),
        default_params={
            "program": "wastral",
            "gene_tree_builder": "iqtree",
            "min_gene_length": 0,
        },
        # ASTER prints "Version: v1.25.3.8" in its help. 1.19: the wastral
        # support presets (-B, -L) and the astral4 and astral-pro3 binaries.
        required_binaries=(BinarySpec("wastral", version_args=("-h",), min_version="1.19"),),
        # Declared, not measured.
        recommended_max_genomes=5000,
        # Supports come from the gene trees, not from a bootstrap of the tree.
        ignored_params=frozenset({"bootstrap"}),
        # Only meaningful for the gene partitions of a core-gene alignment.
        auto_eligible=False,
    )
    input_kind = InputKind.MSA_FASTA
    # Local posterior probabilities are always written.
    supports_without_bootstrap = True
    supports_partitions = True

    def writes_terminal_lengths(self, extra: Mapping[str, object]) -> bool:
        """wastral leaves the terminal branches without a length."""
        return str(extra.get("program", "wastral")).strip() != "wastral"

    def preflight(self) -> dict[str, str]:
        """ASTER, and the gene-tree builders that can run here.

        Which gene-tree builder is used is an extra, which preflight does not
        see; build() checks the one in use before it starts.
        """
        versions = preflight(self.capabilities)
        found = [caps for caps in (_IQTREE, _FASTTREE) if tool_available(caps)]
        if not found:
            # Neither can run: report the default's missing binary.
            found = [_IQTREE]
        for caps in found:
            versions.update(preflight(caps))
        return versions

    def build(
        self,
        msa_or_genomes: Path | Sequence[Path],
        out_dir: Path,
        params: TreeParams,
        logger: logging.Logger,
    ) -> Path:
        msa = as_msa_path(msa_or_genomes)
        program = str(params.extra.get("program", "wastral")).strip()
        if program not in PROGRAMS:
            raise UserInputError(
                f"aster: program={program!r} is not supported; use one of {', '.join(PROGRAMS)}."
            )
        gene_builder = str(params.extra.get("gene_tree_builder", "iqtree")).strip().lower()
        if gene_builder not in GENE_TREE_BUILDERS:
            raise UserInputError(
                f"aster: gene_tree_builder={gene_builder!r} is not supported; use one of "
                f"{', '.join(GENE_TREE_BUILDERS)}."
            )
        min_length = parse_extra_int(params.extra, "min_gene_length", 0)
        rounds = _rounds(params)
        model = str(params.extra.get("model", "")).strip()
        if model and gene_builder == "fasttree":
            logger.warning(
                "aster: model=%s applies to gene_tree_builder=iqtree only; FastTree "
                "fits GTR (nucleotide) or LG (protein)",
                model,
            )
        if params.partitions is None:
            raise UserInputError(
                "aster needs the gene partitions of a core-gene alignment; use "
                "--msa-source coregene --partitioned"
            )
        gene_caps = _IQTREE if gene_builder == "iqtree" else _FASTTREE
        # Fail on a missing gene-tree builder or ASTER program before any work.
        preflight(gene_caps)
        version = _program_version(self.capabilities, program)

        out_dir.mkdir(parents=True, exist_ok=True)
        protein = params.alphabet == "protein"
        genes, skipped_short, skipped_few, n_records = _write_genes(
            msa, params.partitions, min_length, out_dir / "genes", protein=protein
        )
        logger.info(
            "aster: %d of %d gene partitions used; %d skipped as shorter than %d columns, "
            "%d with fewer than %d genomes holding a called character",
            len(genes),
            len(params.partitions),
            skipped_short,
            min_length,
            skipped_few,
            MIN_RECORDS,
        )
        if not genes:
            raise UserInputError(
                f"aster: no gene partition has {MIN_RECORDS} or more genomes with a called "
                "character (and passes min_gene_length); a species tree needs gene trees "
                f"of at least {MIN_RECORDS} genomes."
            )

        tree_dir = out_dir / "gene_trees"
        tree_dir.mkdir(parents=True, exist_ok=True)
        if gene_builder == "iqtree":
            steps = [_iqtree_step(g, tree_dir, model, protein) for g in genes]
            tool_env = None
        else:
            steps = [_fasttree_step(g, tree_dir, protein) for g in genes]
            # One thread per gene tree; FastTreeMP reads it from OpenMP.
            tool_env = {"OMP_NUM_THREADS": "1"}
        workers = max(1, params.threads)
        n_batches = min(len(steps), workers)
        batches = [steps[i::n_batches] for i in range(n_batches)]

        def run(batch: list[tuple[str, list[str | Path]]]) -> None:
            run_chain(gene_caps, batch, logger=logger, cwd=out_dir, tool_env=tool_env)

        logger.info(
            "aster: inferring %d gene trees with %s in %d batches, one thread each",
            len(genes),
            gene_builder,
            n_batches,
        )
        # The first failure cancels the batches not yet started.
        parallel_map(run, batches, workers, logger=logger)

        gene_trees = out_dir / GENE_TREES
        taxa = _collect_gene_trees(genes, tree_dir, gene_builder, gene_trees)
        logger.info(
            "aster: %d gene trees in %s over %d of %d genomes",
            len(genes),
            GENE_TREES,
            taxa,
            n_records,
        )

        species = out_dir / "aster.nwk"
        cmd: list[str | Path] = [program]
        if program == "wastral":
            # The support scale of the gene trees: aBayes (IQ-TREE) or SH-like
            # local supports (FastTree), both in [0, 1].
            cmd.append("-B" if gene_builder == "iqtree" else "-L")
        if rounds is not None:
            cmd += ["-r", str(rounds), "-s", str(rounds)]
        cmd += ["-i", gene_trees, "-o", species, "-t", str(max(1, params.threads))]
        units = (
            "coalescent units (internal branches; no terminal lengths)"
            if program == "wastral"
            else "substitutions per site"
        )
        logger.info(
            "aster: species tree with %s%s; branch lengths in %s, node labels local "
            "posterior probabilities",
            program,
            f" {version}" if version else "",
            units,
        )
        run_tool(self.capabilities, cmd, logger=logger, cwd=out_dir, log_prefix=program)
        if not species.is_file() or species.stat().st_size == 0:
            raise WorkdirError(f"{program} wrote no species tree ({species})")
        tree = out_dir / "tree.nwk"
        with atomic_path(tree) as tmp:
            shutil.copy2(species, tmp)
        return tree


def _rounds(params: TreeParams) -> int | None:
    if str(params.extra.get("rounds", "")).strip() == "":
        return None
    rounds = parse_extra_int(params.extra, "rounds", 0)
    if rounds < 1:
        raise UserInputError(f"aster: rounds must be 1 or more, got {rounds}.")
    return rounds


def _read_alignment(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    current: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().split(maxsplit=1)[0] if line[1:].strip() else ""
                if name in seqs:
                    raise WorkdirError(f"The alignment {path} names record {name!r} twice")
                current = seqs.setdefault(name, [])
            else:
                current.append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def _called(segment: str, protein: bool) -> bool:
    if protein:
        return any(c.isalpha() and c not in _UNCALLED_PROTEIN for c in segment)
    return any(c in _CALLED_DNA for c in segment)


def _write_genes(
    msa: Path,
    partitions: Sequence[tuple[str, int, int]],
    min_length: int,
    gene_dir: Path,
    *,
    protein: bool = False,
) -> tuple[list[_Gene], int, int, int]:
    """Cut ``msa`` into one FASTA per partition under ``gene_dir``.

    A record with no called character in a gene is left out of that gene; a
    gene shorter than ``min_length`` columns, or with fewer than MIN_RECORDS
    records left, is skipped. Returns the genes written, the counts skipped
    as short and as too few, and the number of records in the alignment.
    """
    seqs = _read_alignment(msa)
    width = max((len(s) for s in seqs.values()), default=0)
    gene_dir.mkdir(parents=True, exist_ok=True)
    names = _casefold_unique(partition_names([gene for gene, _, _ in partitions]))
    genes: list[_Gene] = []
    short = few = 0
    for name, (gene, start, end) in zip(names, partitions, strict=True):
        if start < 1 or end < start or end > width:
            raise WorkdirError(
                f"Partition {gene!r} ({start}-{end}) lies outside the {width} columns of {msa}"
            )
        if end - start + 1 < min_length:
            short += 1
            continue
        kept = [
            (record, seq[start - 1 : end])
            for record, seq in seqs.items()
            if _called(seq[start - 1 : end], protein)
        ]
        if len(kept) < MIN_RECORDS:
            few += 1
            continue
        path = gene_dir / f"{name}.fasta"
        with open(path, "w", encoding="utf-8") as fo:
            for record, segment in kept:
                fo.write(f">{record}\n{segment}\n")
        genes.append(_Gene(name, path, len(kept)))
    return genes, short, few, len(seqs)


def _casefold_unique(names: Sequence[str]) -> list[str]:
    """``names`` made unique without regard to case.

    The names become file names, and "dnaA" and "DnaA" would be one file on
    a case-insensitive file system; a later one gains a numeric suffix.
    """
    taken: set[str] = set()
    out: list[str] = []
    for base in names:
        name, n = base, 1
        while name.casefold() in taken:
            n += 1
            name = f"{base}_{n}"
        taken.add(name.casefold())
        out.append(name)
    return out


def _iqtree_step(
    gene: _Gene, tree_dir: Path, model: str, protein: bool
) -> tuple[str, list[str | Path]]:
    default_model = "LG+G" if protein else "GTR+G"
    return (
        "iqtree",
        [
            "iqtree",
            "-s",
            gene.fasta,
            "-st",
            "AA" if protein else "DNA",
            "-m",
            model or default_model,
            "--fast",
            "--abayes",
            "-T",
            "1",
            "-pre",
            tree_dir / gene.name,
            "-redo",
            "--quiet",
        ],
    )


def _fasttree_step(gene: _Gene, tree_dir: Path, protein: bool) -> tuple[str, list[str | Path]]:
    from ..core.bindirs import host_which

    # The image's canonical name is lowercase; on the host either spelling.
    if get_config().active or not host_which("FastTree", _FASTTREE.name):
        binary = "fasttree"
    else:
        binary = "FastTree"
    model = ["-lg"] if protein else ["-nt", "-gtr"]
    return (
        "fasttree",
        [binary, *model, "-quiet", "-out", tree_dir / f"{gene.name}.nwk", gene.fasta],
    )


def _collect_gene_trees(genes: Sequence[_Gene], tree_dir: Path, builder: str, out: Path) -> int:
    """Write one Newick per gene to ``out``; return the number of genomes in any."""
    taxa: set[str] = set()
    leaf = re.compile(r"[(,]\s*([^():,;\s]+)")
    with atomic_path(out) as tmp, open(tmp, "w", encoding="utf-8") as fo:
        for gene in genes:
            suffix = ".treefile" if builder == "iqtree" else ".nwk"
            path = tree_dir / f"{gene.name}{suffix}"
            if not path.is_file() or path.stat().st_size == 0:
                raise WorkdirError(f"{builder} wrote no tree for gene {gene.name} ({path})")
            newick = " ".join(path.read_text(encoding="utf-8").split())
            newick = _SLASH_LABEL.sub(r")\1", newick)
            taxa.update(leaf.findall(newick))
            fo.write(newick + "\n")
    return len(taxa)


def _program_version(caps: ToolCapabilities, program: str) -> str:
    """Check the ASTER program on the host and return its version.

    Under a container backend the image holds all three programs and the
    version is not read ('').
    """
    if not runs_on_host(caps):
        return ""
    from ..core.bindirs import tool_path

    found = check_binaries((BinarySpec(program, version_args=("-h",)),), path=tool_path(caps.name))
    return found.get(program, "")
