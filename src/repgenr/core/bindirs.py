"""Per-tool host directories: ``--bin-dir TOOL=DIR`` and ``REPGENR_BIN_DIRS``.

A tool from a satellite conda environment (envs/gubbins.yml and others) can be
given its environment's ``bin`` directory. For that tool only, the directory is
put first on PATH: when its binaries are looked up at preflight, when its
version is queried, and in the environment of every host subprocess it runs,
so a wrapper script calls the helpers of its own environment. Other tools keep
the inherited PATH, so a satellite never shadows the core environment.

The setting is process-global, like the container backend. With no directory
configured every lookup is unchanged. Under a container backend a tool that
runs in an image is not affected. The directories are recorded per tool in
the stage record and kept out of the resume fingerprint: they choose where a
tool is found, and its version is already recorded.
"""

from __future__ import annotations

import os
import shutil
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path

from .errors import UserInputError
from .plugins import ToolCapabilities

ENV_VAR = "REPGENR_BIN_DIRS"


@dataclass
class HostConfig:
    bin_dirs: dict[str, Path] = field(default_factory=dict)


_HOST = HostConfig()
# Tools whose directory was used since the last reset (the stage harness
# resets it before a stage body runs and records it afterwards).
_USED: dict[str, str] = {}


def configure_bin_dirs(mapping: Mapping[str, str | os.PathLike[str]]) -> HostConfig:
    """Set the process-global per-tool directories (an empty mapping clears them)."""
    global _HOST
    _HOST = HostConfig({name: Path(path) for name, path in mapping.items()})
    _USED.clear()
    return _HOST


def get_bin_dirs() -> dict[str, Path]:
    return dict(_HOST.bin_dirs)


def tool_path(name: str, base: str | None = None) -> str | None:
    """PATH for tool ``name``: its directory, then ``base`` (default: PATH).

    None when no directory is configured for the tool, which means "use the
    inherited PATH unchanged".
    """
    directory = _HOST.bin_dirs.get(name)
    if directory is None:
        return None
    rest = base if base is not None else os.environ.get("PATH", os.defpath)
    return os.pathsep.join([str(directory), rest]) if rest else str(directory)


def note_used(name: str) -> None:
    directory = _HOST.bin_dirs.get(name)
    if directory is not None:
        _USED[name] = str(directory)


def used() -> dict[str, str]:
    return dict(_USED)


def reset_used() -> None:
    _USED.clear()


def host_which(binary: str, tool: str) -> str | None:
    """``shutil.which`` for a binary that belongs to adapter ``tool``."""
    return shutil.which(binary, path=tool_path(tool))


def parse_entries(entries: Iterable[str], source: str = "--bin-dir") -> dict[str, str]:
    """``["gubbins=/p/bin", ...]`` -> ``{"gubbins": "/p/bin"}``; blank entries are skipped."""
    out: dict[str, str] = {}
    for raw in entries:
        entry = raw.strip()
        if not entry:
            continue
        name, sep, directory = entry.partition("=")
        name, directory = name.strip(), directory.strip()
        if not sep or not name or not directory:
            raise UserInputError(
                f"Bad {source} entry {entry!r}: expected TOOL=DIR, e.g. gubbins=/path/env/bin."
            )
        out[name] = directory
    return out


def known_capabilities() -> dict[str, list[ToolCapabilities]]:
    """Tool name -> capabilities, for every adapter and the stage-level tools."""
    import importlib

    from ..polishers.racon import _MINIMAP2
    from ..stages.assemble_qc import CHECKM2_CAPS
    from ..stages.genome import DATASETS_CAPS

    caps: dict[str, list[ToolCapabilities]] = {}
    for family in (
        "aligners",
        "annotators",
        "assemblers",
        "classifiers",
        "coregenes",
        "dereplicators",
        "maskers",
        "polishers",
        "snptypers",
        "treebuilders",
    ):
        registry = importlib.import_module(f"repgenr.{family}.base").registry
        for name in registry.names():
            if registry.is_broken(name):
                caps.setdefault(name, [])
                continue
            caps.setdefault(name, []).append(registry.get(name).capabilities)
    from ..tree.nj import RAPIDNJ  # lazy: tree.nj imports core

    for extra in (CHECKM2_CAPS, DATASETS_CAPS, _MINIMAP2, RAPIDNJ):
        caps.setdefault(extra.name, []).append(extra)
    return caps


def validate(mapping: Mapping[str, str], source: str = "--bin-dir") -> dict[str, Path]:
    """Check tool names and directories; return absolute directories.

    Raises :class:`UserInputError` (exit 2 at the CLI) for an unknown tool or
    a directory that does not exist.
    """
    if not mapping:
        return {}
    known = known_capabilities()
    unknown = sorted(set(mapping) - set(known))
    if unknown:
        raise UserInputError(
            f"{source} names unknown tool(s): {', '.join(unknown)}. "
            f"Known tools: {', '.join(sorted(known))}."
        )
    resolved: dict[str, Path] = {}
    for name, directory in mapping.items():
        path = Path(directory).expanduser()
        if not path.is_dir():
            raise UserInputError(f"{source} {name}={directory}: not a directory.")
        resolved[name] = path.resolve()
    return resolved


def ineffective_under_backend(*, wave_enabled: bool) -> list[str]:
    """Configured tools that run in an image under an active container backend."""
    known = known_capabilities()
    return sorted(
        name
        for name in _HOST.bin_dirs
        if any(
            cap.container is not None or (wave_enabled and bool(cap.conda))
            for cap in known.get(name, [])
        )
    )
