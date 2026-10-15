"""Shared plugin infrastructure: capability metadata + entry-point registry.

Each tool family (dereplicators, aligners, snptypers, treebuilders) defines an
ABC and instantiates a :class:`Registry` bound to its entry-point group. In-tree
adapters and third-party packages are discovered identically through
``importlib.metadata`` entry points, so core never imports a concrete adapter
and external tools need no core edits.
"""

from __future__ import annotations

import logging
import shutil
from collections.abc import Collection, Mapping
from dataclasses import dataclass, field, fields, is_dataclass
from importlib.metadata import entry_points

from .binaries import BinarySpec, check_binaries
from .errors import PluginError, UserInputError


def parse_extra_int(extra: Mapping[str, object], key: str, default: int) -> int:
    """Read an integer tool override from ``--extra`` with a clear error.

    A non-integer value (e.g. ``--extra ksize=big``) raises
    :class:`UserInputError` naming the key instead of a raw ``ValueError``
    traceback from ``int()``.
    """
    raw = extra.get(key, default)
    try:
        return int(raw)  # type: ignore[call-overload]
    except (TypeError, ValueError) as exc:
        raise UserInputError(f"Extra parameter '{key}' must be an integer, got {raw!r}.") from exc


def parse_extra_bool(extra: Mapping[str, object], key: str, default: bool) -> bool:
    """Read a boolean tool override (``true``/``false``, ``1``/``0``, ``yes``/``no``);
    an absent key or an empty value gives ``default``."""
    raw = extra.get(key, default)
    if isinstance(raw, bool):
        return raw
    text = str(raw).strip().lower()
    if text == "":
        return default
    if text in ("1", "true", "yes"):
        return True
    if text in ("0", "false", "no"):
        return False
    raise UserInputError(f"Extra parameter '{key}' must be true or false, got {raw!r}.")


@dataclass(frozen=True)
class ToolCapabilities:
    """Declarative metadata for one tool adapter.

    ``recommended_max_genomes`` drives auto-selection and scale warnings.
    ``supports_native_scaling`` marks a dereplicator that scales to large sets in
    one pass (e.g. skDER), so it is single-pass by default; it no longer gates
    chunking -- ``--process-size`` opts any tool into the two-stage chunked path.
    """

    name: str
    required_binaries: tuple[BinarySpec, ...] = ()
    default_params: dict = field(default_factory=dict)
    recommended_max_genomes: int | None = None
    supports_native_scaling: bool = False
    # Container execution: a pinned image URI (BioContainer or Wave-minted) for
    # single-tool adapters, or a conda spec resolved to an image via Wave when no
    # explicit image is set and Wave is enabled.
    container: str | None = None
    conda: tuple[str, ...] = ()
    # Extra-dict keys this adapter actually reads; used to warn on (and avoid
    # injecting) tuning that a tool would silently ignore.
    accepted_extras: frozenset[str] = frozenset()
    # Names of the family's standard parameters (``primary_ani``, ``threads``,
    # ``bootstrap``...) the adapter never passes to its tool. The stage warns
    # when a user sets one of them to a non-default value.
    ignored_params: frozenset[str] = frozenset()
    # The tool reads gzipped FASTA genomes itself. When False (the default, the
    # safe side for a new adapter), the stages hand an aligner, SNP typer or
    # genome-input tree builder decompressed copies (core.process.stage_plain_inputs).
    reads_gzip: bool = False


def warn_unconsumed_extras(
    caps: ToolCapabilities,
    extra: Mapping[str, object],
    logger: logging.Logger,
    *,
    family: str,
) -> None:
    """Warn once, by name, about extra keys the adapter never reads.

    Keys that are also in ``default_params`` are the adapter's own defaults and
    are not reported.
    """
    unread = sorted(set(extra) - caps.accepted_extras - set(caps.default_params))
    if unread:
        logger.warning(
            "%s '%s' ignores extra parameter(s): %s", family, caps.name, ", ".join(unread)
        )


def warn_ignored_params(
    caps: ToolCapabilities, params: object, logger: logging.Logger, *, family: str
) -> None:
    """Warn, by name, about standard parameters the adapter declares it ignores.

    ``params`` is the family's parameter dataclass. Only values that differ
    from the field default are reported: a default the tool happens not to
    read is not a user request that went missing.
    """
    if not caps.ignored_params or not is_dataclass(params):
        return
    defaults = {f.name: f.default for f in fields(params)}
    for name in sorted(caps.ignored_params):
        if name not in defaults:
            continue
        value = getattr(params, name)
        if value == defaults[name]:
            continue
        logger.warning(
            "%s '%s' does not use %s; the value %r is ignored.", family, caps.name, name, value
        )


class Registry[T]:
    """Lazily-loaded registry of adapter classes for one entry-point group."""

    def __init__(self, group: str):
        self.group = group
        self._classes: dict[str, type[T]] = {}
        self._loaded = False
        # Broken plugins auto-selection has already warned about: a stage
        # selects twice (precheck, then run) and should warn once.
        self._skip_warned: set[str] = set()

    def _load(self) -> None:
        if self._loaded:
            return
        for ep in entry_points(group=self.group):
            try:
                self._classes[ep.name] = ep.load()
            except Exception as exc:  # a broken third-party plugin must not kill the run
                # Warned here and raised only if the broken name is requested, so
                # one broken plugin does not disable the rest of the family.
                logging.getLogger("repgenr").warning(
                    "Plugin %r (group %s) failed to load and is unavailable: %s",
                    ep.name,
                    self.group,
                    exc,
                )
                self._classes.setdefault(ep.name, _BrokenPlugin(ep.name, exc))  # type: ignore[arg-type]
        self._loaded = True

    def names(self) -> list[str]:
        self._load()
        return sorted(self._classes)

    def register(self, name: str, cls: type[T], *, replace: bool = False) -> None:
        """Register an adapter class under ``name`` (tests, embedding callers).

        Entry points remain the discovery path for installed packages; this is
        the public API for programmatic registration, replacing direct mutation
        of the private class table.
        """
        self._load()
        if name in self._classes and not replace:
            raise PluginError(
                f"Tool '{name}' is already registered for {self.group} "
                "(pass replace=True to override)."
            )
        self._classes[name] = cls

    def unregister(self, name: str) -> None:
        self._load()
        self._classes.pop(name, None)

    def is_broken(self, name: str) -> bool:
        """True when ``name`` is registered but its adapter failed to import."""
        self._load()
        return isinstance(self._classes.get(name), _BrokenPlugin)

    def load_error(self, name: str) -> Exception | None:
        """The import error of a broken plugin, or None."""
        self._load()
        cls = self._classes.get(name)
        return cls.error if isinstance(cls, _BrokenPlugin) else None

    def get(self, name: str) -> type[T]:
        self._load()
        if name not in self._classes:
            available = ", ".join(self.names()) or "none"
            raise PluginError(f"Unknown tool '{name}' for {self.group}. Available: {available}")
        cls = self._classes[name]
        if isinstance(cls, _BrokenPlugin):
            raise PluginError(f"Plugin '{name}' failed to load: {cls.error}") from cls.error
        return cls

    def create(self, name: str) -> T:
        return self.get(name)()  # type: ignore[call-arg]


class _BrokenPlugin:
    """Placeholder for an entry point that failed to import."""

    def __init__(self, name: str, error: Exception):
        self.name = name
        self.error = error


def preflight(capabilities: ToolCapabilities) -> dict[str, str]:
    """Check the adapter's required binaries; return resolved versions.

    When a container backend is active and an image resolves for this tool, the
    tool lives in the image (not on the host): check the engine binary, and
    that the Docker daemon answers, instead, and record the image reference in
    place of host tool versions, with the engine's name and version beside it
    (``{"iqtree": "<image>", "docker": "29.5.3"}``).
    """
    from pathlib import Path

    from .containers import check_engine_ready, get_config, resolve_image  # avoids a cycle

    config = get_config()
    if config.active:
        image = resolve_image(capabilities, config)
        if image:
            engine = config.engine_binary()
            found = check_binaries((BinarySpec(engine, version_args=("--version",)),))
            check_engine_ready(config)
            return {capabilities.name: image, Path(engine).name: found[engine]}
    from .bindirs import note_used, tool_path

    note_used(capabilities.name)
    return check_binaries(capabilities.required_binaries, path=tool_path(capabilities.name))


AUTO = "auto"


def tool_choices_help(registry: Registry, *, auto: bool = True, prefix: str = "") -> str:
    """Help text naming the currently registered tools for a CLI option.

    Generated from the registry so third-party adapters appear in --help and
    the text can never drift from the actual choices.
    """
    names = registry.names()
    choices = (["auto"] if auto else []) + names
    # Comma-space separators so Rich can wrap the list in a narrow help
    # panel instead of ellipsizing one long slash-joined token.
    return f"{prefix}{', '.join(choices) or '(none registered)'}."


def _capabilities_of(registry: Registry, name: str) -> ToolCapabilities | None:
    try:
        return registry.get(name).capabilities  # type: ignore[attr-defined]
    except PluginError:
        return None


def _tool_available(cap: ToolCapabilities) -> bool:
    """Is this adapter runnable in the CURRENT execution environment?

    This mirrors :func:`repgenr.core.containers.run_tool`. Under an active
    container backend a tool runs in an image when it pins one, or when
    ``--wave`` is on and it declares a conda spec; such a tool counts as
    available (resolution itself happens at preflight). Any other tool runs on
    the host, so it is available only when its required binaries are on PATH.
    """
    from .bindirs import tool_path
    from .containers import get_config

    config = get_config()
    if config.active and (cap.container is not None or (config.wave_enabled and bool(cap.conda))):
        return True
    path = tool_path(cap.name)
    return all(shutil.which(spec.name, path=path) is not None for spec in cap.required_binaries)


def tool_available(caps: ToolCapabilities) -> bool:
    """Whether an adapter can run here: on the host, or in an image under a backend."""
    return _tool_available(caps)


# Tie-break order for auto-selection, matching the documented per-family
# defaults; unlisted tools rank after these, alphabetically.
_PREFERRED_ORDER = (
    "skder",
    "iqtree",
    "progressivemauve",
    "simple",
    "skesa",
    "shovill",
    "flye",
    "medaka",
    "racon",
)


def _preference_rank(name: str) -> int:
    try:
        return _PREFERRED_ORDER.index(name)
    except ValueError:
        return len(_PREFERRED_ORDER)


def auto_select(
    registry: Registry, n_items: int, allowed: Collection[str] | None = None
) -> str | None:
    """Pick the best *available* registered tool for ``n_items`` inputs.

    ``allowed`` restricts the candidates to the named tools.

    Preference order: runnable in the current environment (container-aware),
    then fits the recommended scale, then the TIGHTEST fitting limit --
    higher-quality tools declare tighter recommended scales, so the loosest
    tool must not win small inputs -- with unbounded tools last among the
    fitting ones. Ties follow the documented defaults, then alphabetical.
    When nothing fits, the largest-capacity tool is chosen.
    """
    inf = float("inf")
    best: tuple[tuple, str] | None = None
    for name in registry.names():
        if allowed is not None and name not in allowed:
            continue
        cap = _capabilities_of(registry, name)
        if cap is None:
            if name not in registry._skip_warned:
                registry._skip_warned.add(name)
                logging.getLogger("repgenr").warning(
                    "auto-select skipping '%s' (%s): the plugin failed to load.",
                    name,
                    registry.group,
                )
            continue
        limit = cap.recommended_max_genomes
        limit_value = inf if limit is None else float(limit)
        fits = limit is None or limit >= n_items
        # Among fitting tools smaller limits sort first (tightest fit); when
        # nothing fits, larger capacity sorts first.
        key = (
            0 if _tool_available(cap) else 1,
            0 if fits else 1,
            limit_value if fits else -limit_value,
            _preference_rank(name),
            name,
        )
        if best is None or key < best[0]:
            best = (key, name)
    return best[1] if best else None


def scale_warning(registry: Registry, tool: str, n_items: int) -> tuple[int, list[str]] | None:
    """If ``tool`` is over its recommended scale, return (limit, alternatives).

    Alternatives are registered tools whose recommended scale accommodates
    ``n_items``. Returns None when the tool is within its recommended scale.
    """
    cap = _capabilities_of(registry, tool)
    if cap is None or cap.recommended_max_genomes is None:
        return None
    if n_items <= cap.recommended_max_genomes:
        return None
    alternatives = []
    for name in registry.names():
        if name == tool:
            continue
        other = _capabilities_of(registry, name)
        if other is None:
            continue
        if other.recommended_max_genomes is None or other.recommended_max_genomes >= n_items:
            alternatives.append(name)
    return cap.recommended_max_genomes, sorted(alternatives)
