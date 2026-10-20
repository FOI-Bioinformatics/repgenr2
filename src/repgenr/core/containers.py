"""Container execution backend for external tools.

RepGenR2 drives bioinformatics tools via :func:`repgenr.core.process.run`. When a
container backend is active, each tool's argument vector is wrapped so it runs
inside a pinned image (a BioContainer, the Cactus image, or a Wave-built image),
which unblocks Linux-only tools on any host and pins tool versions.

Only the external-tool subprocess is containerized; RepGenR2's own Python (e.g.
the format converters) keeps running on the host and reads tool outputs from the
bind-mounted working directory.

The backend is process-global: configure it once (CLI/env), then adapters call
:func:`run_tool` with their capabilities.
"""

from __future__ import annotations

import logging
import os
import re
import shlex
import shutil
import socket
import subprocess
import tempfile
import time
import uuid
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from . import process
from .errors import MissingBinaryError, ToolExecutionError, UserInputError
from .plugins import ToolCapabilities
from .resources import usable_cpus

NATIVE = "none"
DOCKER = "docker"
SINGULARITY = "singularity"


@dataclass
class ContainerConfig:
    backend: str = NATIVE
    engine: str | None = None  # explicit engine binary; defaults per backend
    platform: str | None = None  # e.g. "linux/amd64" for emulated BioContainers
    cache_dir: Path | None = None  # where Singularity .sif images and caches live
    wave_enabled: bool = False
    extra_mounts: tuple[Path, ...] = field(default_factory=tuple)

    @property
    def active(self) -> bool:
        return self.backend != NATIVE

    def cache_key(self) -> tuple:
        """Value-based identity for memoization keyed on this config.

        ``id()`` is unsuitable: a reconfigured global config frees the old
        object, whose address can be reused, yielding stale cache hits.
        """
        return (
            self.backend,
            self.engine,
            self.platform,
            str(self.cache_dir) if self.cache_dir else None,
            self.wave_enabled,
            tuple(str(m) for m in self.extra_mounts),
        )

    def engine_binary(self) -> str:
        if self.engine:
            return self.engine
        return "docker" if self.backend == DOCKER else "singularity"


_CONFIG = ContainerConfig()
# Keyed by (conda spec, platform): the same spec mints a different image per
# architecture, so an arm64 and an amd64 request must not share an entry.
_WAVE_CACHE: dict[tuple[tuple[str, ...], str | None], str] = {}


def configure_container(
    backend: str = NATIVE,
    *,
    engine: str | None = None,
    platform: str | None = None,
    cache_dir: str | os.PathLike[str] | None = None,
    wave_enabled: bool = False,
) -> ContainerConfig:
    """Set the process-global container backend. ``backend='none'`` = native."""
    global _CONFIG
    if backend not in (NATIVE, DOCKER, SINGULARITY):
        raise UserInputError(
            f"Unknown container backend '{backend}'; choose {NATIVE}, {DOCKER} or "
            f"{SINGULARITY} (use --container-engine for podman/apptainer binaries)."
        )
    _CONFIG = ContainerConfig(
        backend=backend,
        engine=engine,
        platform=platform,
        cache_dir=Path(cache_dir).resolve() if cache_dir else None,
        wave_enabled=wave_enabled,
    )
    return _CONFIG


def get_config() -> ContainerConfig:
    return _CONFIG


# Configurations whose engine answered: the check costs a round trip to the
# daemon, and a stage preflights several adapters.
_ENGINE_READY: set[tuple] = set()
_ENGINE_TIMEOUT = 60


def check_engine_ready(config: ContainerConfig | None = None) -> None:
    """Raise :class:`MissingBinaryError` when the Docker daemon cannot be reached.

    ``docker --version`` answers without a daemon, so a stopped Docker Desktop
    passed the binary check and every tool then failed with the engine's error
    hidden in the run log. ``docker info`` (also answered by podman) needs the
    daemon. Singularity and Apptainer have no daemon and are not asked.

    The check runs through :func:`subprocess.run`, not :func:`process.run`, so
    a termination signal is not forwarded to it. None is needed: the client
    stays in repgenr's process group, which a terminal Ctrl-C or a group-wide
    scheduler signal reaches directly, and on the main thread the handler's
    SystemExit interrupts the wait, after which subprocess.run kills the child.
    """
    config = config or _CONFIG
    if config.backend != DOCKER:
        return
    key = config.cache_key()
    if key in _ENGINE_READY:
        return
    engine = config.engine_binary()
    cmd = [engine, "info"]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=_ENGINE_TIMEOUT,
            stdin=subprocess.DEVNULL,
        )
    except subprocess.TimeoutExpired as exc:
        raise MissingBinaryError(
            f"The {engine} daemon did not answer '{engine} info' within "
            f"{_ENGINE_TIMEOUT} s. Start it (Docker Desktop on macOS) or run "
            "without --container."
        ) from exc
    except OSError as exc:
        raise MissingBinaryError(f"Could not run '{engine} info': {exc}") from exc
    if proc.returncode != 0:
        text = (proc.stderr or "") + (proc.stdout or "")
        detail = next((line.strip() for line in text.splitlines() if line.strip()), "")
        raise MissingBinaryError(
            f"The {engine} daemon is not reachable"
            + (f" ({detail})" if detail else "")
            + ". Start it (Docker Desktop on macOS) or run without --container."
        )
    _ENGINE_READY.add(key)


# Docker CPU counts per configuration (Docker Desktop runs tools in a VM whose
# CPU count can be lower than the host's).
_ENGINE_CPUS: dict[tuple, int] = {}


def available_cpus(caps: ToolCapabilities) -> int:
    """CPUs the tool sees where it runs: the host's, or the Docker engine's.

    Singularity shares the host's CPUs. When Docker does not answer, the host
    count is used.
    """
    # The CPUs this process may use: the affinity mask (a scheduler's cpuset,
    # taskset) capped by a cgroup CPU quota (docker --cpus, a Kubernetes
    # limit). A Docker engine's containers do not inherit repgenr's cgroup,
    # so the engine branch below asks the engine instead.
    host = usable_cpus()
    config = _CONFIG
    if config.backend != DOCKER or runs_on_host(caps):
        return host
    key = config.cache_key()
    if key not in _ENGINE_CPUS:
        count = host
        try:
            proc = subprocess.run(
                [config.engine_binary(), "info", "--format", "{{.NCPU}}"],
                capture_output=True,
                text=True,
                timeout=_ENGINE_TIMEOUT,
                stdin=subprocess.DEVNULL,
            )
            if proc.returncode == 0 and proc.stdout.strip().isdigit():
                count = max(1, int(proc.stdout.strip()))
        except (OSError, subprocess.SubprocessError):
            pass
        _ENGINE_CPUS[key] = count
    return _ENGINE_CPUS[key]


def resolve_image(caps: ToolCapabilities, config: ContainerConfig | None = None) -> str | None:
    """Return the image URI for an adapter, or None to run natively.

    With ``--wave`` and a ``caps.conda`` spec, an image is minted (and cached)
    from the spec via the Wave CLI: that is what the flag asks for, and it is
    the way to a native-architecture build. Otherwise the pinned
    ``caps.container`` is used. An adapter with neither runs on the host.
    """
    config = config or _CONFIG
    if config.wave_enabled and caps.conda:
        return _wave_image(caps.conda, config)
    if caps.container:
        return caps.container
    return None


def result_env_fragment() -> dict[str, list]:
    """Execution environment that can change a tool's result, for fingerprints.

    The backend, platform and Wave selection decide which tool builds run; the
    engine binary, cache directory and extra mounts are plumbing and stay out,
    so re-running elsewhere on the same machine still resumes.
    """
    config = get_config()
    return {"container": [config.backend, config.platform, config.wave_enabled]}


def runs_on_host(caps: ToolCapabilities) -> bool:
    """True when :func:`run_tool` would execute this adapter's commands natively.

    Adapters use it to gate checks that only make sense for the host PATH,
    such as looking for a particular build of a helper binary.
    """
    config = _CONFIG
    if not config.active:
        return True
    return resolve_image(caps, config) is None


_IMAGE_INSPECT_TIMEOUT = 30
# How docker ("No such image") and podman ("image not known") report an image
# that is not in the local store.
_IMAGE_ABSENT_PHRASES = ("no such image", "image not known")


def image_present(image: str, config: ContainerConfig | None = None) -> bool | None:
    """Is ``image`` available locally, so a run would not pull it first?

    Docker (or podman) is asked with ``image inspect``, which neither pulls
    nor contacts a registry. For Singularity the image is present when it is
    a local ``.sif`` file or its cached ``<cache_dir>/<name>.sif`` exists.
    None means the answer is not known: no cache directory for Singularity,
    or the engine could not be asked or did not say the image is absent
    (a daemon that is down).
    """
    config = config or _CONFIG
    if config.backend == SINGULARITY:
        if image.endswith(".sif"):
            return Path(image).is_file()
        if config.cache_dir is None:
            return None
        return (config.cache_dir / f"{_sanitize(image)}.sif").is_file()
    if config.backend != DOCKER:
        return None
    try:
        proc = subprocess.run(
            [config.engine_binary(), "image", "inspect", image],
            capture_output=True,
            text=True,
            timeout=_IMAGE_INSPECT_TIMEOUT,
            stdin=subprocess.DEVNULL,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode == 0:
        return True
    # Only an answer that the image is absent means "not pulled"; a daemon that
    # is down or refuses the request says nothing about the image.
    text = ((proc.stderr or "") + (proc.stdout or "")).lower()
    if any(phrase in text for phrase in _IMAGE_ABSENT_PHRASES):
        return False
    return None


def _wave_image(conda_spec: tuple[str, ...], config: ContainerConfig) -> str:
    cache_key = (conda_spec, config.platform)
    if cache_key in _WAVE_CACHE:
        return _WAVE_CACHE[cache_key]
    if shutil.which("wave") is None:
        raise MissingBinaryError(
            "Wave is enabled but the 'wave' CLI is not on PATH. Install it or pin "
            f"an explicit image for conda spec {' '.join(conda_spec)}."
        )
    cmd = ["wave"]
    for pkg in conda_spec:
        cmd += ["--conda-package", pkg]
    if config.platform:
        cmd += ["--platform", config.platform]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired as exc:
        raise ToolExecutionError(cmd, -1, output="wave timed out after 600s", timeout=600) from exc
    if proc.returncode != 0:
        raise ToolExecutionError(cmd, proc.returncode, output=proc.stderr.strip())
    lines = proc.stdout.strip().splitlines()
    if not lines or not lines[-1].strip():
        raise ToolExecutionError(
            cmd, proc.returncode, output="wave returned no image URI on stdout"
        )
    image = lines[-1].strip()
    _WAVE_CACHE[cache_key] = image
    return image


def _engine_env(config: ContainerConfig) -> dict[str, str]:
    """Env additions so Singularity/Apptainer caches live in cache_dir."""
    if config.backend != SINGULARITY or config.cache_dir is None:
        return {}
    cache = str(config.cache_dir)
    tmp = str(config.cache_dir / "tmp")
    Path(tmp).mkdir(parents=True, exist_ok=True)
    return {
        "SINGULARITY_CACHEDIR": cache,
        "APPTAINER_CACHEDIR": cache,
        "SINGULARITY_TMPDIR": tmp,
        "APPTAINER_TMPDIR": tmp,
    }


# Host variables that point a containerized tool at reference data (dRep runs
# CheckM, which reads CHECKM_DATA_PATH). When set, each is passed into the
# container and the directory it names is bound at the same path.
_FORWARDED_DATA_ENV = ("CHECKM_DATA_PATH",)


def _forwarded_data_env() -> dict[str, str]:
    # Absolute, because the container's working directory is not the host's.
    return {
        name: os.path.abspath(os.environ[name])
        for name in _FORWARDED_DATA_ENV
        if os.environ.get(name)
    }


def _default_mounts(
    config: ContainerConfig,
    cwd: str | os.PathLike[str] | None,
    argv: Sequence[str] = (),
    extra_mounts: Sequence[str | os.PathLike[str]] = (),
) -> list[Path]:
    # Use absolute (NOT symlink-resolved) paths: tools are given paths like
    # macOS /var/folders/... and /Users/...; resolving to /private/var or
    # /System/Volumes/Data would mount at paths that do not exist inside the
    # Linux container. Docker resolves host-side symlinks itself.
    def absp(p: str | os.PathLike[str]) -> Path:
        return Path(os.path.abspath(p))

    mounts: list[Path] = []
    # No explicit cwd means "the host process's working directory": relative
    # paths in argv (a stateless step run with `-o .`) must resolve there,
    # not in the temp mount.
    mounts.append(absp(cwd) if cwd is not None else absp(os.getcwd()))
    mounts.append(absp(tempfile.gettempdir()))
    mounts.extend(absp(m) for m in config.extra_mounts)
    for value in _forwarded_data_env().values():
        if Path(value).exists():
            mounts.append(absp(value))
    # Per-call mounts for inputs referenced indirectly (e.g. genome paths listed
    # inside a manifest file rather than passed as argv tokens).
    for m in extra_mounts:
        d = Path(m)
        d = d if d.is_dir() else d.parent
        if d.exists():
            mounts.append(absp(d))
    # Bind the directories referenced by absolute-path arguments (genome inputs,
    # output dirs, references), so the tool sees its files at identical paths.
    for token in argv:
        if token.startswith("/"):
            p = Path(token)
            d = p if p.is_dir() else p.parent
            if d.exists():
                mounts.append(absp(d))
    # Genomes staged by `ingest` (and outgroup files) are symlinks whose
    # targets can live anywhere on the host; a bind of the link's directory
    # alone leaves a dangling link inside the container. Bind each target's
    # directory as well (one level: the entries of every mounted directory).
    for d in list(mounts):
        if not d.is_dir():
            continue
        try:
            entries = list(d.iterdir())
        except OSError:
            continue
        for entry in entries:
            if entry.is_symlink():
                target = Path(os.path.abspath(os.path.join(d, os.readlink(entry))))
                if target.parent.exists():
                    mounts.append(target.parent)
    # A path may pass through a symlinked directory (Nextflow stages an input
    # directory as a link inside the task dir: chunks/c0 -> /work/ab/c0). The
    # container resolves that link to the real directory, which the mount of
    # the link's path does not cover, so bind the real directory at its own
    # path as well. The un-resolved mount stays, for tools given the link path.
    for d in list(mounts):
        real = Path(os.path.realpath(d))
        if real != d and real.exists():
            mounts.append(real)
    # de-duplicate, dropping any mount nested under another
    unique: list[Path] = []
    for m in sorted(set(mounts), key=lambda p: len(str(p))):
        if not any(_is_relative_to(m, u) for u in unique):
            unique.append(m)
    return unique


def _is_relative_to(child: Path, parent: Path) -> bool:
    try:
        child.relative_to(parent)
        return True
    except ValueError:
        return False


def _sanitize(image: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", image)


def _singularity_source(image: str, config: ContainerConfig, logger: logging.Logger) -> str:
    """Return a Singularity image reference, pre-pulling a .sif into cache_dir.

    A local ``.sif`` is used as-is. Otherwise the image is treated as a registry
    reference (``docker://`` is assumed when no scheme is present). With a cache
    dir set, the image is pulled once to ``<cache>/<name>.sif`` and reused.
    """
    if image.endswith(".sif"):
        return image
    source = image if "://" in image else f"docker://{image}"
    if config.cache_dir is None:
        return source
    config.cache_dir.mkdir(parents=True, exist_ok=True)
    sif = config.cache_dir / f"{_sanitize(image)}.sif"
    if not sif.exists():
        logger.info("Pulling %s -> %s", source, sif)
        # Pull to a temp path and rename on success: an interrupted pull must not
        # leave a corrupt .sif that sif.exists() then treats as cached forever.
        tmp = sif.with_name(sif.name + ".partial")
        tmp.unlink(missing_ok=True)
        try:
            process.run(
                [config.engine_binary(), "pull", str(tmp), source],
                logger=logger,
                env=_engine_env(config),
                log_prefix="singularity",
            )
        except Exception:
            tmp.unlink(missing_ok=True)
            raise
        tmp.replace(sif)
    return str(sif)


def wrap_command(
    image: str,
    argv: Sequence[str],
    *,
    config: ContainerConfig,
    cwd: str | os.PathLike[str] | None,
    logger: logging.Logger,
    extra_mounts: Sequence[str | os.PathLike[str]] = (),
    stdin: bool = False,
    tool_env: Mapping[str, str] | None = None,
) -> list[str]:
    """Build the engine command that runs ``argv`` inside ``image``.

    ``stdin`` keeps the client's standard input open for the tool (the
    consumer of :func:`run_pipeline`); Docker closes it otherwise.
    ``tool_env`` holds variables the tool itself must see (a thread count
    such as ``OMP_NUM_THREADS``): Docker gets them as ``-e``; Singularity
    passes the client's environment on, so the runners put them there.
    """
    mounts = _default_mounts(config, cwd, argv, extra_mounts)
    workdir = str(Path(os.path.abspath(cwd if cwd is not None else os.getcwd())))

    if config.backend == DOCKER:
        # --init: an init process becomes process 1 and passes signals on, so
        # the SIGTERM the client forwards when repgenr stops a tool ends it
        # (a tool running as process 1 ignores a signal it has no handler for).
        cmd = [config.engine_binary(), "run", "--rm", "--init", "--entrypoint", ""]
        if stdin:
            cmd.append("-i")
        # A name lets repgenr stop the container when it is stopped itself
        # (see _stop_container).
        cmd += ["--name", f"repgenr-{os.getpid()}-{uuid.uuid4().hex[:12]}"]
        # Labels find the containers of a repgenr that was killed with SIGKILL,
        # which no handler sees (a label filter matches exactly, a name
        # filter matches substrings):
        #   docker ps -q --filter label=repgenr.pid=<pid> | xargs -r docker stop
        cmd += ["--label", f"repgenr.pid={os.getpid()}"]
        cmd += ["--label", f"repgenr.host={socket.gethostname()}"]
        cmd += ["-u", f"{os.getuid()}:{os.getgid()}"]
        # Run as an arbitrary host UID with no passwd entry, so HOME defaults to
        # "/" and is not writable. Point it at the mounted, writable workdir so
        # tools that touch HOME (e.g. Toil/Cactus creating its config dir) work.
        cmd += ["-e", f"HOME={workdir}"]
        for name, value in _forwarded_data_env().items():
            cmd += ["-e", f"{name}={value}"]
        for name, value in (tool_env or {}).items():
            cmd += ["-e", f"{name}={value}"]
        if config.platform:
            cmd += ["--platform", config.platform]
        for m in mounts:
            cmd += ["-v", f"{m}:{m}"]
        cmd += ["-w", workdir, image, *argv]
        return cmd

    # singularity / apptainer
    source = _singularity_source(image, config, logger)
    cmd = [config.engine_binary(), "exec"]
    for m in mounts:
        cmd += ["--bind", str(m)]
    cmd += ["--pwd", workdir, source, *argv]
    return cmd


def _stop_container(wrapped: Sequence[str], config: ContainerConfig) -> None:
    """Stop the named Docker container of an interrupted ``docker run``.

    When repgenr is stopped it forwards SIGTERM to the docker client, which
    passes it through --init to the tool. A tool that ignores SIGTERM kept its
    container running after repgenr and the client had exited. The forwarded
    SIGTERM was the tool's notice, so ``docker stop -t 0`` kills it at
    once: a scheduler that sends SIGKILL a few seconds after SIGTERM must not
    find the container still waiting out a second grace period.
    Errors are ignored: the container may already be gone.

    ``docker stop`` runs outside the registry of :func:`process.run`; it is
    the cleanup itself, and when a second signal ends repgenr while it runs,
    the client runs on and completes the stop.
    """
    name = _container_name(wrapped, config)
    if name is None:
        return
    try:
        subprocess.run(
            [config.engine_binary(), "stop", "-t", "0", name],
            capture_output=True,
            timeout=20,
            stdin=subprocess.DEVNULL,
        )
    except (OSError, subprocess.SubprocessError):
        pass


def _container_name(wrapped: Sequence[str], config: ContainerConfig) -> str | None:
    """The ``--name`` of a wrapped ``docker run`` command, else None."""
    if config.backend != DOCKER or "--name" not in wrapped:
        return None
    return wrapped[list(wrapped).index("--name") + 1]


# Names of the Docker containers whose ``docker run`` client is running now,
# guarded by process._live_lock (reentrant, also taken by the signal handler).
# A second signal ends repgenr without Python's cleanup, so the final stop hook
# stops these containers instead.
_live_names: dict[str, str] = {}  # container name -> engine binary


def _register_container(wrapped: Sequence[str], config: ContainerConfig) -> str | None:
    name = _container_name(wrapped, config)
    if name is not None:
        with process._live_lock:
            _live_names[name] = config.engine_binary()
    return name


def _unregister_container(name: str | None) -> None:
    if name is not None:
        with process._live_lock:
            _live_names.pop(name, None)


def _stop_live_containers() -> None:
    """Start one detached ``docker stop -t 0`` for every live container.

    Runs from the termination handler on a second signal, just before repgenr
    exits: it starts the client in its own session and does not wait, so the
    stop completes after repgenr has gone.
    """
    with process._live_lock:
        by_engine: dict[str, list[str]] = {}
        for name, engine in _live_names.items():
            by_engine.setdefault(engine, []).append(name)
    for engine, names in by_engine.items():
        try:
            subprocess.Popen(
                [engine, "stop", "-t", "0", *names],
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
        except OSError:
            pass


process.register_final_stop_hook(_stop_live_containers)


def _stop_after_failure(
    exc: ToolExecutionError, wrapped: Sequence[str], config: ContainerConfig
) -> None:
    """Stop the container of a ``docker run`` whose client was killed.

    The client is killed while repgenr stops, and on a timeout (the ``timeout``
    argument or REPGENR_SUBPROCESS_TIMEOUT). Killing the client does not stop
    the container, so a tool that ignores the forwarded SIGTERM, or never got
    it, would keep running. On a pool thread the termination handler does not
    interrupt the wait: the escalation timer kills the client and process.run
    raises ToolExecutionError rather than SystemExit. A tool that was never
    started has no container to stop.
    """
    if exc.timeout is not None or (
        process.stop_requested.is_set() and exc.output != process.NOT_STARTED
    ):
        _stop_container(wrapped, config)


# `docker run` exits 125 when the engine itself fails (image not found, pull
# denied, bad option): the tool never started.
_DOCKER_ENGINE_FAILURE = 125
_ENGINE_ERROR = re.compile(r"^(docker|podman): |Error response from daemon|Unable to find image")


def _engine_failure(
    exc: ToolExecutionError, image: str, caps: ToolCapabilities, config: ContainerConfig
) -> ToolExecutionError:
    """Name a failure of the engine itself instead of blaming the tool."""
    if (
        config.backend != DOCKER
        or exc.returncode != _DOCKER_ENGINE_FAILURE
        or exc.timeout is not None
    ):
        return exc
    lines = [line.strip() for line in (exc.output or "").splitlines() if line.strip()]
    # A tool may exit 125 itself; only the engine's own error marks a failure
    # to start the image.
    engine_lines = [line for line in lines if _ENGINE_ERROR.search(line)]
    if not engine_lines:
        return exc
    detail = engine_lines[-1]
    message = f"{config.engine_binary()} could not start image {image} for {caps.name} (exit 125)"
    if detail:
        message += f": {detail}"
    return ToolExecutionError(
        exc.command, exc.returncode, output=exc.output, tool=exc.tool, message=message
    )


def engine_failure_line(returncode: int, lines: Sequence[str]) -> str | None:
    """The engine's own error when a Docker client exited 125 without
    starting its tool (image not found, daemon down), else None.

    For callers that read exit codes instead of catching
    :class:`ToolExecutionError` (:func:`run_pipeline`): such a failure is not
    the tool's and not one run's, so a batch stops on it instead of
    recording it against every run.
    """
    if _CONFIG.backend != DOCKER or returncode != _DOCKER_ENGINE_FAILURE:
        return None
    hits = [line.strip() for line in lines if _ENGINE_ERROR.search(line.strip())]
    return hits[-1] if hits else None


def ensure_image(caps: ToolCapabilities, logger: logging.Logger) -> str | None:
    """Pull the image of ``caps`` once, before a batch runs it many times.

    Returns the image reference, or None when no backend is active or the
    tool has no image (it then runs on the host). Under Docker an absent
    image is pulled with ``docker pull``; under Singularity the pull into
    the cache directory happens here instead of at the first run. A batch of
    concurrent runs otherwise each pulls the image inside its own measured
    time, and a pull failure is reported once, by name, instead of once per
    run.
    """
    config = _CONFIG
    if not config.active:
        return None
    image = resolve_image(caps, config)
    if image is None:
        return None
    if config.backend == SINGULARITY:
        _singularity_source(image, config, logger)
        return image
    if config.backend == DOCKER and image_present(image, config) is False:
        logger.info("Pulling %s for %s", image, caps.name)
        cmd = [config.engine_binary(), "pull"]
        if config.platform:
            # The platform of wrap_command, so the variant pulled is the one run.
            cmd += ["--platform", config.platform]
        cmd.append(image)
        try:
            process.run(
                cmd,
                logger=logger,
                env=_engine_env(config),
                log_prefix=config.engine_binary(),
                console_command=cmd,
            )
        except ToolExecutionError as exc:
            raise ToolExecutionError(
                exc.command,
                exc.returncode,
                output=exc.output,
                tool=caps.name,
                message=f"{config.engine_binary()} could not pull image {image} for {caps.name}",
            ) from exc
    return image


def run_tool(
    caps: ToolCapabilities,
    command: Sequence[str | os.PathLike[str]],
    *,
    logger: logging.Logger,
    cwd: str | os.PathLike[str] | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = True,
    stdout_path: str | os.PathLike[str] | None = None,
    log_prefix: str | None = None,
    extra_mounts: Sequence[str | os.PathLike[str]] = (),
    timeout: float | None = None,
    tool_env: Mapping[str, str] | None = None,
) -> int:
    """Run an adapter's tool command, containerized when a backend is active.

    ``env`` is the client's environment (a --bin-dir PATH, the engine's cache
    variables); ``tool_env`` reaches the tool itself in every backend, which
    ``env`` does not under Docker (it goes to the ``docker run`` client).
    """
    config = _CONFIG
    image = resolve_image(caps, config) if config.active else None
    if image is None:
        env = _with_tool_env(env, tool_env)
        if config.active:
            hint = (
                "pass --wave to mint one from its conda spec"
                if caps.conda
                else "it declares neither a container image nor a conda spec"
            )
            logger.warning(
                "Tool '%s' has no container image to run in (%s); running on the "
                "host despite --container %s.",
                caps.name,
                hint,
                config.backend,
            )
        from .bindirs import note_used, tool_path

        path = tool_path(caps.name, base=(env or {}).get("PATH"))
        if path is not None:
            # --bin-dir: this tool's directory comes first for it and its helpers.
            env = {**(env or {}), "PATH": path}
            note_used(caps.name)
        return process.run(
            command,
            logger=logger,
            cwd=cwd,
            env=env,
            check=check,
            stdout_path=stdout_path,
            log_prefix=log_prefix,
            timeout=timeout,
        )

    argv = [str(part) for part in command]
    wrapped = wrap_command(
        image,
        argv,
        config=config,
        cwd=cwd,
        logger=logger,
        extra_mounts=extra_mounts,
        tool_env=tool_env,
    )
    merged_env = _client_env(config, env, tool_env)
    name = _register_container(wrapped, config)
    try:
        returncode = process.run(
            wrapped,
            logger=logger,
            cwd=cwd,
            env=merged_env,
            check=check,
            stdout_path=stdout_path,
            log_prefix=log_prefix or caps.name,
            timeout=timeout,
            # The console shows the tool's argv; the engine preamble (mounts,
            # labels, user) is in the run log.
            console_command=argv,
        )
        if returncode != 0 and process.stop_requested.is_set():
            # check=False: the killed client returned instead of raising.
            _stop_container(wrapped, config)
        return returncode
    except ToolExecutionError as exc:
        _stop_after_failure(exc, wrapped, config)
        raise _engine_failure(exc, image, caps, config) from exc
    except (SystemExit, KeyboardInterrupt):
        _stop_container(wrapped, config)
        raise
    finally:
        _unregister_container(name)


def _with_tool_env(
    env: Mapping[str, str] | None, tool_env: Mapping[str, str] | None
) -> Mapping[str, str] | None:
    """``env`` with the tool's own variables added (a host or Singularity run,
    where the process environment is the tool's)."""
    if not tool_env:
        return env
    return {**(env or {}), **tool_env}


def _client_env(
    config: ContainerConfig,
    env: Mapping[str, str] | None,
    tool_env: Mapping[str, str] | None,
) -> Mapping[str, str] | None:
    """The environment of the engine client: the engine's cache variables and
    ``env``; under Singularity also ``tool_env``, which the container inherits
    (Docker gets it as ``-e`` from :func:`wrap_command` instead)."""
    merged = {**_engine_env(config), **(dict(env) if env else {})}
    if config.backend != DOCKER and tool_env:
        merged.update(tool_env)
        if config.backend == SINGULARITY:
            # Also as the engine's own injection variable, which survives a
            # --cleanenv set through the environment (one prefix, the
            # engine's own, so it does not warn about a duplicate).
            prefix = "APPTAINERENV_" if "apptainer" in config.engine_binary() else "SINGULARITYENV_"
            for name, value in tool_env.items():
                merged[f"{prefix}{name}"] = value
    return merged or None


def _host_env(
    caps: ToolCapabilities,
    env: Mapping[str, str] | None,
    tool_env: Mapping[str, str] | None = None,
) -> Mapping[str, str] | None:
    """``env`` with this tool's --bin-dir first on PATH, when one is set, and
    its ``tool_env``."""
    from .bindirs import note_used, tool_path

    env = _with_tool_env(env, tool_env)
    path = tool_path(caps.name, base=(env or {}).get("PATH"))
    if path is None:
        return env
    note_used(caps.name)
    return {**(env or {}), "PATH": path}


def run_pipeline(
    producer: tuple[ToolCapabilities, Sequence[str | os.PathLike[str]]],
    consumer: tuple[ToolCapabilities, Sequence[str | os.PathLike[str]]],
    *,
    logger: logging.Logger,
    cwd: str | os.PathLike[str] | None = None,
    env: Mapping[str, str] | None = None,
    extra_mounts: Sequence[str | os.PathLike[str]] = (),
    timeout: float | None = None,
    tool_envs: tuple[Mapping[str, str] | None, Mapping[str, str] | None] = (None, None),
) -> process.PipelineResult:
    """Run ``producer | consumer``, two tools joined by a pipe, and return both
    exit codes with the last stderr lines of each (nothing is raised for a
    non-zero exit; see :func:`process.run_pipeline`).

    Each half is ``(capabilities, argv)``. On the host the pipe is the plain
    one; each half gets its own --bin-dir on PATH. Under a container backend
    each half runs in its own image with the host pipe between them, the
    consumer's client keeping stdin open, so the single image per tool rule
    holds and no combined image is needed. A half without an image runs on
    the host, as :func:`run_tool` does. ``tool_envs`` holds the variables
    each tool itself must see, as ``tool_env`` of :func:`run_tool`.
    """
    config = _CONFIG
    halves: list[list[str]] = []
    shown: list[list[str]] = []
    names: list[str | None] = []
    # One environment per half: a --bin-dir applies to its own tool only.
    envs: list[Mapping[str, str] | None] = []
    for (caps, command), stdin, tool_env in (
        (producer, False, tool_envs[0]),
        (consumer, True, tool_envs[1]),
    ):
        argv = [str(part) for part in command]
        image = resolve_image(caps, config) if config.active else None
        if image is None:
            halves.append(argv)
            envs.append(_host_env(caps, env, tool_env))
            names.append(None)
        else:
            wrapped = wrap_command(
                image,
                argv,
                config=config,
                cwd=cwd,
                logger=logger,
                extra_mounts=extra_mounts,
                stdin=stdin,
                tool_env=tool_env,
            )
            halves.append(wrapped)
            envs.append(_client_env(config, env, tool_env))
            names.append(_register_container(wrapped, config))
        shown.append(argv)
    try:
        result = process.run_pipeline(
            halves[0],
            halves[1],
            logger=logger,
            cwd=cwd,
            envs=(envs[0], envs[1]),
            log_prefixes=(producer[0].name, consumer[0].name),
            timeout=timeout,
            console_commands=(shown[0], shown[1]),
        )
        if not result.ok:
            # A killed or broken client leaves its container running (the
            # producer after a consumer failure, both after a timeout or a
            # stop request); stop them as run_tool does.
            for wrapped in halves:
                _stop_container(wrapped, config)
        return result
    except (SystemExit, KeyboardInterrupt):
        for wrapped in halves:
            _stop_container(wrapped, config)
        raise
    finally:
        for name in names:
            _unregister_container(name)


def run_chain(
    caps: ToolCapabilities,
    steps: Sequence[tuple[str, Sequence[str | os.PathLike[str]]]],
    *,
    logger: logging.Logger,
    cwd: str | os.PathLike[str] | None = None,
    env: Mapping[str, str] | None = None,
    extra_mounts: Sequence[str | os.PathLike[str]] = (),
    timeout: float | None = None,
    tool_env: Mapping[str, str] | None = None,
) -> None:
    """Run a sequence of commands that form one unit of work.

    Each step is ``(log_prefix, argv)`` and they run in order, stopping at the
    first failure. On the host this is exactly a loop over :func:`run_tool`.
    Containerized, the whole sequence runs inside a single container: starting
    one per command costs seconds each under an emulated engine, and a
    per-genome chain is half a dozen commands.

    The steps must be plain argument vectors. Redirection and pipes are not
    supported, so a tool that writes to stdout needs its own output flag.
    """
    config = _CONFIG
    image = resolve_image(caps, config) if config.active else None
    if image is None:
        for prefix, command in steps:
            run_tool(
                caps,
                command,
                logger=logger,
                cwd=cwd,
                env=env,
                log_prefix=prefix,
                extra_mounts=extra_mounts,
                timeout=timeout,
                tool_env=tool_env,
            )
        return

    argvs = [[str(part) for part in command] for _, command in steps]
    for (prefix, _), argv in zip(steps, argvs, strict=True):
        process.log_command(logger, argv, f"[{prefix}] ")
    # The commands live inside a script string, so the mount scan cannot see
    # their paths in argv; hand it every absolute path they name.
    paths = [token for argv in argvs for token in argv if token.startswith("/")]
    script = "set -e\n" + "\n".join(shlex.join(argv) for argv in argvs) + "\n"
    wrapped = wrap_command(
        image,
        ["sh", "-c", script],
        config=config,
        cwd=cwd,
        logger=logger,
        extra_mounts=[*extra_mounts, *paths],
        tool_env=tool_env,
    )
    merged_env = _client_env(config, env, tool_env)
    name = _register_container(wrapped, config)
    try:
        process.run(
            wrapped,
            logger=logger,
            cwd=cwd,
            env=merged_env,
            log_prefix=caps.name,
            timeout=timeout,
            # Each step was logged above; the console names the script only.
            console_command=["sh", "-c", f"<{len(argvs)} commands above>"],
        )
    except ToolExecutionError as exc:
        _stop_after_failure(exc, wrapped, config)
        raise _engine_failure(exc, image, caps, config) from exc
    except (SystemExit, KeyboardInterrupt):
        _stop_container(wrapped, config)
        raise
    finally:
        _unregister_container(name)


def run_tool_with_retries(
    caps: ToolCapabilities,
    command: Sequence[str | os.PathLike[str]],
    *,
    logger: logging.Logger,
    attempts: int = 3,
    retry_delay: float = 5.0,
    permanent: Callable[[ToolExecutionError], bool] | None = None,
    **kwargs,
) -> int:
    """Run a network-dependent tool with exponential-backoff retries.

    For tools like the NCBI datasets CLI that perform their own transfers and
    have no built-in retry: a transient failure is retried (backoff doubles per
    attempt); the last failure propagates unchanged. ``permanent`` marks a
    failure that a retry cannot fix (e.g. "no such accession"), which then
    propagates at once.
    """
    for attempt in range(1, attempts + 1):
        try:
            return run_tool(caps, command, logger=logger, **kwargs)
        except ToolExecutionError as exc:
            if attempt == attempts or (permanent is not None and permanent(exc)):
                raise
            delay = retry_delay * 2 ** (attempt - 1)
            logger.warning(
                "%s failed (attempt %d/%d): %s; retrying in %.0fs",
                caps.name,
                attempt,
                attempts,
                exc,
                delay,
            )
            time.sleep(delay)
    raise RuntimeError("unreachable")
