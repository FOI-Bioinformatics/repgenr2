"""External binary discovery and version preflight.

Every adapter declares the binaries it needs via :class:`BinarySpec` in its
``ToolCapabilities``. Before a tool runs, :func:`check_binaries` confirms each
binary is on PATH and, when a minimum version is declared, that it is new
enough. Resolved versions are returned so the stage can record them in
``repgenr.yaml`` provenance.
"""

from __future__ import annotations

import contextlib
import json
import logging
import os
import re
import shutil
import subprocess
from collections.abc import Iterator
from contextvars import ContextVar
from dataclasses import dataclass, field
from pathlib import Path

from . import process
from .errors import MissingBinaryError

# Seconds a version query may take. A cold start (a Java tool, a conda
# wrapper on a network disk) can take several seconds; a query that has not
# answered by then is stopped and its version is recorded as unknown.
VERSION_TIMEOUT = 30.0
_timeout_override: ContextVar[float | None] = ContextVar("version_timeout", default=None)
# Seconds a stopped version query gets between SIGTERM and SIGKILL.
_STOP_GRACE = 1.0

# A version is not part of a longer token: "python3.12" in an interpreter path
# or "GLIBC_2.17" in a loader error is not the tool's version. A leading "v"
# ("dRep v3.4.5", "hal2maf v2.2") is allowed.
_VERSION_RE = re.compile(r"(?<![0-9A-UW-Za-uw-z_./])(\d+)\.(\d+)(?:\.(\d+))?")
# Terminal colour and style codes. OrthoFinder 3 colours its version even
# when the output is not a terminal, splitting "3.1.5" into "3." and "1.5".
_ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


@dataclass(frozen=True)
class BinarySpec:
    """A required external executable.

    ``version_args`` is the argument vector that prints a version (e.g.
    ``("--version",)``). ``min_version`` is an optional ``(major, minor, patch)``
    or dotted string requirement.
    """

    name: str
    version_args: tuple[str, ...] = field(default=("--version",))
    min_version: str | None = None
    # When the version cannot be parsed, lenient specs pass (avoid false rejects);
    # strict specs fail -- use this for tools where a legacy build that does not
    # answer its version flag can shadow a modern one on PATH (samtools/bcftools).
    strict_version: bool = False


def _parse_version(text: str) -> tuple[int, int, int] | None:
    match = _VERSION_RE.search(_ANSI_RE.sub("", text))
    if not match:
        return None
    major, minor, patch = match.group(1), match.group(2), match.group(3)
    return (int(major), int(minor), int(patch or 0))


@contextlib.contextmanager
def version_timeout(seconds: float) -> Iterator[None]:
    """Use ``seconds`` as the version-query timeout of :func:`check_binaries` calls
    made inside the block (``list-tools --check`` asks with a short one)."""
    token = _timeout_override.set(seconds)
    try:
        yield
    finally:
        _timeout_override.reset(token)


def _ask(argv: list[str], timeout: float, path: str | None = None) -> tuple[int, str] | None:
    """Run a version query; return its exit status and combined output.

    ``path`` (a per-tool PATH, see core.bindirs) is also the query's PATH, so
    a wrapper script answers with the helpers of its own environment.

    The query runs in its own session, so on timeout the whole process group
    is stopped: a shell wrapper's helper that holds the output pipe open
    would otherwise keep the read waiting after the wrapper was killed.
    Returns None when it could not be run or did not answer in time.
    """
    try:
        proc = subprocess.Popen(
            argv,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            # A tool queried without arguments (FastTree) may otherwise wait
            # on an inherited terminal for input until the timeout.
            stdin=subprocess.DEVNULL,
            text=True,
            start_new_session=True,
            env={**os.environ, "PATH": path} if path is not None else None,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    # Registered while it runs, so a termination signal to repgenr stops the
    # query's group too (it is in its own session and gets no terminal SIGINT).
    process.register_live(proc)
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        process.stop_group(proc, grace=_STOP_GRACE)
        try:
            proc.communicate(timeout=5)
        except (subprocess.TimeoutExpired, OSError, ValueError):
            # A helper that left the group still holds the pipe; stop reading.
            for stream in (proc.stdout, proc.stderr):
                if stream is not None:
                    stream.close()
        logging.getLogger("repgenr").warning(
            "%s did not answer within %g s; its version is unknown unless its "
            "conda package record names one.",
            " ".join([Path(argv[0]).name, *argv[1:]]),
            timeout,
        )
        return None
    except BaseException:
        # Ctrl-C or the termination handler's SystemExit: the query and its
        # helpers are in their own session, so stop them before re-raising.
        process.stop_group(proc, grace=0)
        raise
    finally:
        process.unregister_live(proc)
    return proc.returncode, (out or "") + (err or "")


def _query_version(
    name: str,
    version_args: tuple[str, ...],
    timeout: float | None = None,
    path: str | None = None,
) -> str | None:
    answer = _ask([name, *version_args], VERSION_TIMEOUT if timeout is None else timeout, path)
    if answer is None:
        return None
    returncode, blob = answer
    # A tool that crashed while answering (a Python traceback, exit != 0) has
    # not reported a version; numbers in its traceback are file paths.
    if returncode != 0 and "Traceback (most recent call last)" in blob:
        return None
    parsed = _parse_version(blob)
    if parsed:
        return ".".join(map(str, parsed))
    # Without a version number, keep the first line only when the tool accepted
    # the flag: a non-zero exit means the line is an error message (sibeliaz
    # has no version flag and answers `-v` with "illegal option"), not a version.
    if returncode != 0:
        return None
    return blob.strip().splitlines()[0] if blob.strip() else None


def _owns_binary(data: object, name: str) -> bool:
    """True when a conda-meta record lists ``bin/<name>`` among its files."""
    return isinstance(data, dict) and f"bin/{name}" in (data.get("files") or ())


def _metadata_version(name: str, path: str | None = None) -> str | None:
    """The version of the conda package that installed ``name``, if any.

    Used when the binary does not report a version itself (SibeliaZ has no
    version flag). The package records sit in ``<prefix>/conda-meta`` next to
    the ``bin`` directory that holds the binary; the path is not resolved,
    since conda links files into the prefix. The record named after the binary
    is tried first, then any record whose file list contains ``bin/<name>``.
    """
    found = path or shutil.which(name)
    if found is None:
        return None
    meta = Path(found).parent.parent / "conda-meta"
    try:
        if not meta.is_dir():
            return None
        named = sorted(meta.glob(f"{name}-*.json"))
        others = [p for p in sorted(meta.glob("*.json")) if p not in set(named)]
        for record in (*named, *others):
            try:
                text = record.read_text(encoding="utf-8")
                if f'"bin/{name}"' not in text:
                    continue
                data = json.loads(text)
            except (OSError, ValueError):
                continue
            if _owns_binary(data, name):
                version = data.get("version")
                if isinstance(version, str) and version.strip():
                    return version.strip()
    except OSError:
        return None
    return None


def check_binaries(
    specs: tuple[BinarySpec, ...], *, timeout: float | None = None, path: str | None = None
) -> dict[str, str]:
    """Confirm all ``specs`` are present (and new enough). Return name -> version.

    Each version query may take ``timeout`` seconds; None means the value set
    by :func:`version_timeout`, or :data:`VERSION_TIMEOUT`. A query that does
    not answer in time is stopped and logged, and its version is unknown.
    ``path`` replaces PATH for the lookup and the version query (a tool given
    its own directory with ``--bin-dir``); None uses PATH.
    Raises :class:`MissingBinaryError` listing every missing or too-old binary.
    """
    if timeout is None:
        override = _timeout_override.get()
        timeout = VERSION_TIMEOUT if override is None else override
    versions: dict[str, str] = {}
    problems: list[str] = []

    for spec in specs:
        found = shutil.which(spec.name, path=path)
        if found is None:
            problems.append(f"{spec.name}: not found on PATH")
            continue

        if path is None:
            reported = _query_version(spec.name, spec.version_args, timeout=timeout)
        else:
            reported = _query_version(found, spec.version_args, timeout=timeout, path=path)
        if _parse_version(reported or "") is None:
            # No version number from the tool itself: the conda package
            # record, when there is one, still names the installed version.
            reported = _metadata_version(spec.name, found) or reported
        versions[spec.name] = reported or "unknown"

        if spec.min_version is not None:
            have = _parse_version(reported or "")
            want = _parse_version(spec.min_version)
            if have is None and want is not None and spec.strict_version:
                got = f"got {reported!r}" if reported else "no version reported"
                problems.append(
                    f"{spec.name}: could not read a version ({got}); need "
                    f">= {spec.min_version}. A legacy build that does not answer its "
                    "version flag may be shadowing it on PATH -- use a modern environment."
                )
            elif have is not None and want is not None and have < want:
                problems.append(f"{spec.name}: version {reported} < required {spec.min_version}")

    if problems:
        raise MissingBinaryError(
            "Required external tools are missing or outdated:\n  " + "\n  ".join(problems)
        )
    return versions
