"""Safe subprocess execution.

Replaces the old ``subprocess.call(' '.join(cmd), shell=True)`` pattern and the
shell globs (``genomes/*.fasta``) that break past ``ARG_MAX`` for large genome
sets. Commands are always passed as argument vectors (no shell), output is
streamed to the logger, and a non-zero exit raises :class:`ToolExecutionError`.

For tools that genuinely need a large list of input files, write the list to a
file-of-filenames with :func:`write_fofn` and pass that path instead of a glob.
"""

from __future__ import annotations

import gzip
import logging
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import zipfile
import zlib
from collections import deque
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

from .contracts import atomic_path, record_name
from .errors import ToolExecutionError, UserInputError, WorkdirError
from .logging import console_extra, console_room

if TYPE_CHECKING:
    from .plugins import ToolCapabilities

_DEFAULT_TAIL = 50

_module_logger = logging.getLogger(__name__)

# A malformed REPGENR_SUBPROCESS_TIMEOUT is reported once, not on every
# subprocess launch.
_warned_bad_timeout = False

# Every tool starts as the leader of its own session and process group
# (POSIX), so it can be stopped together with the helpers it starts (RAxML
# under run_gubbins.py, skani under skDER, minimap2 under a typer). The tool
# then no longer receives the terminal's Ctrl-C or hangup, so repgenr's own
# handlers (install_termination_handler) forward them.
_OWN_GROUP = os.name == "posix"

# After SIGTERM, seconds a tool group gets to exit before it is sent SIGKILL.
STOP_GRACE_SECONDS = 5.0

# Tools running now, so a termination signal can stop them before repgenr
# exits. Several run at once under parallel_map. Reentrant: the signal handler
# runs on the main thread, which may hold the lock inside run() when the
# signal arrives.
_live_lock = threading.RLock()
_live: set[subprocess.Popen[Any]] = set()
# Set by the termination handler: run() then starts no further tool, so tasks
# still queued in a thread pool do not launch while repgenr shuts down.
stop_requested = threading.Event()
# Output of the ToolExecutionError that run() raises for a tool it did not
# start because repgenr is stopping.
NOT_STARTED = "not started: repgenr is stopping"
# Called by the termination handler on a second signal, just before repgenr
# takes the default action (see register_final_stop_hook).
_final_stop_hooks: list[Callable[[], None]] = []


def register_final_stop_hook(hook: Callable[[], None]) -> None:
    """Call ``hook`` when a second signal ends repgenr at once.

    A second signal skips Python's cleanup (``finally`` blocks, the container
    stop in :mod:`repgenr.core.containers`), so a resource that outlives the
    process, such as a container managed by a daemon, needs this hook. The hook
    runs inside the signal handler on the main thread: it must be quick and
    must not wait for a child. An exception it raises is ignored. Registering
    the same hook twice has no effect.
    """
    if hook not in _final_stop_hooks:
        _final_stop_hooks.append(hook)


def _run_final_stop_hooks() -> None:
    for hook in list(_final_stop_hooks):
        try:
            hook()
        except Exception:  # repgenr is exiting; there is nothing to report to
            pass


def _signal_group(proc: subprocess.Popen[Any], sig: int) -> bool:
    """Send ``sig`` to the tool's process group (or the tool alone off POSIX).

    The group is signalled even when the tool itself has exited, because a
    helper it started can outlive it. Returns False when nothing was there.
    """
    try:
        if _OWN_GROUP:
            os.killpg(proc.pid, sig)
        elif proc.poll() is None:
            proc.send_signal(sig)
        else:
            return False
    except (ProcessLookupError, PermissionError, OSError):
        return False
    return True


def stop_group(proc: subprocess.Popen[Any], grace: float | None = None) -> None:
    """Stop a tool and its helpers: SIGTERM, wait up to ``grace`` s, then SIGKILL.

    The tool must lead its own process group (``start_new_session=True``), as
    tools started by :func:`run` and version queries do.
    """
    _signal_group(proc, signal.SIGTERM)
    try:
        proc.wait(timeout=STOP_GRACE_SECONDS if grace is None else grace)
    except subprocess.TimeoutExpired:
        pass
    # Also reaches helpers that ignored SIGTERM or outlived the tool.
    _signal_group(proc, signal.SIGKILL)
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass


def register_live(proc: subprocess.Popen[Any]) -> None:
    """Track a process started outside :func:`run` (a version query), so
    :func:`stop_running_tools` reaches its process group too."""
    with _live_lock:
        _live.add(proc)


def unregister_live(proc: subprocess.Popen[Any]) -> None:
    with _live_lock:
        _live.discard(proc)


def stop_running_tools(sig: int = signal.SIGTERM) -> int:
    """Send ``sig`` to every tool started by :func:`run` that is still running.

    Each tool is signalled with its process group, so helpers it started are
    reached too. Returns the number of tool groups signalled.
    """
    with _live_lock:
        running = list(_live)
    return sum(1 for proc in running if _signal_group(proc, sig))


def install_termination_handler() -> None:
    """Stop running tools when repgenr receives SIGTERM, SIGHUP or SIGINT.

    Python's default action for SIGTERM and SIGHUP ends the interpreter at
    once, leaving the external tool running without its parent (seen live:
    FastTree kept running after its phylo run was terminated). Because each
    tool runs in its own session, it does not receive the terminal's Ctrl-C
    either. The handler sends SIGTERM to every running tool group, sends
    SIGKILL to the groups still present after :data:`STOP_GRACE_SECONDS`, and
    then exits with 128 + the signal number (SIGTERM, SIGHUP) or raises
    :class:`KeyboardInterrupt` (SIGINT), so the stage record stays marked as
    interrupted and a temporary output is removed. A second signal kills the
    remaining tool groups at once, calls the hooks given to
    :func:`register_final_stop_hook`, and takes the default action. Only the main
    thread can install handlers. A signal that repgenr inherited as ignored
    (``nohup`` ignores SIGHUP; a shell ignores SIGINT for a background job of
    a script) stays ignored.

    Ctrl-Z (SIGTSTP) suspends the running tool groups together with repgenr,
    and they resume when repgenr is continued. This applies only when repgenr
    runs as the foreground job of a terminal; elsewhere SIGTSTP keeps its
    default action, which the kernel discards for a process without a
    terminal.
    """
    if threading.current_thread() is not threading.main_thread():
        return

    def _on_signal(signum: int, _frame: object) -> None:
        if stop_requested.is_set():
            # Second signal: do not wait any longer.
            stop_running_tools(signal.SIGKILL)
            _run_final_stop_hooks()
            signal.signal(signum, signal.SIG_DFL)
            os.kill(os.getpid(), signum)
            return
        stop_requested.set()
        stopped = stop_running_tools(signal.SIGTERM)
        # Tools on pool threads are not waited on by the main thread; a timer
        # kills the groups that ignore SIGTERM. Daemon, so it never delays exit.
        escalate = threading.Timer(STOP_GRACE_SECONDS, stop_running_tools, args=(signal.SIGKILL,))
        escalate.daemon = True
        escalate.start()
        # os.write is safe in a signal handler; logging is not.
        os.write(
            2,
            f"repgenr: received signal {signum}; stopped {stopped} running tool(s)\n".encode(),
        )
        if signum == signal.SIGINT:
            raise KeyboardInterrupt
        raise SystemExit(128 + signum)

    def _on_suspend(_signum: int, _frame: object) -> None:
        # Ctrl-Z: suspend the tool groups with repgenr, resume them with it.
        # Only as the foreground job of a terminal, where a shell will send
        # SIGCONT; otherwise the signal is ignored, as the kernel would do.
        if not _foreground_of_terminal():
            return
        # repgenr stops itself with SIGSTOP: a re-sent SIGTSTP is discarded in
        # an orphaned process group, which would leave the tools stopped while
        # repgenr waits for them.
        stop_running_tools(signal.SIGSTOP)
        os.kill(os.getpid(), signal.SIGSTOP)
        # Execution continues here once repgenr receives SIGCONT (fg, bg).
        stop_running_tools(signal.SIGCONT)

    for signum in (signal.SIGTERM, signal.SIGHUP, signal.SIGINT):
        if signal.getsignal(signum) != signal.SIG_IGN:
            signal.signal(signum, _on_signal)
    if (
        hasattr(signal, "SIGTSTP")
        and signal.getsignal(signal.SIGTSTP) == signal.SIG_DFL
        and _terminal_fd() is not None
    ):
        signal.signal(signal.SIGTSTP, _on_suspend)


def _terminal_fd() -> int | None:
    """stdin or stderr when it is a terminal, else None."""
    for fd in (0, 2):
        try:
            if os.isatty(fd):
                return fd
        except OSError:
            continue
    return None


def _foreground_of_terminal() -> bool:
    """True when repgenr's process group is the foreground job of its terminal."""
    fd = _terminal_fd()
    if fd is None:
        return False
    try:
        return os.getpgrp() == os.tcgetpgrp(fd)
    except OSError:  # no controlling terminal
        return False


def default_timeout() -> float | None:
    """The global subprocess time limit, for a caller that sets a limit of
    its own and must not exceed the operator's."""
    return _default_timeout()


def _default_timeout() -> float | None:
    """Global subprocess timeout (seconds) from ``REPGENR_SUBPROCESS_TIMEOUT``.

    Unset (the default) means no timeout, preserving prior behavior; operators
    can cap every external tool with one environment variable. A value that is
    not a positive number is ignored (no timeout) with a single warning.
    """
    global _warned_bad_timeout
    raw = os.environ.get("REPGENR_SUBPROCESS_TIMEOUT")
    if not raw:
        return None
    try:
        value = float(raw)
    except ValueError:
        value = -1.0
    if value > 0:
        return value
    if not _warned_bad_timeout:
        _warned_bad_timeout = True
        _module_logger.warning(
            "Ignoring REPGENR_SUBPROCESS_TIMEOUT=%r: expected a positive "
            "number of seconds; running without a subprocess timeout.",
            raw,
        )
    return None


# Console budget for a tool command line with no prefix; the run log keeps the
# full line. log_command subtracts the prefix so the whole line fits.
CONSOLE_COMMAND_CHARS = console_room("INFO") - len("$ ")
# The console form never gets less room than this, so a long prefix cannot
# reduce it to nothing; such a line runs past CONSOLE_COLUMNS instead.
_MIN_COMMAND_ROOM = 40
# Longest adapter prefix shown on the console, brackets included.
_MAX_CONSOLE_PREFIX = 28
# A run of at least this many consecutive path arguments is counted, not listed.
_MIN_PATH_RUN = 3
_OPTION = re.compile(r"-{1,2}[A-Za-z]")

_Segment = tuple[str, str | list[str]]


def _path_like(arg: str) -> bool:
    return not arg.startswith("-") and os.sep in arg


def _abbreviate(arg: str) -> str:
    name = Path(arg).name
    return f".../{name}" if name and name != arg else arg


def _segments(parts: Sequence[str]) -> list[_Segment]:
    """Split argv into ("word", text), ("value", path), ("path", text) and ("run", paths).

    A path directly after an option (``--reference /r/ref.fasta``) is that
    option's value; a run is three or more consecutive other paths.
    """
    segments: list[_Segment] = []
    i = 0
    while i < len(parts):
        arg = parts[i]
        if not _path_like(arg):
            segments.append(("word", arg))
            i += 1
            continue
        if i > 0 and _OPTION.match(parts[i - 1]):
            segments.append(("value", arg))
            i += 1
            continue
        j = i
        while j < len(parts) and _path_like(parts[j]):
            j += 1
        if j - i >= _MIN_PATH_RUN:
            segments.append(("run", list(parts[i:j])))
        else:
            segments.extend(("path", p) for p in parts[i:j])
        i = j
    return segments


def _render_tokens(segments: list[_Segment], level: int) -> list[str]:
    """One console form of ``segments`` as tokens; a higher level is shorter.

    1: a run is its first path and a count; 2: paths are cut to their last
    component; 3: a run is only its count; 4: a path that is not an option
    value is "...", an option value its last component; 5: every path is "...".
    """
    out: list[str] = []
    for kind, value in segments:
        if kind == "word":
            out.append(str(value))
        elif kind == "run":
            paths = list(value)
            first = _abbreviate(paths[0]) if level >= 2 else paths[0]
            count = f"({len(paths)} paths)"
            out.append(count if level >= 3 else f"{first} ... {count}")
        else:
            path = str(value)
            if level >= 5 or (level >= 4 and kind == "path"):
                out.append("...")
            elif level >= 4:
                out.append(Path(path).name or path)
            elif level >= 2:
                out.append(_abbreviate(path))
            else:
                out.append(path)
    return out


def _head_size(segments: list[_Segment]) -> int:
    """Leading tokens always shown: the program and the token after it."""
    if len(segments) > 1 and segments[1][0] == "word":
        return 2
    return min(len(segments), 1)


def shorten_command(
    command: Sequence[str | os.PathLike[str]], limit: int = CONSOLE_COMMAND_CHARS
) -> str:
    """A console form of ``command``, at most ``limit`` characters where possible.

    A command that fits is returned as it is. Otherwise it is shortened in
    steps until it fits: a run of input paths becomes its first path and a
    count, paths are cut to their last component, a run becomes only its
    count, other paths become "..." (an option's value keeps its last
    component), and then option values too. A line still too long ends at
    the last whole token that fits, followed by "...". The program and the
    token after it are always kept, even past ``limit``.
    """
    parts = [str(p) for p in command]
    full = " ".join(parts)
    if len(full) <= limit:
        return full
    segments = _segments(parts)
    tokens: list[str] = []
    for level in range(1, 6):
        tokens = _render_tokens(segments, level)
        if len(" ".join(tokens)) <= limit:
            return " ".join(tokens)
    keep = _head_size(segments)
    while len(tokens) > keep and len(" ".join([*tokens, "..."])) > limit:
        tokens.pop()
    return " ".join([*tokens, "..."])


def _console_prefix(prefix: str) -> str:
    if len(prefix) <= _MAX_CONSOLE_PREFIX:
        return prefix
    if prefix.startswith("[") and prefix.rstrip().endswith("]"):
        return prefix[: _MAX_CONSOLE_PREFIX - 6] + "...] "
    return prefix[: _MAX_CONSOLE_PREFIX - 4] + "... "


def log_command(
    logger: logging.Logger,
    command: Sequence[str | os.PathLike[str]],
    prefix: str = "",
    *,
    console_command: Sequence[str | os.PathLike[str]] | None = None,
) -> None:
    """Log a tool command line at INFO: in full in the run log, shortened on the console.

    The console form is sized so that the whole line, timestamp, level and
    run-log pointer included, fits ``core.logging.CONSOLE_COLUMNS``; a long
    prefix is capped, and the form keeps at least the program and its first
    option. ``console_command`` is what the console shortens instead of
    ``command``: a container run passes the tool's own argv, without the
    engine preamble.
    """
    full = " ".join(str(p) for p in command)
    shown = command if console_command is None else console_command
    console_prefix = _console_prefix(prefix)
    room = max(console_room("INFO") - len(console_prefix) - len("$ "), _MIN_COMMAND_ROOM)
    short = shorten_command(shown, limit=room)
    if short == full and console_prefix == prefix:
        logger.info("%s$ %s", prefix, full)
        return
    logger.info("%s$ %s", prefix, full, extra=console_extra(f"{console_prefix}$ {short}"))


def run(
    command: Sequence[str | os.PathLike[str]],
    *,
    logger: logging.Logger,
    cwd: str | os.PathLike[str] | None = None,
    env: Mapping[str, str] | None = None,
    check: bool = True,
    stdout_path: str | os.PathLike[str] | None = None,
    log_prefix: str | None = None,
    timeout: float | None = None,
    console_command: Sequence[str | os.PathLike[str]] | None = None,
) -> int:
    """Run ``command`` (an argument vector) without a shell.

    The command line is logged at INFO (in full in the run log, shortened on the
    console; see :func:`log_command`, whose ``console_command`` this passes
    on); the tool's own output is line-streamed
    at DEBUG (progress bars and per-file chatter would otherwise dominate the
    log), and the last lines are kept for the error message on failure. A line
    a tool redraws in place with carriage returns is logged once, in its final
    state. When ``stdout_path`` is given, stdout is written there instead
    (stderr still goes to the logger) -- use this for tools that emit their
    result on stdout (e.g. ``mashtree``).

    The tool starts in its own session and process group, so stopping it
    (timeout, signal, or a failure here) also stops the helpers it started.

    ``timeout`` (seconds) caps the run: on expiry the whole process group is
    killed and :class:`ToolExecutionError` is raised, so a hung tool (or a stuck
    ``docker pull``) cannot wedge the pipeline. It defaults to the
    ``REPGENR_SUBPROCESS_TIMEOUT`` environment variable (unset = no timeout).

    Returns the process exit code. Raises :class:`ToolExecutionError` on a
    non-zero exit when ``check`` is True.
    """
    cmd = [str(part) for part in command]
    prefix = f"[{log_prefix}] " if log_prefix else ""
    if stop_requested.is_set():
        raise ToolExecutionError(cmd, -signal.SIGTERM, output=NOT_STARTED, tool=log_prefix)
    log_command(logger, cmd, prefix, console_command=console_command)

    full_env = {**os.environ, **env} if env else None
    # Refuse a command the operating system would not start; the container
    # backends add their own tokens, so the final argv is measured here.
    check_argv_size(cmd, full_env, tool=log_prefix)
    tail: deque[str] = deque(maxlen=_DEFAULT_TAIL)
    limit = timeout if timeout is not None else _default_timeout()

    # stdout goes to a temp sibling first and is published only when the tool
    # did not fail: several adapters point stdout_path at a final deliverable
    # (e.g. tree/tree.nwk), and opening that in "w" mode would truncate the
    # previous good file the moment a doomed rebuild starts. os.devnull is
    # written directly (it has no meaningful siblings).
    out_target: Path | None = None
    out_tmp: Path | None = None
    if stdout_path is not None:
        if str(stdout_path) == os.devnull:
            out_tmp = Path(os.devnull)
        else:
            out_target = Path(stdout_path)
            out_tmp = out_target.with_name(out_target.name + ".part")
    out_handle = open(out_tmp, "w", encoding="utf-8") if out_tmp is not None else None
    timer: threading.Timer | None = None
    timed_out = False
    proc: subprocess.Popen[bytes] | None = None
    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(cwd) if cwd is not None else None,
            env=full_env,
            # No tool is fed through stdin. A tool that asks a question
            # (skDER below 80 percent ANI) would otherwise wait on the
            # terminal with its prompt hidden in the log; with no stdin it
            # reads end-of-file and fails instead of hanging.
            stdin=subprocess.DEVNULL,
            stdout=(out_handle if out_handle is not None else subprocess.PIPE),
            stderr=subprocess.STDOUT if out_handle is None else subprocess.PIPE,
            # Binary so that a "\r" inside a line is ours to interpret; text
            # mode's universal newlines would split every progress-bar redraw
            # into its own line.
            # Own session and process group: a timeout or a signal stops the
            # tool together with the helpers it starts.
            start_new_session=_OWN_GROUP,
        )
        with _live_lock:
            _live.add(proc)
        if stop_requested.is_set():
            # Started while the handler was stopping the others.
            _signal_group(proc, signal.SIGTERM)

        if limit is not None:

            def _kill() -> None:
                nonlocal timed_out
                timed_out = True
                _signal_group(proc, signal.SIGKILL)

            timer = threading.Timer(limit, _kill)
            timer.start()

        stream = proc.stdout if out_handle is None else proc.stderr
        if stream is not None:
            for raw in stream:
                # A progress bar redraws one line with "\r"; keep its last frame.
                text = raw.decode("utf-8", errors="replace")
                line = text.rstrip("\r\n").rsplit("\r", 1)[-1]
                if not line:
                    continue
                tail.append(line)
                logger.debug("%s%s", prefix, line)
        returncode = proc.wait()
        if stop_requested.is_set() or returncode != 0 or timed_out:
            # The tool ended abnormally (or repgenr is stopping): a helper that
            # outlived it, or ignored SIGTERM, would be left behind once the
            # tool is reaped and leaves the registry.
            _signal_group(proc, signal.SIGKILL)
    except BaseException:
        # A failure in this function (a logging handler, a KeyboardInterrupt,
        # the SystemExit of the termination handler) must not orphan the tool
        # or its helpers: stop the whole group before re-raising.
        if proc is not None:
            stop_group(proc)
        if out_target is not None and out_tmp is not None:
            out_tmp.unlink(missing_ok=True)
        raise
    finally:
        if proc is not None:
            with _live_lock:
                _live.discard(proc)
        if timer is not None:
            timer.cancel()
        if out_handle is not None:
            out_handle.close()

    failed = timed_out or (check and returncode != 0)
    if out_target is not None and out_tmp is not None:
        if failed:
            # Discard the partial capture; a previous good file stays intact.
            out_tmp.unlink(missing_ok=True)
        else:
            os.replace(out_tmp, out_target)
    if timed_out:
        tail.append(f"[killed after {limit}s timeout]")
        raise ToolExecutionError(
            cmd, returncode, output="\n".join(tail), tool=log_prefix, timeout=limit
        )
    if check and returncode != 0:
        raise ToolExecutionError(cmd, returncode, output="\n".join(tail), tool=log_prefix)
    return returncode


def unzip(zip_path: str | os.PathLike[str], dest: str | os.PathLike[str]) -> None:
    """Extract a zip, turning a truncated/corrupt archive into a clear error.

    A `datasets` download that was cut short (network reset, full disk) leaves a
    bad zip; ``zipfile`` then raises ``BadZipFile`` which would surface as a raw
    traceback. Map it to :class:`WorkdirError` naming the file so a re-run retries.
    """
    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(dest)
    except zipfile.BadZipFile as exc:
        raise WorkdirError(
            f"Corrupt or truncated download: {zip_path} ({exc}). Re-run to retry."
        ) from exc


# Used when the platform does not report ARG_MAX: the common value on macOS
# and on Linux with the default 8 MiB stack.
_FALLBACK_ARG_MAX = 1 << 20
# Bytes kept free below ARG_MAX: execve also stores the program path and the
# auxiliary vector, and the kernel's accounting differs slightly between
# platforms.
_ARG_MAX_MARGIN = 16 * 1024
# Linux refuses a single argument or environment string longer than this
# (MAX_ARG_STRLEN, 32 pages), whatever ARG_MAX is.
_LINUX_MAX_ARG_STRLEN = 32 * 4096
# warn_argv_bytes warns once a command uses this fraction of the limit.
_ARGV_WARN_FRACTION = 0.75
# What a user can do about an oversized command line, named in the warning
# and in the error.
_ARGV_ADVICE = (
    "Reduce the number of genomes per call (for dereplication, --process-size "
    "bounds the genomes in one tool call), or dereplicate before this stage."
)


def arg_max() -> int:
    """The operating system's limit on argv plus environment, in bytes."""
    try:
        value = os.sysconf("SC_ARG_MAX")
    except (AttributeError, ValueError, OSError):
        return _FALLBACK_ARG_MAX
    return value if value > 0 else _FALLBACK_ARG_MAX


def _strings_bytes(strings: Sequence[str]) -> int:
    """Bytes execve needs for ``strings``: each one with its NUL and pointer."""
    return sum(len(os.fsencode(s)) + 1 + 8 for s in strings)


def _env_strings(env: Mapping[str, str] | None) -> list[str]:
    return [f"{k}={v}" for k, v in (os.environ if env is None else env).items()]


def argv_bytes(argv: Sequence[str | os.PathLike[str]], env: Mapping[str, str] | None = None) -> int:
    """Bytes ``argv`` plus the environment occupy in an execve call.

    ``env`` is the child's full environment (None: this process's).
    """
    return _strings_bytes([os.fspath(a) for a in argv]) + _strings_bytes(_env_strings(env))


def check_argv_size(
    argv: Sequence[str | os.PathLike[str]],
    env: Mapping[str, str] | None = None,
    *,
    tool: str | None = None,
) -> None:
    """Refuse a command that the operating system would not start (E2BIG).

    Measures the argument vector and the environment against ARG_MAX less a
    margin and, on Linux, every string against the 128 KiB single-argument
    limit. Raises :class:`UserInputError` naming --process-size before the
    tool is started, instead of an opaque "Argument list too long" from exec.
    """
    args = [os.fspath(a) for a in argv]
    name = tool or (Path(args[0]).name if args else "tool")
    total = argv_bytes(args, env)
    limit = arg_max() - _ARG_MAX_MARGIN
    if total > limit:
        raise UserInputError(
            f"{name}: the command line and environment take {total // 1024} KiB "
            f"across {len(args)} arguments, more than the operating system "
            f"accepts ({limit // 1024} KiB of ARG_MAX). {_ARGV_ADVICE}"
        )
    if sys.platform.startswith("linux"):
        for s in [*args, *_env_strings(env)]:
            if len(os.fsencode(s)) + 1 > _LINUX_MAX_ARG_STRLEN:
                raise UserInputError(
                    f"{name}: one command-line or environment string is "
                    f"{len(os.fsencode(s)) // 1024} KiB long, more than the "
                    f"{_LINUX_MAX_ARG_STRLEN // 1024} KiB Linux accepts for a "
                    f"single argument ({s[:60]}...). {_ARGV_ADVICE}"
                )


def warn_argv_bytes(
    tool: str, argv: Sequence[str | os.PathLike[str]], logger: logging.Logger
) -> None:
    """Warn when an argv approaches the operating system's ARG_MAX limit.

    For tools that take every genome path on argv and have no list-file
    alternative (SibeliaZ, snippy-core). The size is measured
    (argv plus this process's environment) against the platform's ARG_MAX, not
    estimated from a genome count. :func:`run` refuses a command over the
    limit; this warns earlier, from three quarters of it.
    """
    total = argv_bytes(argv)
    limit = arg_max() - _ARG_MAX_MARGIN
    if total >= _ARGV_WARN_FRACTION * limit:
        logger.warning(
            "%s receives %d KiB across %d command-line arguments (with the "
            "environment), %.0f%% of the operating system's ARG_MAX limit (%d KiB usable). %s",
            tool,
            total // 1024,
            len(argv),
            100.0 * total / limit,
            limit // 1024,
            _ARGV_ADVICE,
        )


def write_fofn(paths: Sequence[str | os.PathLike[str]], dest: str | os.PathLike[str]) -> Path:
    """Write a file-of-filenames (one absolute path per line) and return its path.

    Use this instead of shell globs when handing a large genome set to a tool.

    Paths are absolute but NOT symlink-resolved: the container backend binds
    un-resolved abspaths (macOS firmlinks resolve /Users -> /System/Volumes/Data,
    which is outside Docker's shared directories), so a tool reading this fofn
    inside a container must see the same un-resolved paths.
    """
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(dest_path, "w", encoding="utf-8") as fo:
        for p in paths:
            fo.write(f"{os.path.abspath(os.fspath(p))}\n")
    return dest_path


def link_or_copy(src: str | os.PathLike[str], dst: str | os.PathLike[str]) -> bool:
    """Stage ``src`` at ``dst`` cheaply: hardlink it, copying only as a fallback.

    Staging genomes into representatives/cluster dirs copies tens of GB at 1000s
    of genomes. A hardlink is instant and uses no extra disk; it shares the inode
    with the source, which is safe because these staged files are only read by
    downstream stages, never modified in place. Falls back to a real copy when
    the filesystem can't hardlink (cross-device, or exFAT/NTFS on the dev box).
    Returns True when the file was linked, False when it was copied.
    """
    # Resolve symlinks to the real file first. Tools such as skDER emit their
    # representative genomes as symlinks (often into a Nextflow-/container-staged
    # input tree); hardlinking the symlink itself -- which os.link does on macOS --
    # produces a broken, 0-byte staged file. Linking the real target instead keeps
    # the content and still shares the inode (no extra disk). realpath() walks the
    # whole path with lstat, so only pay it when src is actually a symlink (the
    # common case -- staging plain genome files -- skips it).
    src_path = os.fspath(src)
    src_s = os.path.realpath(src_path) if os.path.islink(src_path) else src_path
    dst_s = os.fspath(dst)
    try:
        os.link(src_s, dst_s)
    except OSError:
        shutil.copy2(src_s, dst_s)
        return False
    return True


def _ignore_vanished(func, path, exc):  # noqa: ANN001
    """rmtree error handler: an entry that disappeared mid-walk is not an error.

    macOS removes a file's AppleDouble twin (``._name``) on volumes without
    native extended attributes the moment the data file goes, so the walk
    then fails to unlink a name that is already gone.
    """
    if isinstance(exc, FileNotFoundError):
        return
    raise exc


def remove_tree(path: str | os.PathLike[str]) -> None:
    """Remove a directory tree, tolerating entries that vanish during the walk."""
    target = Path(path)
    if not target.exists():
        return
    shutil.rmtree(target, onexc=_ignore_vanished)
    if target.exists():  # a second pass catches what the first walk skipped
        shutil.rmtree(target, onexc=_ignore_vanished)


@contextmanager
def staged_dir(final: str | os.PathLike[str]) -> Iterator[Path]:
    """Build a directory's new contents beside it; publish them only when complete.

    Yields a sibling staging directory to fill. On clean exit the previous
    ``final`` (if any) is removed and the staging directory takes its place;
    on an exception the staging directory is removed and ``final`` is left as
    it was. A crash mid-write therefore never leaves a partial set behind.
    """
    final = Path(final)
    staging = final.parent / f".{final.name}.staging"
    if staging.exists():
        remove_tree(staging)
    staging.mkdir(parents=True)
    try:
        yield staging
    except BaseException:
        remove_tree(staging)
        raise
    if final.exists():
        remove_tree(final)
    os.replace(staging, final)


# Hard floor for any stage about to write a lot: below this a download or an
# assembly is refused rather than allowed to fill the volume mid-run.
MIN_FREE_BYTES = 1_000_000_000


def check_free_disk(path: str | os.PathLike[str], estimate: int, logger, *, what: str) -> None:
    """Refuse to start ``what`` with almost no free disk under ``path``; warn when tight.

    ``estimate`` is the caller's rough byte need (a genome count times a
    typical size, the FASTQ sizes an archive reported); it decides the warning
    only, the floor is absolute.
    """
    free = shutil.disk_usage(path).free
    if free < MIN_FREE_BYTES:
        raise WorkdirError(
            f"Only {free / 1e9:.1f} GB free under {path}; refusing to {what}. "
            "Free disk space or point the working directory at a larger volume."
        )
    if free < estimate:
        logger.warning(
            "Low disk: ~%.1f GB free, up to ~%.1f GB may be needed to %s.",
            free / 1e9,
            estimate / 1e9,
            what,
        )


_GZIP_MAGIC = b"\x1f\x8b"


def is_gzip(path: str | os.PathLike[str]) -> bool:
    """True when ``path`` starts with the gzip magic bytes, whatever its name."""
    with open(path, "rb") as fh:
        return fh.read(2) == _GZIP_MAGIC


def copy_plain_fasta(src: str | os.PathLike[str], dest: str | os.PathLike[str]) -> None:
    """Copy ``src`` to ``dest`` uncompressed (gzip judged by magic bytes).

    A truncated or corrupt gzip file is a UserInputError naming ``src``.
    """
    opener = gzip.open if is_gzip(src) else open
    try:
        with opener(src, "rb") as fi, open(dest, "wb") as fo:
            shutil.copyfileobj(fi, fo, 1 << 20)
    except (EOFError, gzip.BadGzipFile, zlib.error) as exc:
        raise UserInputError(
            f"Genome {src} is a truncated or corrupt gzip file ({exc}); "
            "replace it with a complete copy."
        ) from exc


def stage_plain_inputs(
    paths: Sequence[Path],
    caps: ToolCapabilities,
    dest_dir: Path,
    logger: logging.Logger,
) -> dict[Path, Path]:
    """Map each input genome to a path the tool of ``caps`` can read.

    A tool that reads gzip (``caps.reads_gzip``) gets every path unchanged.
    Otherwise each gzipped input (by magic bytes, not by name) is decompressed
    to ``dest_dir/<record name>.fasta``, written through a temporary sibling,
    and plain inputs map to themselves. The copy keeps the record name, so
    alignments and trees name the genome as they would the original. Two
    inputs with one record name (``x.fasta`` and ``x.fasta.gz``) are refused
    for every tool, since their records and tree leaves could not be told
    apart. ``dest_dir`` is created only when a copy is written; the caller
    removes it when the tool has finished.
    """
    unique = list(dict.fromkeys(paths))
    by_name: dict[str, Path] = {}
    for p in unique:
        other = by_name.setdefault(record_name(p), p)
        if other != p:
            raise UserInputError(
                f"Two input genomes share the record name '{record_name(p)}': {other} "
                f"and {p}. Remove or rename one of them."
            )
    if caps.reads_gzip:
        return {p: p for p in unique}
    gzipped = [p for p in unique if is_gzip(p)]
    if gzipped:
        # Genome FASTA compresses about three- to fourfold; refuse before the
        # first copy rather than fail part-way through.
        check_free_disk(
            _existing_parent(dest_dir),
            4 * sum(p.stat().st_size for p in gzipped),
            logger,
            what=f"decompress {len(gzipped)} genome(s) for {caps.name}",
        )
    staged: dict[Path, Path] = {}
    for p in unique:
        if p not in gzipped:
            staged[p] = p
            continue
        dest = dest_dir / f"{record_name(p)}.fasta"
        with atomic_path(dest) as tmp:
            copy_plain_fasta(p, tmp)
        staged[p] = dest
    copies = sum(1 for p, d in staged.items() if p != d)
    if copies:
        logger.info(
            "%s does not read gzipped FASTA; decompressed %d genome(s) to %s",
            caps.name,
            copies,
            dest_dir,
        )
    return staged


def _existing_parent(path: Path) -> Path:
    """``path`` or its nearest existing ancestor, for a free-space query."""
    for candidate in (path, *path.parents):
        if candidate.exists():
            return candidate
    return Path(".")


@dataclass
class PipelineResult:
    """The outcome of :func:`run_pipeline`: exit codes and the last stderr lines
    of the producer and the consumer."""

    returncodes: tuple[int, int]
    tails: tuple[list[str], list[str]]
    timed_out: bool = False

    @property
    def ok(self) -> bool:
        return self.returncodes == (0, 0) and not self.timed_out


def run_pipeline(
    producer: Sequence[str | os.PathLike[str]],
    consumer: Sequence[str | os.PathLike[str]],
    *,
    logger: logging.Logger,
    cwd: str | os.PathLike[str] | None = None,
    envs: tuple[Mapping[str, str] | None, Mapping[str, str] | None] = (None, None),
    log_prefixes: tuple[str, str] = ("producer", "consumer"),
    timeout: float | None = None,
    console_commands: tuple[Sequence[str], Sequence[str]] | None = None,
) -> PipelineResult:
    """Run ``producer | consumer`` without a shell and return both exit codes.

    ``envs`` holds the environment additions of each half (a --bin-dir PATH
    applies to its own tool only).

    The producer's stdout is the consumer's stdin; the stderr of both is
    line-streamed at DEBUG under its prefix and the last lines of each are
    returned for parsing and for error messages. Each process leads its own
    session and process group, is registered as a live tool (a termination
    signal stops both) and is killed when the other ends abnormally, when the
    stop request is set, or when ``timeout`` expires. A consumer that exits
    early leaves the producer with a broken pipe, reported as the producer's
    exit code; the caller decides what that means.

    Nothing is raised for a non-zero exit: a pipeline over many runs treats a
    failure as one failed run.
    """
    first = [str(part) for part in producer]
    second = [str(part) for part in consumer]
    if stop_requested.is_set():
        raise ToolExecutionError(first, -signal.SIGTERM, output=NOT_STARTED, tool=log_prefixes[0])
    shown = console_commands or (first, second)
    log_command(logger, first, f"[{log_prefixes[0]}] ", console_command=shown[0])
    log_command(logger, second, f"[{log_prefixes[1]}] ", console_command=shown[1])
    full_envs = tuple({**os.environ, **env} if env else None for env in envs)
    check_argv_size(first, full_envs[0], tool=log_prefixes[0])
    check_argv_size(second, full_envs[1], tool=log_prefixes[1])
    limit = timeout if timeout is not None else _default_timeout()
    tails: tuple[deque[str], deque[str]] = (
        deque(maxlen=_DEFAULT_TAIL),
        deque(maxlen=_DEFAULT_TAIL),
    )
    procs: list[subprocess.Popen[bytes]] = []
    timer: threading.Timer | None = None
    timed_out = False

    def drain(stream: Any, prefix: str, tail: deque[str]) -> None:
        # The pipe is read to its end whatever happens to a line: a reader
        # that stopped would let the pipe fill and the tool block on it.
        for raw in stream:
            try:
                text = raw.decode("utf-8", errors="replace")
                line = text.rstrip("\r\n").rsplit("\r", 1)[-1]
                if line:
                    tail.append(line)
                    logger.debug("[%s] %s", prefix, line)
            except Exception:  # noqa: BLE001 - a logging failure must not stall the pipe
                continue

    head: subprocess.Popen[bytes] | None = None
    try:
        head = subprocess.Popen(
            first,
            cwd=str(cwd) if cwd is not None else None,
            env=full_envs[0],
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=_OWN_GROUP,
        )
        procs.append(head)
        tail_proc = subprocess.Popen(
            second,
            cwd=str(cwd) if cwd is not None else None,
            env=full_envs[1],
            stdin=head.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            start_new_session=_OWN_GROUP,
        )
        procs.append(tail_proc)
        # This is the parent's copy of the pipe's read end; the consumer
        # holds its own. Closing it here lets the producer see the broken
        # pipe (EPIPE) once the consumer ends.
        assert head.stdout is not None
        head.stdout.close()
        with _live_lock:
            _live.update(procs)
        if stop_requested.is_set():
            for proc in procs:
                _signal_group(proc, signal.SIGTERM)
        if limit is not None:

            def _kill() -> None:
                nonlocal timed_out
                timed_out = True
                for proc in procs:
                    _signal_group(proc, signal.SIGKILL)

            timer = threading.Timer(limit, _kill)
            timer.start()
        readers = [
            threading.Thread(
                target=drain, args=(head.stderr, log_prefixes[0], tails[0]), daemon=True
            ),
            threading.Thread(
                target=drain, args=(tail_proc.stderr, log_prefixes[1], tails[1]), daemon=True
            ),
            threading.Thread(
                target=drain, args=(tail_proc.stdout, log_prefixes[1], tails[1]), daemon=True
            ),
        ]
        for reader in readers:
            reader.start()
        # The consumer ends when its input ends or when it has seen enough;
        # a producer still writing then gets a broken pipe and ends too.
        tail_rc = tail_proc.wait()
        if tail_rc != 0:
            # SIGTERM first: a container client forwards it to its tool,
            # which SIGKILL would not reach. stop_group escalates.
            stop_group(head)
        head_rc = head.wait()
        for reader in readers:
            reader.join()
        if stop_requested.is_set() or head_rc != 0 or tail_rc != 0 or timed_out:
            for proc in procs:
                _signal_group(proc, signal.SIGKILL)
    except BaseException:
        for proc in procs:
            stop_group(proc)
        if head is not None and len(procs) == 1:
            # The consumer never started: the producer's pipes are ours to close.
            for stream in (head.stdout, head.stderr):
                if stream is not None:
                    stream.close()
        raise
    finally:
        with _live_lock:
            for proc in procs:
                _live.discard(proc)
        if timer is not None:
            timer.cancel()
    if timed_out:
        tails[1].append(f"[killed after {limit}s timeout]")
    return PipelineResult((head_rc, tail_rc), (list(tails[0]), list(tails[1])), timed_out=timed_out)
