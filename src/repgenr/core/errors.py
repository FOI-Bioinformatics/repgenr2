"""Exception hierarchy for RepGenR2.

Stages and plugins raise these instead of calling ``sys.exit``. The CLI layer
(:mod:`repgenr.cli.main`) catches :class:`RepGenRError`, logs it cleanly and
sets a non-zero exit code, so library callers and tests can handle failures
without a process exit.
"""

from __future__ import annotations

from pathlib import Path

# `repgenr doctor` exits with this status when at least one finding is a
# failure, so a caller can tell a damaged workdir from a crashed check (1).
DOCTOR_FAILURES_EXIT = 7


class RepGenRError(Exception):
    """Base class for all expected, user-facing RepGenR2 failures.

    ``exit_code`` is the process status the CLI exits with for the class, so
    a caller can tell a usage error from a missing tool or a failed run
    without parsing the log. 1 is reserved for unexpected exceptions.
    """

    exit_code: int = 1


class UserInputError(RepGenRError):
    """Invalid or missing user input (bad arguments, missing target, etc.)."""

    exit_code = 2


class WorkdirError(RepGenRError):
    """A working directory is missing required files or is in a bad state."""

    exit_code = 3


class MissingBinaryError(RepGenRError):
    """A required external tool was not found on PATH or failed its version check."""

    exit_code = 4


class ToolExecutionError(RepGenRError):
    """An external tool exited with a non-zero status.

    The CLI exits 6, or with the tool's own status under
    ``REPGENR_PROPAGATE_TOOL_EXIT=1`` (see ``cli.base``).
    """

    exit_code = 6

    def __init__(
        self,
        command: list[str],
        returncode: int,
        output: str | None = None,
        *,
        tool: str | None = None,
        timeout: float | None = None,
        message: str | None = None,
    ):
        self.command = command
        self.returncode = returncode
        self.output = output
        self.timeout = timeout
        # The adapter's tool name when known: a containerized command starts
        # with the engine, which says nothing about what failed.
        self.tool = tool or (Path(command[0]).name if command else "tool")
        if timeout is not None:
            msg = f"{self.tool} timed out after {timeout:g}s (killed, exit {returncode})"
        elif returncode == 0:
            # Some tools report a fatal problem and still exit 0 (dRep without
            # CheckM); the adapter raises with returncode 0 when the expected
            # output is absent.
            msg = f"{self.tool} exited 0 without writing its results"
        else:
            msg = f"{self.tool} failed (exit {returncode})"
        super().__init__(message or msg)

    @property
    def output_tail(self) -> str:
        return self.output or ""

    def details(self) -> str:
        """The command line and the output tail, for the run log."""
        text = "command: " + " ".join(self.command)
        if self.output:
            text += f"\n--- output tail ---\n{self.output}"
        return text


class PluginError(RepGenRError):
    """A requested plugin/tool could not be found or loaded."""

    exit_code = 5
