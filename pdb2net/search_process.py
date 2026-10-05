"""Local cancellation of direct native search processes, without new sessions."""

from __future__ import annotations

import subprocess
import threading
import time

TERMINATION_GRACE_SECONDS = 5.0
_POLL_SECONDS = 0.1


class SearchAborted(RuntimeError):
    """A timeout or another cancelled search stopped this search group."""


class SearchProcesses:
    """One search invocation's cancellation state; each worker owns its Popen."""

    def __init__(self, timeout_seconds: float | None):
        self.timeout_seconds = timeout_seconds
        self._lock = threading.Lock()
        self._reason: str | None = None
        self._cancelled_at: float | None = None

    def check(self) -> None:
        with self._lock:
            if self._reason is not None:
                raise SearchAborted(self._reason)

    def cancel(self, reason: str = "Sequence search interrupted.") -> None:
        with self._lock:
            if self._reason is None:
                self._reason = reason
                self._cancelled_at = time.monotonic()

    def _stop(self, process: subprocess.Popen) -> None:
        if process.poll() is None:
            try:
                process.terminate()
            except ProcessLookupError:
                pass
        # All concurrent searches share the same grace deadline, not five
        # seconds per child. communicate drains pipes before temporary cleanup.
        with self._lock:
            deadline = self._cancelled_at + TERMINATION_GRACE_SECONDS
        try:
            process.communicate(timeout=max(0.0, deadline - time.monotonic()))
        except subprocess.TimeoutExpired:
            process.kill()
            process.communicate()

    def run(self, command: list[str], *, tool: str) -> subprocess.CompletedProcess:
        # Cancellation cannot race a launch: a launch already holding this lock
        # will immediately see cancellation in its polling loop after it starts.
        process = None
        try:
            with self._lock:
                if self._reason is not None:
                    raise SearchAborted(self._reason)
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            started = time.monotonic()
            while True:
                self.check()
                interval = _POLL_SECONDS
                if self.timeout_seconds is not None:
                    remaining = self.timeout_seconds - (time.monotonic() - started)
                    if remaining <= 0:
                        self.cancel(
                            f"{tool} exceeded the configured time limit of {self.timeout_seconds:g} seconds; "
                            "sequence search aborted."
                        )
                        self.check()
                    interval = min(interval, remaining)
                try:
                    stdout, stderr = process.communicate(timeout=interval)
                except subprocess.TimeoutExpired:
                    continue
                self.check()
                return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
        except BaseException as exc:
            if process is not None or isinstance(exc, KeyboardInterrupt):
                self.cancel()
                if process is not None:
                    self._stop(process)
            raise
        finally:
            if process is not None:
                if process.stdout is not None:
                    process.stdout.close()
                if process.stderr is not None:
                    process.stderr.close()
