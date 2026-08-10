"""Shared input-contract errors for folder-based processing."""

from __future__ import annotations


class InputValidationError(ValueError):
    """Raised when the configured folder-based input contract is violated."""

    def __init__(self, code: str, message: str) -> None:
        self.code = str(code)
        self.message = str(message)
        super().__init__(f"{self.code}: {self.message}")

    def __reduce__(self):
        """Keep the structured code when an error crosses a process boundary."""
        return self.__class__, (self.code, self.message)
