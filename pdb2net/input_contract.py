"""Shared input-contract errors for folder-based processing."""

from __future__ import annotations


class InputValidationError(ValueError):
    """Raised when the configured folder-based input contract is violated."""

    def __init__(self, code: str, message: str) -> None:
        self.code = code
        super().__init__(f"{code}: {message}")
