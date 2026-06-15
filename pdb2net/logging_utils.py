"""Small logging helpers for server-friendly diagnostics."""

from __future__ import annotations

import logging
import os


def get_logger(name: str) -> logging.Logger:
    """Return a project logger with a quiet default handler."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(levelname)s:%(name)s:%(message)s"))
        logger.addHandler(handler)
    logger.propagate = False
    logger.setLevel(logging.DEBUG if _verbose_enabled() else logging.INFO)
    return logger


def _verbose_enabled() -> bool:
    return os.environ.get("PDB2NET_VERBOSE", "").strip().lower() in {"1", "true", "yes", "on"}
