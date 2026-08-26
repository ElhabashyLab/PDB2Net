"""Stable machine-readable interface versions shared across PDB2Net modules."""

from __future__ import annotations


CAPABILITIES_SCHEMA_VERSION = "3"
CLI_CONTRACT_VERSION = "2"
OUTPUT_CONTRACT_VERSION = "2.0"
PRECOMPUTED_ARTIFACT_SCHEMA_VERSION = "3"


__all__ = [
    "CAPABILITIES_SCHEMA_VERSION",
    "CLI_CONTRACT_VERSION",
    "OUTPUT_CONTRACT_VERSION",
    "PRECOMPUTED_ARTIFACT_SCHEMA_VERSION",
]
