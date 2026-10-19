"""Stable machine-readable interface versions shared across PDB2Net modules."""

from __future__ import annotations


CAPABILITIES_SCHEMA_VERSION = "4"
CLI_CONTRACT_VERSION = "3"
OUTPUT_CONTRACT_VERSION = "3.0"
STRUCTURE_FACTS_SCHEMA_VERSION = "1"


__all__ = [
    "CAPABILITIES_SCHEMA_VERSION",
    "CLI_CONTRACT_VERSION",
    "OUTPUT_CONTRACT_VERSION",
    "STRUCTURE_FACTS_SCHEMA_VERSION",
]
