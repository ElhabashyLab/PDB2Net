"""Configuration-free description of PDB2Net's server-facing interfaces."""

from __future__ import annotations

from typing import Any

from . import __version__
from .contracts import (
    CAPABILITIES_SCHEMA_VERSION,
    CLI_CONTRACT_VERSION,
    OUTPUT_CONTRACT_VERSION,
    PRECOMPUTED_ARTIFACT_SCHEMA_VERSION,
    WEB_CONFIG_SCHEMA_VERSION,
)
from .server_interface import server_interface_document
from .server_interface import WEB_SERVER_INPUT_SUFFIXES


def capability_document() -> dict[str, Any]:
    """Return stable capabilities without loading machine-local configuration."""
    return {
        "capabilities_schema_version": CAPABILITIES_SCHEMA_VERSION,
        "pdb2net_version": __version__,
        "cli_contract_version": CLI_CONTRACT_VERSION,
        "output_contract_version": OUTPUT_CONTRACT_VERSION,
        "precomputed_artifact_schema_version": PRECOMPUTED_ARTIFACT_SCHEMA_VERSION,
        "web_config_schema_version": WEB_CONFIG_SCHEMA_VERSION,
        "commands": ["run", "precompute", "assemble", "capabilities"],
        "input_formats": list(WEB_SERVER_INPUT_SUFFIXES),
        "network_outputs": [
            "chain_per_pdb",
            "protein_per_pdb",
            "combined_chain_network",
            "combined_protein_network",
        ],
        "network_annotation_databases": ["uniprot", "pfam", "cath", "scop2"],
        "structure_model_policies": ["first", "all"],
        "features": [
            "headless_cx2",
            "detailed_interactions",
            "embedded_sifts",
            "diamond_uniref90",
            "precomputed_store",
        ],
        "server_interface": server_interface_document(),
    }


__all__ = ["capability_document"]
