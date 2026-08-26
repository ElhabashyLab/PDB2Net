"""Small public interface for offline precompute and read-only assembly."""

from .assemble import load_entry, load_manifest, run_assemble_pipeline
from .build import precompute_directory, precompute_sources
from .io import entry_path
from .schema import (
    normalize_pdb_id,
    normalize_requested_ids,
    profile_id,
    scientific_profile,
)

__all__ = [
    "entry_path",
    "load_entry",
    "load_manifest",
    "normalize_pdb_id",
    "normalize_requested_ids",
    "precompute_directory",
    "precompute_sources",
    "profile_id",
    "run_assemble_pipeline",
    "scientific_profile",
]
