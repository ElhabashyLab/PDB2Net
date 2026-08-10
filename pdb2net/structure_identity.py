"""Stable identities for official PDB entries and user-supplied structures.

The wwPDB extended identifier is the canonical representation for official
entries.  Legacy four-character identifiers remain available as display and
reference-data aliases while they exist.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Mapping


LEGACY_PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
EXTENDED_PDB_ID_RE = re.compile(r"^pdb_[A-Za-z0-9]{8}$", re.IGNORECASE)


def extended_pdb_id(legacy_id: str) -> str:
    """Convert a legacy four-character ID to its wwPDB extended alias."""
    value = str(legacy_id or "").strip()
    if not LEGACY_PDB_ID_RE.fullmatch(value):
        raise ValueError(f"Invalid legacy PDB ID: {legacy_id!r}")
    return f"pdb_0000{value.lower()}"


def legacy_pdb_id(extended_id: str) -> str | None:
    """Return the legacy alias encoded by an extended ID, when one exists."""
    value = str(extended_id or "").strip().lower()
    if not EXTENDED_PDB_ID_RE.fullmatch(value):
        return None
    payload = value[4:]
    candidate = payload[-4:]
    if payload[:4] == "0000" and LEGACY_PDB_ID_RE.fullmatch(candidate):
        return candidate.lower()
    return None


def canonical_pdb_id(value: str) -> str:
    """Normalize either official PDB ID form to the extended representation."""
    raw = str(value or "").strip()
    if LEGACY_PDB_ID_RE.fullmatch(raw):
        return extended_pdb_id(raw)
    if EXTENDED_PDB_ID_RE.fullmatch(raw):
        return raw.lower()
    raise ValueError(f"Invalid PDB ID: {value!r}")


def pdb_archive_shard(value: str) -> str:
    """Return the two-character wwPDB hash directory for an official ID."""
    canonical = canonical_pdb_id(value)
    return canonical[-3:-1]


@dataclass(frozen=True)
class StructureIdentity:
    """Minimal stable identity carried across parser/cache boundaries."""

    source: str
    canonical_id: str
    legacy_id: str | None = None

    def __post_init__(self) -> None:
        source = str(self.source or "").strip().lower()
        canonical = str(self.canonical_id or "").strip()
        legacy = str(self.legacy_id).strip().lower() if self.legacy_id else None
        if source not in {"pdb", "local", "alphafold"}:
            raise ValueError("Structure identity source is unsupported.")
        if not canonical:
            raise ValueError("Structure canonical ID must not be empty.")
        if source == "pdb":
            canonical = canonical_pdb_id(canonical)
            encoded_legacy = legacy_pdb_id(canonical)
            if legacy is not None and legacy != encoded_legacy:
                raise ValueError("Legacy and extended PDB IDs do not describe the same entry.")
            legacy = encoded_legacy
        else:
            if (
                len(canonical.encode("utf-8")) > 512
                or "/" in canonical
                or "\\" in canonical
                or any(ord(character) < 32 or ord(character) == 127 for character in canonical)
            ):
                raise ValueError("Local structure identity is unsafe or oversized.")
            if legacy is not None:
                raise ValueError("Local structure identities do not have legacy PDB aliases.")
        object.__setattr__(self, "source", source)
        object.__setattr__(self, "canonical_id", canonical)
        object.__setattr__(self, "legacy_id", legacy)

    @property
    def key(self) -> str:
        """Stable key suitable for caches, jobs, and cross-file joins."""
        return self.canonical_id if self.source == "pdb" else f"{self.source}:{self.canonical_id}"

    @property
    def display_id(self) -> str:
        """Compact human-readable label without weakening the canonical key."""
        if self.legacy_id:
            return self.legacy_id.upper()
        return self.canonical_id

    def as_dict(self) -> dict[str, str | None]:
        return {
            "source": self.source,
            "canonical_id": self.canonical_id,
            "legacy_id": self.legacy_id,
            "key": self.key,
            "display_id": self.display_id,
        }

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "StructureIdentity":
        return cls(
            source=str(value.get("source") or "local"),
            canonical_id=str(value.get("canonical_id") or value.get("key") or ""),
            legacy_id=str(value["legacy_id"]) if value.get("legacy_id") else None,
        )


def identity_from_official_id(value: str) -> StructureIdentity:
    canonical = canonical_pdb_id(value)
    return StructureIdentity("pdb", canonical, legacy_pdb_id(canonical))


@dataclass(frozen=True)
class ChainIdentity:
    """Structured chain identity used internally instead of parsing node labels."""

    structure_key: str
    structure_display_id: str
    chain_id: str
    model_index: int = 1

    def __post_init__(self) -> None:
        structure_key = str(self.structure_key or "").strip()
        display = str(self.structure_display_id or "").strip()
        chain = str(self.chain_id or "").strip()
        if not structure_key or not display or not chain:
            raise ValueError("Chain identity fields must not be empty.")
        if isinstance(self.model_index, bool) or not isinstance(self.model_index, int):
            raise ValueError("Chain model index must be a positive integer.")
        if self.model_index < 1:
            raise ValueError("Chain model index must be a positive integer.")
        object.__setattr__(self, "structure_key", structure_key)
        object.__setattr__(self, "structure_display_id", display)
        object.__setattr__(self, "chain_id", chain)

    def node_id(self, *, include_model: bool) -> str:
        """Return the backwards-compatible public node label."""
        if include_model:
            return f"{self.structure_display_id}:model{self.model_index}:{self.chain_id}"
        return f"{self.structure_display_id}:{self.chain_id}"

    def as_dict(self) -> dict[str, Any]:
        return {
            "structure_key": self.structure_key,
            "structure_display_id": self.structure_display_id,
            "chain_id": self.chain_id,
            "model_index": self.model_index,
        }

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "ChainIdentity":
        return cls(
            structure_key=str(value.get("structure_key") or ""),
            structure_display_id=str(value.get("structure_display_id") or ""),
            chain_id=str(value.get("chain_id") or ""),
            model_index=value.get("model_index", 1),
        )


def chain_identity_from_record(chain: Mapping[str, Any]) -> ChainIdentity:
    """Read the canonical structured identity carried on a chain record."""
    raw = chain.get("chain_identity")
    if isinstance(raw, Mapping):
        return ChainIdentity.from_mapping(raw)
    return ChainIdentity(
        structure_key=str(chain.get("structure_key") or chain.get("_parent_structure_key") or ""),
        structure_display_id=str(
            chain.get("structure_display_id") or chain.get("_parent_pdb_id") or ""
        ),
        chain_id=str(chain.get("chain_id") or ""),
        model_index=chain.get("model_index", 1),
    )


def edge_structure_key(edge: Mapping[str, Any]) -> str:
    """Return the structure key explicitly attached to a scientific edge."""
    return str(edge.get("structure_key") or "")


__all__ = [
    "EXTENDED_PDB_ID_RE",
    "ChainIdentity",
    "LEGACY_PDB_ID_RE",
    "StructureIdentity",
    "canonical_pdb_id",
    "chain_identity_from_record",
    "edge_structure_key",
    "extended_pdb_id",
    "identity_from_official_id",
    "legacy_pdb_id",
    "pdb_archive_shard",
]
