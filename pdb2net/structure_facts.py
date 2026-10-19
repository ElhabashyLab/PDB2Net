"""Reference-independent structure facts shared with external input adapters."""
from __future__ import annotations

import copy
from datetime import datetime
import hashlib
import json
import math
import re
from collections import Counter

import gemmi

from .contracts import STRUCTURE_FACTS_SCHEMA_VERSION
from .data_processor import PARSER_SEMANTICS
from .distances import INTERACTION_PIPELINE_VERSION
from .input_contract import InputValidationError
from .network_annotations import validate_embedded_annotations
from .structure_identity import ChainIdentity, StructureIdentity

MAX_CHAINS = 50_000
MAX_PAIRS = 500_000
MAX_REQUEST_CHAINS = 200_000
MAX_REQUEST_PAIRS = 2_000_000
CHAIN_FIELDS = {
    "chain_id", "unique_chain_id", "chain_identity", "structure_key",
    "structure_display_id", "model_index", "structure_sequence",
    "structure_molecule_type", "aa_len", "nt_len", "embedded_annotations",
}


def canonical_json(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")


def geometry_profile():
    """Only changes that affect prepared structural information invalidate it."""
    return {
        "facts_schema": STRUCTURE_FACTS_SCHEMA_VERSION,
        "parser": PARSER_SEMANTICS, "gemmi": gemmi.__version__,
        "contacts": INTERACTION_PIPELINE_VERSION,
        "sequence_type_policy": "polymer-residues-v1", "length_policy": "polymer-length-v4",
        "source_scope": "asymmetric_unit", "structure_model_policy": "first",
        "ca_radius": 12.0, "all_atoms_radius": 5.0,
    }


def profile_id(profile=None):
    return hashlib.sha256(canonical_json(geometry_profile() if profile is None else profile)).hexdigest()


def _require(condition, message):
    if not condition:
        raise InputValidationError("INVALID_STRUCTURE_FACTS", message)


def _text(value, maximum=512, *, empty=False):
    return isinstance(value, str) and (empty or bool(value)) and len(value) <= maximum and not any(ord(c) < 32 or ord(c) == 127 for c in value)


def _integer(value, minimum=0):
    return type(value) is int and value >= minimum


def validate_facts(value):
    """Validate the complete scientific payload before annotation mutates a copy."""
    fields = {"schema_version", "producer_version", "profile", "profile_id", "structure_identity",
              "source_sha256", "source_expanded_bytes", "input_format", "input_kind", "warnings", "chains", "pairs"}
    _require(isinstance(value, dict) and set(value) == fields, "Structure fields are invalid.")
    _require(value["schema_version"] == STRUCTURE_FACTS_SCHEMA_VERSION, "Unsupported structure facts schema.")
    _require(value["profile"] == geometry_profile() and value["profile_id"] == profile_id(), "Structure profile is incompatible.")
    _require(_text(value["producer_version"], 64), "Producer version is invalid.")
    _require(isinstance(value["source_sha256"], str) and re.fullmatch(r"[0-9a-f]{64}", value["source_sha256"]), "Source hash is invalid.")
    _require(_integer(value["source_expanded_bytes"], 1), "Source size is invalid.")
    _require(value["input_format"] in ("pdb", "mmcif") and value["input_kind"] in (
        ("mmcif", "nextgen_enriched_mmcif") if value["input_format"] == "mmcif" else ("pdb",)
    ), "Source format is invalid.")
    try:
        identity = StructureIdentity.from_mapping(value["structure_identity"])
        _require(identity.source == "pdb" and value["structure_identity"] == identity.as_dict(), "An official canonical PDB identity is required.")
    except (TypeError, ValueError, AttributeError) as exc:
        raise InputValidationError("INVALID_STRUCTURE_FACTS", "Structure identity is invalid.") from exc
    warnings = value["warnings"]
    _require(isinstance(warnings, list) and len(warnings) <= 100_000, "Warnings are invalid.")
    for warning in warnings:
        _require(isinstance(warning, dict) and set(warning) == {"code", "message"}
                 and _text(warning["code"], 128) and _text(warning["message"], 4000), "Warning is invalid.")
    chains = value["chains"]
    _require(isinstance(chains, list) and 0 < len(chains) <= MAX_CHAINS, "Chain count is invalid.")
    lookup = {}
    for chain in chains:
        _require(isinstance(chain, dict) and set(chain) == CHAIN_FIELDS, "Chain fields are invalid.")
        _require(_text(chain["chain_id"]) and chain["chain_id"].strip() == chain["chain_id"], "Chain ID is invalid.")
        expected = ChainIdentity(identity.key, identity.display_id, chain["chain_id"], 1)
        try:
            _require(ChainIdentity.from_mapping(chain["chain_identity"]) == expected, "Chain identity is inconsistent.")
        except (TypeError, ValueError, AttributeError) as exc:
            raise InputValidationError("INVALID_STRUCTURE_FACTS", "Chain identity is invalid.") from exc
        _require(chain["chain_identity"] == expected.as_dict() and type(chain["model_index"]) is int
                 and chain["model_index"] == 1 and chain["structure_key"] == identity.key
                 and chain["structure_display_id"] == identity.display_id
                 and chain["unique_chain_id"] == expected.node_id(include_model=False), "Chain identity is inconsistent.")
        _require(chain["unique_chain_id"] not in lookup, "Duplicate chain identity.")
        lookup[chain["unique_chain_id"]] = chain
        _require(_integer(chain["aa_len"]) and _integer(chain["nt_len"]), "Polymer lengths are invalid.")
        seq = chain["structure_sequence"]
        _require(isinstance(seq, str) and len(seq) <= 10_000_000 and re.fullmatch(r"[A-Z]*", seq) is not None, "Structural sequence is invalid.")
        _require(chain["structure_molecule_type"] in ("Unknown", "Protein", "DNA", "RNA", "DNA/RNA", "Nucleic Acid"), "Structural type is invalid.")
        _require(isinstance(chain["embedded_annotations"], dict), "Embedded evidence must be an object.")
        validate_embedded_annotations(chain, chain["chain_id"], error_code="INVALID_STRUCTURE_FACTS")
    pairs = value["pairs"]
    _require(isinstance(pairs, list) and len(pairs) <= MAX_PAIRS, "Pair count is invalid.")
    seen = set()
    for pair in pairs:
        _require(isinstance(pair, dict) and set(pair) == {"chain_a", "chain_b", "all_atoms_count", "ca_neighbors"}, "Pair fields are invalid.")
        a, b = pair["chain_a"], pair["chain_b"]
        _require(isinstance(a, str) and isinstance(b, str) and a in lookup and b in lookup and a != b, "Pair endpoint is invalid.")
        key = tuple(sorted((a, b)))
        _require(key not in seen, "Duplicate chain pair.")
        seen.add(key)
        _require(_integer(pair["all_atoms_count"], 1) and _integer(pair["ca_neighbors"]), "Pair counts are invalid.")
    return value


def materialize_facts(value, *, use_embedded_sifts=True):
    from .network_annotations import apply_embedded_annotations, embedded_annotation_counts
    value = validate_facts(value)
    identity = value["structure_identity"]
    chains = copy.deepcopy(value["chains"])
    for chain in chains:
        chain.update(molecule_name="UNKNOWN", molecule_type="UNKNOWN", uniprot_id=None)
    apply_embedded_annotations(chains, {"by_chain": {c["chain_id"]: c["embedded_annotations"] for c in chains}},
                               use_embedded_sifts=use_embedded_sifts)
    return {
        "pdb_id": identity["display_id"], "file_path": "pdb:" + identity["canonical_id"],
        "structure_identity": dict(identity), "atom_data": chains,
        "input_format": value["input_format"], "input_kind": value["input_kind"],
        "input_warnings": copy.deepcopy(value["warnings"]),
        "embedded_annotation_counts": embedded_annotation_counts(chains),
    }


def validate_contacts(facts, rows):
    """Validate ordered contact rows in bounded memory, including all pair counts.

    Order is candidate-pair order, then partner B site, then partner A site,
    matching contact_atoms.unique_contacts. Strict ordering detects duplicates.
    """
    pairs = {(p["chain_a"], p["chain_b"]): (i, p["all_atoms_count"]) for i, p in enumerate(facts["pairs"])}
    counts = Counter()
    previous = None
    for row in rows:
        _require(isinstance(row, (list, tuple)) and len(row) == 11, "Contact row is invalid.")
        _require(all(_text(row[i], empty=i in (3, 7)) for i in range(10)), "Contact labels are invalid.")
        pair = (row[0], row[1])
        _require(pair in pairs, "Contact endpoints are invalid.")
        for number, insertion, label in ((row[2], row[3], row[4]), (row[6], row[7], row[8])):
            _require(re.fullmatch(r"-?[0-9]+", number) is not None and label.startswith(number + insertion + ":"), "Contact residue identity is inconsistent.")
        distance = row[10]
        _require(type(distance) in (int, float) and math.isfinite(distance) and 0 <= distance <= facts["profile"]["all_atoms_radius"], "Contact distance is invalid.")
        key = (pairs[pair][0], (row[6], row[7], row[9]), (row[2], row[3], row[5]))
        _require(previous is None or key > previous, "Contact rows are duplicated or out of order.")
        previous = key
        counts[pair] += 1
        _require(counts[pair] <= pairs[pair][1], "Contact row count exceeds the measured count.")
        yield row
    _require(all(counts[pair] == count for pair, (_, count) in pairs.items()), "Contact rows are incomplete.")


def facts_id(facts):
    """Identity of prepared source contents under one scientific profile."""
    return hashlib.sha256(canonical_json([facts["schema_version"], facts["profile_id"], facts["source_sha256"]])).hexdigest()


def prepared_provenance(facts, supplied=None):
    entries = [{"pdb_id": f["structure_identity"]["canonical_id"], "revision_id": facts_id(f), "source_sha256": f["source_sha256"]}
               for f in facts]
    value = supplied if supplied is not None else {
        "schema_version": STRUCTURE_FACTS_SCHEMA_VERSION, "profile_id": profile_id(),
        "catalog_id": None, "checked_at": None, "entries": entries,
    }
    validate_prepared_provenance(value)
    _require(value["entries"] == entries and value["profile_id"] == profile_id(), "Prepared provenance differs from analyzed inputs.")
    return copy.deepcopy(value)


def validate_prepared_provenance(value):
    """Validate path-free public provenance; raw inputs use null."""
    from .structure_identity import canonical_pdb_id
    if value is None:
        return None
    _require(isinstance(value, dict) and set(value) == {"schema_version", "profile_id", "catalog_id", "checked_at", "entries"}, "Prepared provenance fields are invalid.")
    _require(value["schema_version"] == STRUCTURE_FACTS_SCHEMA_VERSION, "Prepared provenance schema is invalid.")
    _require(isinstance(value["profile_id"], str) and re.fullmatch(r"[a-f0-9]{64}", value["profile_id"]), "Prepared profile identity is invalid.")
    if value["catalog_id"] is None:
        _require(value["checked_at"] is None, "A check timestamp requires a catalog identity.")
    else:
        _require(isinstance(value["catalog_id"], str) and re.fullmatch(r"[a-f0-9]{64}", value["catalog_id"]), "Catalog identity is invalid.")
        try:
            timestamp = datetime.fromisoformat(value["checked_at"])
            _require(timestamp.tzinfo is not None, "Catalog timestamp needs a timezone.")
        except (TypeError, ValueError) as exc:
            raise InputValidationError("INVALID_STRUCTURE_FACTS", "Catalog timestamp is invalid.") from exc
    _require(isinstance(value["entries"], list) and value["entries"], "Prepared source provenance is empty.")
    seen = set()
    for entry in value["entries"]:
        _require(isinstance(entry, dict) and set(entry) == {"pdb_id", "revision_id", "source_sha256"}, "Prepared source fields are invalid.")
        try:
            _require(canonical_pdb_id(entry["pdb_id"]) == entry["pdb_id"], "Prepared source identity is not canonical.")
        except ValueError as exc:
            raise InputValidationError("INVALID_STRUCTURE_FACTS", "Prepared source identity is invalid.") from exc
        _require(entry["pdb_id"] not in seen, "Prepared source identity is duplicated.")
        seen.add(entry["pdb_id"])
        _require(all(isinstance(entry[k], str) and re.fullmatch(r"[a-f0-9]{64}", entry[k]) for k in ("source_sha256", "revision_id")), "Prepared source hashes are invalid.")
        _require(entry["revision_id"] == facts_id({"schema_version": value["schema_version"], "profile_id": value["profile_id"], "source_sha256": entry["source_sha256"]}), "Prepared revision identity is inconsistent.")
    return value
