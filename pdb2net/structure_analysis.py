"""Scientific preparation and analysis APIs, independent of store ownership."""
from __future__ import annotations

from datetime import datetime
import time

from . import __version__, pipeline
from .config_loader import activate_config, config as active_config
from .contact_atoms import ContactView, count_contacts, unique_contacts
from .detailed_results_exporter import DetailedInteractionBudget, export_detailed_interactions
from .distances import chain_pairs, select_interaction
from .file_parser import input_file_signature
from .input_contract import InputValidationError
from .network_annotations import network_annotation_config
from .outputs import (
    collect_generated_outputs, collect_web_outputs, create_run_output_paths,
    reserve_web_output_directory, write_failed_run_manifest, write_run_manifest,
    write_run_summary, write_runtime_analysis,
)
from .residue_types import count_polymer_lengths
from .structure_facts import (
    CHAIN_FIELDS, MAX_CHAINS, MAX_PAIRS, MAX_REQUEST_CHAINS, MAX_REQUEST_PAIRS,
    STRUCTURE_FACTS_SCHEMA_VERSION, geometry_profile, materialize_facts,
    profile_id, prepared_provenance, validate_contacts, validate_facts,
)
from .uniprot_matcher import classify_molecule_type, extract_sequence_from_parsed_data


def prepare_structure(source, profile=None, contact_sink=None):
    """Parse one official PDB structure and optionally stream every atom contact.

    The caller owns file I/O and bounds on streamed contact bytes/rows. An error
    from its sink aborts preparation; no partial facts are returned. Reference
    files and sequence-search databases are never needed by this operation.
    """
    profile = geometry_profile() if profile is None else profile
    if profile != geometry_profile() or active_config.get("structure_model_policy", "first") != "first":
        raise InputValidationError("STRUCTURE_PROFILE_UNSUPPORTED", "Preparation requires the fixed first-model geometry profile.")
    source = str(source)
    signature = input_file_signature(source)
    limits = pipeline._resource_limits()
    maximum = limits["max_single_input_bytes"]
    if maximum is not None and signature[2] > maximum:
        raise InputValidationError("INPUT_FILE_BYTES_LIMIT_EXCEEDED", "Structure source exceeds its byte limit.")
    structure = pipeline.process_single_file(source, expected_signature=signature,
                                             maximum_expanded_bytes=limits["max_single_input_expanded_bytes"])
    if structure["structure_identity"]["source"] != "pdb":
        raise InputValidationError("INVALID_STRUCTURE_FACTS", "Preparation requires an official PDB identity.")
    if not structure["atom_data"]:
        raise InputValidationError("NO_VALID_CHAINS", "No supported contact chains were found.")
    if len(structure["atom_data"]) > MAX_CHAINS:
        raise InputValidationError("STRUCTURE_FACTS_LIMIT_EXCEEDED", "Too many chains to prepare.")
    chains = []
    for raw in structure["atom_data"]:
        chain = {key: raw[key] for key in CHAIN_FIELDS if key in raw}
        aa, nt = count_polymer_lengths(raw["residues"])
        chain.update(structure_sequence=extract_sequence_from_parsed_data(raw),
                     structure_molecule_type=classify_molecule_type(raw), aa_len=aa, nt_len=nt,
                     embedded_annotations=raw.get("embedded_annotations", {}))
        chains.append(chain)
    views = {(c["unique_chain_id"], kind): ContactView(c, kind)
             for c in structure["atom_data"] for kind in ("all_atoms", "ca")}
    pairs = []
    for left, right in chain_pairs(structure["atom_data"]):
        a, b = left["unique_chain_id"], right["unique_chain_id"]
        all_count = count_contacts(views[a, "all_atoms"], views[b, "all_atoms"], profile["all_atoms_radius"])
        if not all_count:
            continue
        ca_count = count_contacts(views[a, "ca"], views[b, "ca"], profile["ca_radius"])
        pairs.append({"chain_a": a, "chain_b": b, "all_atoms_count": all_count, "ca_neighbors": ca_count})
        if len(pairs) > MAX_PAIRS:
            raise InputValidationError("STRUCTURE_FACTS_LIMIT_EXCEEDED", "Too many contact candidates to prepare.")
        if contact_sink is not None:
            emitted = 0
            for atom_a, atom_b, distance in unique_contacts(views[a, "all_atoms"], views[b, "all_atoms"], profile["all_atoms_radius"]):
                contact_sink([a, b, atom_a.site[0], atom_a.site[1], atom_a.residue, atom_a.atom_name,
                              atom_b.site[0], atom_b.site[1], atom_b.residue, atom_b.atom_name, distance])
                emitted += 1
            if emitted != all_count:
                raise InputValidationError("CONTACT_COUNT_MISMATCH", "Detailed contacts differ from measured counts.")
    if input_file_signature(source) != signature:
        raise InputValidationError("INPUT_CHANGED_DURING_PROCESSING", "Structure source changed during preparation.")
    return validate_facts({
        "schema_version": STRUCTURE_FACTS_SCHEMA_VERSION, "producer_version": __version__,
        "profile": dict(profile), "profile_id": profile_id(profile),
        "structure_identity": structure["structure_identity"], "source_sha256": structure["source_sha256"],
        "source_expanded_bytes": structure["source_expanded_bytes"],
        "input_format": structure["input_format"], "input_kind": structure["input_kind"],
        "warnings": [{k: w[k] for k in ("code", "message")} for w in structure["input_warnings"]],
        "chains": chains, "pairs": pairs,
    })


def analyze_prepared(structures, config, output_paths=None, contact_sources=None, provenance=None, *, web_output_dir=None):
    """Annotate validated facts anew and export through the ordinary Core writers.

    Each contact source is a zero-argument iterator factory, keyed by canonical
    PDB ID. Call in a dedicated process with one immutable reference generation.
    The caller must bound file decoding before handing dictionaries to Core.
    """
    activate_config(config)
    started = datetime.now().isoformat(timespec="seconds")
    start = time.monotonic()
    reserved = attempted = False
    timings = pipeline.PipelineTimings()
    try:
        pipeline._validate_output_root()
        pipeline._validate_analysis_config()
        pipeline._validate_web_reference_manifest(web_output_dir)
        annotation_config = network_annotation_config()
        profile = geometry_profile()
        if (active_config.get("structure_model_policy", "first") != "first"
                or any(active_config["distance_thresholds"][k] != profile[k] for k in ("ca_radius", "all_atoms_radius"))):
            raise InputValidationError("STRUCTURE_PROFILE_UNSUPPORTED", "Prepared analysis requires the stored model policy and radii.")
        if web_output_dir:
            attempted = True
            reserve_web_output_directory(web_output_dir)
            reserved = True
        pipeline._combined_graph_limits()
        pipeline._validate_required_reference_files()
        limits = pipeline._resource_limits()
        payloads, materialized, seen = [], [], set()
        total_chains = total_pairs = 0
        for payload in structures:
            payload = validate_facts(payload)
            total_chains += len(payload["chains"])
            total_pairs += len(payload["pairs"])
            if total_chains > MAX_REQUEST_CHAINS or total_pairs > MAX_REQUEST_PAIRS:
                raise InputValidationError("STRUCTURE_FACTS_LIMIT_EXCEEDED", "Prepared inputs exceed the aggregate chain or pair limit.")
            identity = payload["structure_identity"]["canonical_id"]
            if identity in seen:
                raise InputValidationError("DUPLICATE_STRUCTURE_IDENTITY", "Requested structures are duplicated.")
            seen.add(identity)
            if limits["max_input_files"] is not None and len(seen) > limits["max_input_files"]:
                raise InputValidationError("INPUT_FILE_COUNT_LIMIT_EXCEEDED", "Too many prepared inputs.")
            payloads.append(payload)
            materialized.append(materialize_facts(payload, use_embedded_sifts=annotation_config["use_embedded_sifts"]))
        if not payloads:
            raise InputValidationError("NO_PARSEABLE_STRUCTURES", "No prepared inputs were supplied.")
        provenance = prepared_provenance(payloads, provenance)
        output_paths = output_paths or create_run_output_paths(active_config["output_path"])
        before = time.monotonic()
        pipeline.process_molecule_info(materialized)
        timings.sifts = time.monotonic() - before
        before = time.monotonic()
        pipeline._run_blast_annotation(materialized)
        timings.blast = time.monotonic() - before
        interactions = []
        detailed = bool(active_config.get("export_detailed_interactions"))
        budget = DetailedInteractionBudget(max_rows=limits["max_detailed_interaction_rows"],
                                          max_bytes=limits["max_detailed_interaction_bytes"],
                                          min_free_bytes=limits["min_output_free_bytes"])
        for payload, structure in zip(payloads, materialized):
            lookup = {c["unique_chain_id"]: c for c in structure["atom_data"]}
            edges = []
            for pair in payload["pairs"]:
                edge = select_interaction(structure, lookup[pair["chain_a"]], lookup[pair["chain_b"]],
                                          pair["all_atoms_count"], pair["ca_neighbors"])
                if edge is not None:
                    edges.append(edge)
            interactions.extend(edges)
            if detailed:
                factory = (contact_sources or {}).get(payload["structure_identity"]["canonical_id"])
                if not callable(factory):
                    raise InputValidationError("STRUCTURE_CONTACTS_MISSING", "Detailed analysis requires a contact stream for every structure.")
                export_detailed_interactions(structure, edges, output_paths.distances_dir, budget=budget,
                                             contact_rows=validate_contacts(payload, factory()))
        skipped = pipeline._export_network_outputs(materialized, interactions, active_config["networks"], output_paths, timings)
        warnings = pipeline._embedded_annotation_warnings(materialized) + [
            {"code": "COMBINED_GRAPH_LIMIT_EXCEEDED", "message": f"Skipped {item['name']} because configured graph limits were exceeded.", "output": item["name"]}
            for item in skipped
        ]
        total = time.monotonic() - start
        write_runtime_analysis(output_paths.log_file, timings.as_dict(), total)
        write_run_manifest(
            output_paths.manifest_file, input_files=[s["file_path"] for s in materialized], output_paths=output_paths,
            config_snapshot=pipeline._config_snapshot(active_config["networks"]), status="success",
            started_at=started, finished_at=datetime.now().isoformat(timespec="seconds"), total_time=total,
            generated_outputs=collect_generated_outputs(output_paths),
            extra_counts={"structures": len(payloads), "chains": sum(len(s["atom_data"]) for s in materialized),
                          "interactions": len(interactions), "skipped_outputs": len(skipped)},
            annotations=pipeline._annotation_summary(materialized), references=pipeline._reference_summary(),
            resources={"input": {"files": len(payloads)}, "processing": {"mode": "prepared", "parsing_workers": 0},
                       "detailed_interactions": budget.as_dict(enabled=detailed)},
            skipped_outputs=skipped, warnings=warnings, identities=pipeline._identity_summary(materialized),
            structure_inputs=pipeline._structure_input_summary(materialized), prepared_inputs=provenance,
        )
        write_run_summary(output_paths)
        if web_output_dir:
            collect_web_outputs(output_paths, web_output_dir, web_output_prepared=True)
        return output_paths
    except (Exception, KeyboardInterrupt) as exc:
        failure = InputValidationError("RUN_INTERRUPTED", "Run interrupted by user.") if isinstance(exc, KeyboardInterrupt) else exc
        base_output_path = str(active_config.get("output_path") or "").strip()
        if base_output_path:
            try:
                output_paths = write_failed_run_manifest(base_output_path, input_path="prepared",
                                                        config_snapshot=pipeline._config_snapshot(active_config.get("networks", {})),
                                                        error=failure, started_at=started, output_paths=output_paths)
                if web_output_dir and (reserved or not attempted):
                    collect_web_outputs(output_paths, web_output_dir, web_output_prepared=reserved)
            except Exception:
                if not isinstance(exc, KeyboardInterrupt):
                    raise
        raise
    finally:
        pipeline.tree_cache.clear()
        pipeline.coords_cache.clear()
