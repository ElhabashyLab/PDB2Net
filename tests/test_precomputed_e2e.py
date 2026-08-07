from __future__ import annotations

import gzip
import json
from pathlib import Path

from pdb2net import file_parser, pipeline, reference_data
from pdb2net import precomputed_store as store
from pdb2net import unknown_molecule_uniprot
from pdb2net.config_loader import config


def _mini_pdb_text() -> str:
    lines = ["HEADER    MINI PRECOMPUTED TEST                 01-JAN-00   1ABC\n"]
    serial = 1
    for chain, y in (("A", 0.0), ("B", 1.0)):
        for residue in range(1, 11):
            x = (residue - 1) * 3.8
            lines.append(
                f"ATOM  {serial:5d}  CA  ALA {chain}{residue:4d}    "
                f"{x:8.3f}{y:8.3f}{0.0:8.3f}  1.00 20.00           C  \n"
            )
            serial += 1
        lines.append("TER\n")
    lines.append("END\n")
    return "".join(lines)


def _network_view(path: Path) -> tuple[list[dict], list[dict]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    aspects: dict[str, list] = {}
    for block in document:
        if isinstance(block, dict):
            for key, value in block.items():
                if isinstance(value, list):
                    aspects[key] = value
    nodes = aspects["nodes"]
    raw_id_to_uid = {node["id"]: node["v"]["name"] for node in nodes}
    normalized_nodes = sorted([
        {
            "id": node["v"]["name"],
            "molecule_type": node["v"].get("molecule_type"),
            "uniprot_id": node["v"].get("uniprot_id"),
            "aa_len": node["v"].get("aa_len"),
            "nt_len": node["v"].get("nt_len"),
        }
        for node in nodes
    ], key=lambda item: item["id"])
    normalized_edges = sorted([
        {
            "source": raw_id_to_uid[edge["s"]],
            "target": raw_id_to_uid[edge["t"]],
            "interaction": edge["v"]["interaction"],
            "all_atoms_count": edge["v"]["all_atoms_count"],
        }
        for edge in aspects["edges"]
    ], key=lambda item: (item["source"], item["target"]))
    return normalized_nodes, normalized_edges


def test_real_gemmi_to_store_to_headless_cx2_matches_raw_pipeline(
    tmp_path: Path, monkeypatch
) -> None:
    references = tmp_path / "references"
    references.mkdir()
    pdb_fasta = references / "pdb_seqres.txt"
    pdb_fasta.write_text(
        ">1abc_A mol:protein length:10 Chain Alpha\nAAAAAAAAAA\n"
        ">1abc_B mol:protein length:10 Chain Beta\nGGGGGGGGGG\n",
        encoding="utf-8",
    )
    sifts = references / "pdb_chain_uniprot.tsv"
    sifts.write_text("1abc\tA\tP12345\n1abc\tB\tQ99999\n", encoding="utf-8")
    uniprot = references / "uniprot_sprot.fasta"
    uniprot.write_text(
        ">sp|P12345|PROTA Protein Alpha\nAAAAAAAAAA\n"
        ">sp|Q99999|PROTB Protein Beta\nGGGGGGGGGG\n",
        encoding="utf-8",
    )

    archive = tmp_path / "archive" / "divided" / "pdb" / "ab"
    archive.mkdir(parents=True)
    source = archive / "pdb1abc.ent.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(_mini_pdb_text())

    monkeypatch.setitem(config, "pdb_fasta_path", str(pdb_fasta))
    monkeypatch.setitem(config, "sifts_tsv_path", str(sifts))
    monkeypatch.setitem(config, "uniprot_fasta_path", str(uniprot))
    monkeypatch.setitem(config, "reference_manifest_id", "mini-original-layout-v1")
    monkeypatch.setitem(config, "structure_model_policy", "first")
    monkeypatch.setitem(config, "open_in_cytoscape", False)
    monkeypatch.setitem(config, "layout_mode", "python_fast")
    monkeypatch.setitem(config, "export_detailed_interactions", False)
    monkeypatch.setitem(
        config,
        "networks",
        {
            "chain_per_pdb": True,
            "protein_per_pdb": False,
            "combined_chain_network": False,
            "combined_protein_network": False,
        },
    )
    monkeypatch.setitem(
        config,
        "resource_limits",
        {
            "max_input_files": 50,
            "max_total_input_bytes": 10_000_000,
            "max_single_input_bytes": 10_000_000,
            "max_processing_batch_bytes": 10_000_000,
        },
    )

    monkeypatch.setattr(file_parser, "PDB_FASTA_PATH", str(pdb_fasta))
    monkeypatch.setattr(file_parser, "VALID_PDB_IDS", None)
    monkeypatch.setattr(unknown_molecule_uniprot, "PDB_FASTA_PATH", str(pdb_fasta))
    monkeypatch.setattr(unknown_molecule_uniprot, "SIFTS_TSV_PATH", str(sifts))
    monkeypatch.setattr(unknown_molecule_uniprot, "UNIPROT_FASTA_PATH", str(uniprot))
    monkeypatch.setattr(unknown_molecule_uniprot, "_sifts_loaded", False)
    monkeypatch.setattr(unknown_molecule_uniprot, "_uniprot_loaded", False)
    monkeypatch.setattr(unknown_molecule_uniprot, "pdb_to_uniprot", {})
    monkeypatch.setattr(unknown_molecule_uniprot, "uniprot_dict", {})
    reference_data.load_valid_pdb_ids.cache_clear()
    reference_data.load_pdb_fasta_headers.cache_clear()
    reference_data.load_sifts_mapping.cache_clear()
    reference_data.load_uniprot_names.cache_clear()

    # Keep the test deterministic while exercising the real Gemmi parser and
    # data processor rather than process-pool mechanics.
    monkeypatch.setattr(
        pipeline,
        "_parse_input_files",
        lambda paths: [
            parsed
            for path in paths
            if (parsed := pipeline.process_single_file(path)) is not None
        ],
    )

    raw_output = tmp_path / "raw-output"
    monkeypatch.setitem(config, "output_path", str(raw_output))
    raw_run = pipeline.run_pipeline([str(source)])

    cache_root = tmp_path / "precomputed"
    report = store.precompute_directory(cache_root, tmp_path / "archive", recursive=True)
    assert report["written"] == 1
    assert report["failed"] == 0

    cached = store.load_entry(cache_root, "1abc")
    assert cached["counts"] == {"chains": 2, "interactions": 1}
    assert all("residues" not in chain for chain in cached["structure"]["atom_data"])
    assert cached["interactions"][0]["ca_neighbors"] >= 10

    assembled_output = tmp_path / "assembled-output"
    monkeypatch.setitem(config, "output_path", str(assembled_output))
    assembled_run = store.run_assemble_pipeline(cache_root, ["1ABC"])

    raw_network = Path(raw_run.chain_dir) / "Chain_Interaction_Network_1ABC.cx2"
    assembled_network = Path(assembled_run.chain_dir) / "Chain_Interaction_Network_1ABC.cx2"
    assert raw_network.is_file()
    assert assembled_network.is_file()
    assert _network_view(raw_network) == _network_view(assembled_network)

    second = store.precompute_directory(cache_root, tmp_path / "archive", recursive=True)
    assert second["cache_hits"] == 1
    assert second["written"] == 0
