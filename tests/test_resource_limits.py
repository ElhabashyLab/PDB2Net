import json
from pathlib import Path

import pytest

from pdb2net import pipeline
from pdb2net.input_contract import InputValidationError


def _write_bytes(path: Path, size: int) -> None:
    path.write_bytes(b"x" * size)


def _set_limits(monkeypatch, **overrides) -> None:
    limits = {
        "max_input_files": None,
        "max_total_input_bytes": None,
        "max_single_input_bytes": None,
        "max_processing_batch_bytes": None,
        **overrides,
    }
    monkeypatch.setitem(pipeline.config, "resource_limits", limits)


def test_input_limits_fail_before_parsing(tmp_path: Path, monkeypatch) -> None:
    first = tmp_path / "a.cif"
    second = tmp_path / "b.cif"
    _write_bytes(first, 6)
    _write_bytes(second, 7)
    _set_limits(monkeypatch, max_total_input_bytes=12)

    with pytest.raises(InputValidationError) as exc_info:
        pipeline.inspect_input_files([str(first), str(second)])

    assert exc_info.value.code == "INPUT_TOTAL_BYTES_LIMIT_EXCEEDED"


def test_invalid_combined_graph_limit_fails_before_reference_loading(
    tmp_path: Path, monkeypatch
) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    _write_bytes(input_dir / "a.cif", 1)
    monkeypatch.setitem(pipeline.config, "output_path", str(tmp_path / "outputs"))
    monkeypatch.setitem(
        pipeline.config,
        "combined_graph_limits",
        {"max_nodes": "not-an-integer", "max_edges": None},
    )
    reference_checked = []
    monkeypatch.setattr(
        pipeline,
        "_validate_required_reference_files",
        lambda: reference_checked.append(True),
    )

    with pytest.raises(InputValidationError) as exc_info:
        pipeline.run_pipeline(str(input_dir))

    assert exc_info.value.code == "INVALID_COMBINED_GRAPH_LIMIT"
    assert reference_checked == []


def test_processing_batches_are_raw_size_bounded_and_stable(tmp_path: Path, monkeypatch) -> None:
    paths = [tmp_path / f"{name}.cif" for name in "abc"]
    for path, size in zip(paths, [6, 4, 7]):
        _write_bytes(path, size)
    _set_limits(monkeypatch, max_processing_batch_bytes=10)
    file_paths = [str(path) for path in paths]
    inventory = pipeline.inspect_input_files(file_paths)

    batches = list(pipeline.create_processing_batches(file_paths, inventory))

    assert batches == [file_paths[:2], file_paths[2:]]


def test_parser_scheduler_never_submits_more_than_worker_count(tmp_path: Path, monkeypatch) -> None:
    paths = [tmp_path / f"{index}.cif" for index in range(5)]
    for path in paths:
        _write_bytes(path, 1)

    observed_pending = []

    class ImmediateFuture:
        def __init__(self, fn, argument):
            self._result = fn(argument)

        def result(self):
            return self._result

    class RecordingExecutor:
        def __init__(self, max_workers):
            self.max_workers = max_workers

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, traceback):
            return False

        def submit(self, fn, argument):
            return ImmediateFuture(fn, argument)

    def fake_wait(pending, return_when):
        observed_pending.append(len(pending))
        first = next(iter(pending))
        return {first}, set(pending) - {first}

    monkeypatch.setitem(pipeline.config, "workers", {"parsing": 2, "blast_threads": 1})
    monkeypatch.setattr(pipeline, "ProcessPoolExecutor", RecordingExecutor)
    monkeypatch.setattr(pipeline, "wait", fake_wait)
    monkeypatch.setattr(
        pipeline,
        "process_single_file",
        lambda path: {"file_path": path, "pdb_id": Path(path).stem, "atom_data": []},
    )

    parsed = pipeline._parse_input_files([str(path) for path in paths])

    assert [entry["pdb_id"] for entry in parsed] == [str(index) for index in range(5)]
    assert max(observed_pending) == 2


def test_compaction_preserves_network_metadata_but_drops_atoms() -> None:
    structures = [
        {
            "pdb_id": "TST1",
            "file_path": "/inputs/tst1.cif",
            "atom_data": [
                {
                    "chain_id": "A",
                    "unique_chain_id": "TST1:A",
                    "molecule_type": "Protein",
                    "molecule_name": "Example",
                    "uniprot_id": "P12345",
                    "annotation_source": "blast_swissprot",
                    "residues": [
                        {
                            "residue_name": "ALA",
                            "atoms": [{"atom_name": "CA", "coordinates": [0.0, 0.0, 0.0]}],
                        }
                    ],
                }
            ],
        }
    ]

    compacted = pipeline._compact_structure_summaries(structures)

    chain = compacted[0]["atom_data"][0]
    assert chain["uniprot_id"] == "P12345"
    assert chain["annotation_source"] == "blast_swissprot"
    assert chain["aa_len"] == 1
    assert chain["nt_len"] == 0
    assert "residues" not in chain


def test_annotation_summary_counts_sifts_and_database_hits() -> None:
    summary = pipeline._annotation_summary(
        [
            {
                "pdb_id": "1ABC",
                "file_path": "/inputs/1abc.cif",
                "atom_data": [
                    {"chain_id": "A", "uniprot_id": "P11111"},
                    {
                        "chain_id": "B",
                        "uniprot_id": "P22222",
                        "annotation_source": "blast_swissprot",
                        "matched_database": "Swiss-Prot",
                        "annotation_confidence": "high",
                    },
                    {
                        "chain_id": "C",
                        "uniprot_id": None,
                        "annotation_source": "diamond_uniref90",
                        "matched_database": "UniRef90",
                        "matched_id": "UniRef90_P33333",
                        "annotation_confidence": "medium",
                    },
                    {"chain_id": "D", "uniprot_id": None},
                ],
            }
        ]
    )

    assert summary["chains_total"] == 4
    assert summary["chains_annotated"] == 3
    assert summary["chains_unannotated"] == 1
    assert summary["by_source"] == {
        "blast_swissprot": 1,
        "diamond_uniref90": 1,
        "sifts": 1,
    }
    assert summary["by_database"] == {"Swiss-Prot": 1, "UniProtKB": 1, "UniRef90": 1}
    assert summary["by_confidence"] == {"high": 1, "medium": 1, "not_reported": 1}


def test_staged_pipeline_runs_expensive_stages_per_batch_and_exports_once(tmp_path: Path, monkeypatch) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    for name in ["a.cif", "b.cif"]:
        _write_bytes(input_dir / name, 6)

    monkeypatch.setitem(pipeline.config, "output_path", str(tmp_path / "outputs"))
    monkeypatch.setitem(
        pipeline.config,
        "networks",
        {
            "chain_per_pdb": True,
            "combined_chain_network": False,
            "protein_per_pdb": False,
            "combined_protein_network": False,
        },
    )
    monkeypatch.setitem(pipeline.config, "combined_graph_limits", {"max_nodes": None, "max_edges": None})
    monkeypatch.setitem(pipeline.config, "reference_manifest_id", "mini-ref-v1")
    _set_limits(monkeypatch, max_processing_batch_bytes=6)
    monkeypatch.setattr(pipeline, "_validate_required_reference_files", lambda: None)
    monkeypatch.setattr(pipeline, "_reference_summary", lambda: {"manifest_id": "mini-ref-v1"})

    parsed_batches = []
    annotated_batches = []
    distance_batches = []
    exported = []

    def fake_parse(paths):
        parsed_batches.append([Path(path).name for path in paths])
        return [
            {
                "pdb_id": Path(path).stem.upper(),
                "file_path": path,
                "atom_data": [
                    {
                        "chain_id": "A",
                        "unique_chain_id": f"{Path(path).stem.upper()}:A",
                        "molecule_type": "Protein",
                        "molecule_name": "Example",
                        "uniprot_id": "P12345",
                        "residues": [
                            {
                                "residue_name": "ALA",
                                "atoms": [{"atom_name": "CA", "coordinates": [0.0, 0.0, 0.0]}],
                            }
                        ],
                    }
                ],
            }
            for path in paths
        ]

    monkeypatch.setattr(pipeline, "_parse_input_files", fake_parse)
    monkeypatch.setattr(
        pipeline,
        "process_molecule_info",
        lambda data: annotated_batches.append([entry["pdb_id"] for entry in data]),
    )
    monkeypatch.setattr(pipeline, "_run_blast_annotation", lambda data: None)
    monkeypatch.setattr(
        pipeline,
        "calculate_distances_with_ckdtree",
        lambda data: distance_batches.append([entry["pdb_id"] for entry in data]) or [],
    )

    def fake_export(data, results, network_config, output_paths, timings):
        exported.append(data)
        assert all("residues" not in chain for structure in data for chain in structure["atom_data"])
        return []

    monkeypatch.setattr(pipeline, "_export_network_outputs", fake_export)

    output_paths = pipeline.run_pipeline(str(input_dir))

    assert len(parsed_batches) == 2
    assert annotated_batches == [["A"], ["B"]]
    assert distance_batches == [["A"], ["B"]]
    assert len(exported) == 1
    manifest = json.loads(Path(output_paths.manifest_file).read_text(encoding="utf-8"))
    assert manifest["resources"]["processing"]["batches"] == 2
    assert manifest["resources"]["processing"]["parsing_workers"] == 1
    assert manifest["resources"]["input"]["total_bytes"] == 12
    assert manifest["references"]["manifest_id"] == "mini-ref-v1"


def test_combined_chain_network_over_limit_is_reported(monkeypatch) -> None:
    exported = []
    monkeypatch.setattr(
        pipeline,
        "create_cytoscape_network",
        lambda edges, title, output_path, nodes_data=None: exported.append(title),
    )
    combined_data = [
        {
            "pdb_id": pdb_id,
            "file_path": f"/inputs/{pdb_id}.cif",
            "atom_data": [
                {
                    "chain_id": "A",
                    "unique_chain_id": f"{pdb_id}:A",
                    "molecule_type": "Protein",
                    "molecule_name": "Example",
                    "uniprot_id": "P12345",
                    "aa_len": 1,
                    "nt_len": 0,
                }
            ],
        }
        for pdb_id in ["PDB1", "PDB2"]
    ]

    skipped = pipeline._create_linked_identity_network(
        [],
        combined_data,
        "/tmp/combined",
        combined_graph_limits={"max_nodes": 1, "max_edges": 10},
    )

    assert exported == []
    assert skipped[0]["network_kind"] == "combined_chain_network"
    assert skipped[0]["counts"] == {"nodes": 2, "edges": 1}
