from __future__ import annotations

import gzip
import json
from pathlib import Path
import re

import gemmi
import pytest

from pdb2net import precomputed
from pdb2net import network_annotations
from pdb2net import unknown_molecule_uniprot
from pdb2net.config_loader import config
from pdb2net.cytoscape import generate_nodes_from_atom_data
from pdb2net.file_parser import get_structure_identity
from pdb2net.input_contract import InputValidationError
from pdb2net.network_annotations import (
    apply_embedded_annotations,
    annotation_node_metadata,
    extract_embedded_annotations,
    network_annotation_config,
)
from pdb2net.pipeline import process_single_file
from pdb2net.structure_identity import StructureIdentity


NEXTGEN_8ALY = Path(__file__).parent / "fixtures" / "nextgen" / "pdb_00008aly_xyz-enrich.cif"


def _enriched_mmcif(*, multiple_best_uniprot: bool = False) -> str:
    pdb_text = (
        "HEADER    TEST                                    01-JAN-00   1ABC\n"
        "ATOM      1  CA  ALA Z   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ATOM      2  CA  GLY Z   2       1.000   0.000   0.000  1.00 20.00           C  \n"
        "TER\nEND\n"
    )
    base = gemmi.read_pdb_string(pdb_text).make_mmcif_document().as_string().replace(
        "data_string", "data_1abc", 1
    )
    block = gemmi.cif.read_string(base).sole_block()
    label_chain = block.get_mmcif_category("_atom_site.")["label_asym_id"][0]
    extra_uniprot = (
        f"1 {label_chain} Q99999 2 1 10 11 2 2 y 0.75\n"
        if multiple_best_uniprot
        else ""
    )
    return base + (
        "\nloop_\n"
        "_pdbx_sifts_unp_segments.entity_id\n"
        "_pdbx_sifts_unp_segments.asym_id\n"
        "_pdbx_sifts_unp_segments.unp_acc\n"
        "_pdbx_sifts_unp_segments.segment_id\n"
        "_pdbx_sifts_unp_segments.instance_id\n"
        "_pdbx_sifts_unp_segments.unp_start\n"
        "_pdbx_sifts_unp_segments.unp_end\n"
        "_pdbx_sifts_unp_segments.seq_id_start\n"
        "_pdbx_sifts_unp_segments.seq_id_end\n"
        "_pdbx_sifts_unp_segments.best_mapping\n"
        "_pdbx_sifts_unp_segments.identity\n"
        f"1 {label_chain} P12345 1 1 20 21 1 2 y 0.95\n"
        f"{extra_uniprot}#\n"
        "loop_\n"
        "_pdbx_sifts_xref_db_segments.entity_id\n"
        "_pdbx_sifts_xref_db_segments.asym_id\n"
        "_pdbx_sifts_xref_db_segments.xref_db\n"
        "_pdbx_sifts_xref_db_segments.xref_db_acc\n"
        "_pdbx_sifts_xref_db_segments.domain_name\n"
        "_pdbx_sifts_xref_db_segments.segment_id\n"
        "_pdbx_sifts_xref_db_segments.instance_id\n"
        "_pdbx_sifts_xref_db_segments.seq_id_start\n"
        "_pdbx_sifts_xref_db_segments.seq_id_end\n"
        f"1 {label_chain} Pfam PF00001 ? 1 1 1 2\n"
        f"1 {label_chain} CATH 1.10.8.10 1abcZ01 1 1 1 2\n"
        f"1 {label_chain} SCOP2B 123456 SF 1 1 1 2\n#\n"
    )


@pytest.mark.parametrize("compressed", [False, True])
def test_enriched_mmcif_is_detected_and_mapped_to_author_chain(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, compressed: bool
) -> None:
    suffix = ".cif.gz" if compressed else ".cif"
    source = tmp_path / f"pdb_00001abc_xyz-enrich{suffix}"
    content = _enriched_mmcif()
    if compressed:
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write(content)
    else:
        source.write_text(content, encoding="utf-8")
    monkeypatch.setitem(
        config,
        "network_annotations",
        {
            "use_embedded_sifts": True,
            "tooltip_fields": ["uniprot", "pfam", "cath", "scop2"],
            "max_tooltip_segments_per_database": 20,
        },
    )

    parsed = process_single_file(str(source))

    assert parsed is not None
    assert parsed["input_kind"] == "nextgen_enriched_mmcif"
    assert parsed["structure_identity"]["canonical_id"] == "pdb_00001abc"
    chain = parsed["atom_data"][0]
    assert chain["chain_id"] == "Z"
    assert chain["uniprot_id"] == "P12345"
    assert chain["embedded_uniprot_status"] == "unique_best_mapping"
    assert set(chain["embedded_annotations"]) == {"uniprot", "pfam", "cath", "scop2"}

    node = generate_nodes_from_atom_data(parsed["atom_data"], parsed["pdb_id"])[0]
    assert "UniProt mapping: P12345" in node["tooltip"]
    assert "Pfam: PF00001" in node["tooltip"]
    assert "CATH: 1.10.8.10" in node["tooltip"]
    assert "SCOP2: 123456" in node["tooltip"]
    assert json.loads(node["annotation_pfam"])[0]["accession"] == "PF00001"


def test_multiple_best_embedded_uniprot_mappings_remain_ambiguous(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "pdb_00001abc_xyz-enrich.cif"
    source.write_text(_enriched_mmcif(multiple_best_uniprot=True), encoding="utf-8")
    monkeypatch.setitem(config["network_annotations"], "use_embedded_sifts", True)

    parsed = process_single_file(str(source))

    assert parsed is not None
    chain = parsed["atom_data"][0]
    assert chain["embedded_uniprot_status"] == "ambiguous_multi_mapping"
    assert chain["embedded_uniprot_accessions"] == ["P12345", "Q99999"]
    assert chain["uniprot_id"] is None


def test_static_official_8aly_nextgen_fixture_preserves_real_chain_mapping_and_segments(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(config["network_annotations"], "use_embedded_sifts", True)

    parsed = process_single_file(str(NEXTGEN_8ALY))

    assert parsed is not None
    assert parsed["structure_identity"] == StructureIdentity(
        "pdb", "pdb_00008aly", "8aly"
    ).as_dict()
    assert parsed["input_kind"] == "nextgen_enriched_mmcif"
    assert parsed["embedded_annotation_counts"] == {
        "uniprot": 1,
        "pfam": 2,
        "cath": 0,
        "scop2": 0,
    }
    chain = parsed["atom_data"][0]
    assert chain["chain_id"] == "E"
    assert chain["uniprot_id"] == "Q9H2K2"
    assert [segment["accession"] for segment in chain["embedded_annotations"]["pfam"]] == [
        "PF00644",
        "PF07647",
    ]


def test_multiple_best_segments_for_one_accession_remain_unambiguous(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    content = _enriched_mmcif()
    label_match = re.search(r"^1\s+(\S+)\s+P12345\s+1\s+1\s+20\s+21", content, re.MULTILINE)
    assert label_match is not None
    content = content.replace(
        "#\nloop_\n_pdbx_sifts_xref_db_segments.entity_id",
        f"1 {label_match.group(1)} P12345 2 1 30 31 3 4 Y 0.90\n#\n"
        "loop_\n_pdbx_sifts_xref_db_segments.entity_id",
        1,
    )
    source = tmp_path / "1abc.cif"
    source.write_text(content, encoding="utf-8")
    monkeypatch.setitem(config["network_annotations"], "use_embedded_sifts", True)

    parsed = process_single_file(str(source))

    assert parsed is not None
    chain = parsed["atom_data"][0]
    assert chain["embedded_uniprot_status"] == "unique_best_mapping"
    assert chain["uniprot_id"] == "P12345"
    assert len(chain["embedded_annotations"]["uniprot"]) == 2


def test_invalid_embedded_segments_are_aggregated_and_never_form_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    invalid_rows = "".join(
        (
            "? A Q11111 2 1 1 2 1 2 Y 0.9\n",
            "1 A BAD 3 1 1 2 1 2 Y 0.9\n",
            "1 A Q11111 4 1 1 2 1 2 yes 0.9\n",
            "1 A Q11111 5 1 1 2 1 2 Y 1.1\n",
            "1 A Q11111 6 1 1 2 0 2 Y 0.9\n",
            "1 X Q11111 7 1 1 2 1 2 Y 0.9\n",
        )
    )
    content = _enriched_mmcif().replace(
        "#\nloop_\n_pdbx_sifts_xref_db_segments.entity_id",
        invalid_rows + "#\nloop_\n_pdbx_sifts_xref_db_segments.entity_id",
        1,
    )
    source = tmp_path / "1abc.cif"
    source.write_text(content, encoding="utf-8")
    monkeypatch.setitem(config["network_annotations"], "use_embedded_sifts", True)

    parsed = process_single_file(str(source))

    assert parsed is not None
    chain = parsed["atom_data"][0]
    assert chain["uniprot_id"] == "P12345"
    assert {segment["accession"] for segment in chain["embedded_annotations"]["uniprot"]} == {
        "P12345"
    }
    warnings = [
        warning
        for warning in parsed["input_warnings"]
        if warning["code"] == "INVALID_EMBEDDED_SIFTS_SEGMENT"
    ]
    assert sum(warning["count"] for warning in warnings) == 6
    assert {warning["reason"] for warning in warnings} == {
        "missing_required_field",
        "invalid_uniprot_accession",
        "invalid_best_mapping",
        "invalid_identity",
        "invalid_range",
        "unresolvable_label_asym_id",
    }


def test_poly_seq_scheme_ambiguity_is_not_overridden_by_atom_site_fallback() -> None:
    content = NEXTGEN_8ALY.read_text(encoding="utf-8").replace(
        "A 1 13 PHE 13 876 876 PHE PHE E . n",
        "A 1 13 PHE 13 876 876 PHE PHE F . n",
    )
    extracted = extract_embedded_annotations(gemmi.cif.read_string(content).sole_block())

    assert extracted["by_chain"] == {}
    assert extracted["invalid_segment_count"] == 3
    assert {warning["reason"] for warning in extracted["warnings"]} == {
        "unresolvable_label_asym_id"
    }


def test_atom_site_is_used_when_poly_seq_scheme_mapping_is_absent() -> None:
    content = NEXTGEN_8ALY.read_text(encoding="utf-8")
    start = content.index("loop_\n_pdbx_poly_seq_scheme.asym_id")
    end = content.index("loop_\n_pdbx_sifts_unp_segments.entity_id", start)
    content = content[:start] + content[end:]

    extracted = extract_embedded_annotations(gemmi.cif.read_string(content).sole_block())

    assert set(extracted["by_chain"]) == {"E"}
    assert extracted["by_chain"]["E"]["uniprot"][0]["accession"] == "Q9H2K2"


@pytest.mark.parametrize(
    ("headers", "expected_name"),
    [
        ({"1abc_A": {"info": "mol:protein length:2 PDB molecule name"}}, "PDB molecule name"),
        ({}, "UniProt: P12345"),
    ],
)
def test_embedded_external_mismatch_never_borrows_a_uniprot_name(
    monkeypatch: pytest.MonkeyPatch,
    headers: dict[str, dict[str, str]],
    expected_name: str,
) -> None:
    monkeypatch.setattr(
        unknown_molecule_uniprot,
        "_sifts_loaded_path",
        str(config["sifts_tsv_path"]),
    )
    monkeypatch.setattr(
        unknown_molecule_uniprot,
        "_uniprot_loaded_path",
        str(config["uniprot_fasta_path"]),
    )
    monkeypatch.setattr(
        unknown_molecule_uniprot,
        "pdb_to_sifts_segments",
        {"1abc_A": ({"accession": "Q99999"},)},
    )
    monkeypatch.setattr(unknown_molecule_uniprot, "pdb_to_uniprot", {})
    monkeypatch.setattr(
        unknown_molecule_uniprot,
        "uniprot_dict",
        {"P12345": "Embedded UniProt name", "Q99999": "External UniProt name"},
    )
    monkeypatch.setitem(config["network_annotations"], "use_embedded_sifts", True)
    chain = {
        "chain_id": "A",
        "unique_chain_id": "1ABC:A",
        "uniprot_id": "P12345",
        "embedded_uniprot_status": "unique_best_mapping",
    }
    structures = [{"pdb_id": "1ABC", "atom_data": [chain]}]

    unknown_molecule_uniprot.process_molecule_info(
        structures,
        pdb_fasta_headers=headers,
    )

    assert chain["uniprot_id"] == "P12345"
    assert chain["molecule_name"] == expected_name
    assert chain["annotation_warnings"][0]["code"] == "EMBEDDED_EXTERNAL_UNIPROT_MISMATCH"


def test_chain_without_uniprot_segments_does_not_stop_later_chain_annotation() -> None:
    chains = [{"chain_id": "A"}, {"chain_id": "B"}]
    extracted = {
        "by_chain": {
            "A": {"pfam": [{"database": "pfam", "accession": "PF00001"}]},
            "B": {
                "uniprot": [
                    {
                        "database": "uniprot",
                        "accession": "P12345",
                        "best_mapping": True,
                    }
                ]
            },
        }
    }

    apply_embedded_annotations(chains, extracted, use_embedded_sifts=True)

    assert chains[0]["embedded_annotations"]["pfam"][0]["accession"] == "PF00001"
    assert chains[1]["uniprot_id"] == "P12345"
    assert chains[1]["annotation_source"] == "embedded_sifts"


def test_identity_is_content_first_and_extended_canonical(tmp_path: Path) -> None:
    mismatched_name = tmp_path / "2xyz.cif"
    mismatched_name.write_text("data_1abc\n_entry.id 1abc\n", encoding="utf-8")
    identity = get_structure_identity(str(mismatched_name))
    assert identity == StructureIdentity("pdb", "pdb_00001abc", "1abc")
    assert identity.display_id == "1ABC"

    future = tmp_path / "future.cif"
    future.write_text("data_pdb_12345678\n_entry.id pdb_12345678\n", encoding="utf-8")
    future_identity = get_structure_identity(str(future))
    assert future_identity.canonical_id == "pdb_12345678"
    assert future_identity.legacy_id is None
    assert future_identity.display_id == "pdb_12345678"


def test_tooltip_selection_does_not_change_precompute_profiles(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setitem(config, "reference_manifest_id", "profile-test-v1")
    profile_before = precomputed.profile_id()

    monkeypatch.setitem(config["network_annotations"], "tooltip_fields", ["uniprot", "pfam"])

    assert precomputed.profile_id() == profile_before


def test_precompute_keeps_all_embedded_segments_when_tooltip_is_uniprot_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "pdb_00001abc_xyz-enrich.cif"
    source.write_text(_enriched_mmcif(), encoding="utf-8")
    monkeypatch.setitem(config, "reference_manifest_id", "nextgen-test-v1")
    monkeypatch.setitem(config, "structure_model_policy", "first")
    monkeypatch.setitem(config["network_annotations"], "tooltip_fields", ["uniprot"])
    parsed = process_single_file(str(source))
    assert parsed is not None

    from pdb2net import pipeline, reference_data, unknown_molecule_uniprot
    from pdb2net.precomputed import build
    from pdb2net.precomputed.schema import materialize_entry

    monkeypatch.setattr(pipeline, "_validate_required_reference_files", lambda: None)
    monkeypatch.setattr(reference_data, "load_pdb_fasta_headers", lambda *_args: {})
    monkeypatch.setattr(
        pipeline, "_parse_input_files", lambda _paths, **_kwargs: [parsed]
    )
    monkeypatch.setattr(
        unknown_molecule_uniprot, "process_molecule_info", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(pipeline, "_run_blast_annotation", lambda _data: None)
    monkeypatch.setattr(build, "calculate_distances_with_ckdtree", lambda _data: [])

    report = precomputed.precompute_sources(tmp_path / "store", [source])
    assert report["failed"] == 0
    cached = precomputed.load_entry(tmp_path / "store", "1abc")
    structure, _interactions, _references = materialize_entry(cached)
    annotations = structure["atom_data"][0]["embedded_annotations"]
    assert set(annotations) == {"uniprot", "pfam", "cath", "scop2"}


def test_network_annotation_config_rejects_unknown_tooltip_database() -> None:
    with pytest.raises(InputValidationError) as exc_info:
        network_annotation_config({"tooltip_fields": ["interpro"]})
    assert getattr(exc_info.value, "code", None) == "INVALID_NETWORK_ANNOTATIONS_CONFIG"


def _annotation_json_bytes(segments: list[dict[str, object]]) -> int:
    return len(
        json.dumps(
            segments,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    )


def test_raw_annotation_byte_limit_is_independent_of_tooltip_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first = {"database": "pfam", "accession": "PF00001", "pdb_start": 1, "pdb_end": 10}
    second = {"database": "pfam", "accession": "PF00002", "pdb_start": 11, "pdb_end": 20}
    one_segment_limit = _annotation_json_bytes([first])
    monkeypatch.setattr(
        network_annotations, "MAX_ANNOTATION_BYTES_PER_DATABASE", one_segment_limit
    )
    monkeypatch.setattr(network_annotations, "MAX_ANNOTATION_BYTES_PER_NODE", 10_000)
    chain = {"embedded_annotations": {"pfam": [second, first]}}

    one_tooltip = annotation_node_metadata(
        [chain],
        annotation_config={
            "use_embedded_sifts": True,
            "tooltip_fields": ["pfam"],
            "max_tooltip_segments_per_database": 1,
        },
    )
    two_tooltips = annotation_node_metadata(
        [chain],
        annotation_config={
            "use_embedded_sifts": True,
            "tooltip_fields": ["pfam"],
            "max_tooltip_segments_per_database": 2,
        },
    )

    assert one_tooltip["annotation_pfam"] == two_tooltips["annotation_pfam"]
    assert json.loads(one_tooltip["annotation_pfam"]) == [first]
    assert one_tooltip["annotation_pfam_total"] == 2
    assert one_tooltip["annotation_pfam_included"] == 1
    assert one_tooltip["annotation_pfam_truncated"] is True
    assert len(one_tooltip["tooltip_lines"]) == 2  # one segment plus the truncation notice
    assert len(two_tooltips["tooltip_lines"]) == 2  # both complete segment lines

    monkeypatch.setattr(
        network_annotations, "MAX_ANNOTATION_BYTES_PER_DATABASE", one_segment_limit - 1
    )
    below = annotation_node_metadata(
        [chain],
        annotation_config={
            "use_embedded_sifts": True,
            "tooltip_fields": ["pfam"],
            "max_tooltip_segments_per_database": 2,
        },
    )
    assert below["annotation_pfam"] == "[]"
    assert below["annotation_pfam_included"] == 0
    assert below["annotation_pfam_truncated"] is True


def test_raw_annotation_node_budget_accepts_exact_limit_and_truncates_one_byte_below(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pfam = {"database": "pfam", "accession": "PF00001"}
    cath = {"database": "cath", "accession": "1.10.8.10"}
    exact = _annotation_json_bytes([pfam]) + _annotation_json_bytes([cath])
    chain = {"embedded_annotations": {"pfam": [pfam], "cath": [cath]}}
    cfg = {
        "use_embedded_sifts": True,
        "tooltip_fields": ["pfam", "cath"],
        "max_tooltip_segments_per_database": 20,
    }
    monkeypatch.setattr(network_annotations, "MAX_ANNOTATION_BYTES_PER_DATABASE", 10_000)
    monkeypatch.setattr(network_annotations, "MAX_ANNOTATION_BYTES_PER_NODE", exact)
    included = annotation_node_metadata([chain], annotation_config=cfg)
    assert included["annotation_pfam_included"] == 1
    assert included["annotation_cath_included"] == 1
    assert included["annotation_truncated"] is False

    monkeypatch.setattr(network_annotations, "MAX_ANNOTATION_BYTES_PER_NODE", exact - 1)
    truncated = annotation_node_metadata([chain], annotation_config=cfg)
    assert truncated["annotation_pfam_included"] == 1
    assert truncated["annotation_cath_included"] == 0
    assert truncated["annotation_cath_truncated"] is True
    assert truncated["annotation_truncated"] is True

    monkeypatch.setattr(
        network_annotations,
        "MAX_ANNOTATION_BYTES_PER_NODE",
        _annotation_json_bytes([pfam]),
    )
    no_empty_overflow = annotation_node_metadata([chain], annotation_config=cfg)
    raw_values = [
        value
        for key, value in no_empty_overflow.items()
        if key in {"annotation_pfam", "annotation_cath"}
    ]
    assert sum(len(value.encode("utf-8")) for value in raw_values) <= _annotation_json_bytes([pfam])
    assert no_empty_overflow["annotation_pfam_included"] == 1
    assert no_empty_overflow["annotation_cath_included"] == 0
    assert "annotation_cath" not in no_empty_overflow
