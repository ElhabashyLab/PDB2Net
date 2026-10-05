"""Preserve raw residue data while using polymer residues for annotation/contacts."""

import copy

import pytest

from pdb2net import data_processor, distances, file_parser, pipeline, uniprot_matcher


# Author chain A contains four polymer residues, three ligands and three waters.
# The ligand ALA deliberately has a standard amino-acid name.
RECORDS = [
    ("ATOM", "ALA", 1, "CA", 0.0, "P", 1),
    ("HETATM", "MSE", 2, "CA", 3.0, "P", 1),
    ("ATOM", "UNK", 3, "CA", 6.0, "P", 1),
    ("ATOM", "GLY", 4, "CA", 9.0, "P", 1),
    ("HETATM", "ALA", 101, "CA", 2.0, "L1", 2),
    ("HETATM", "GLC", 102, "C1", 2.5, "L2", 3),
    ("HETATM", "DA", 103, "P", 3.0, "L3", 4),
    ("HETATM", "HOH", 201, "O", 1.0, "W", 5),
    ("HETATM", "HOH", 202, "O", 2.0, "W", 5),
    ("HETATM", "HOH", 203, "O", 3.0, "W", 5),
]


def _structure_text(extension, records):
    if extension == "pdb":
        lines = ["HEADER".ljust(62) + "1ABC\n"]
    else:
        lines = [
            "data_1ABC\n_entry.id 1ABC\n"
            "loop_\n_entity.id\n_entity.type\n"
            "1 polymer\n2 non-polymer\n3 non-polymer\n4 non-polymer\n5 water\n"
            "_entity_poly.entity_id 1\n_entity_poly.type 'polypeptide(L)'\n"
            "loop_\n_struct_asym.id\n_struct_asym.entity_id\n"
            "P 1\nL1 2\nL2 3\nL3 4\nW 5\nloop_\n"
        ]
        lines.extend(f"_atom_site.{field}\n" for field in [
            "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
            "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
            "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
            "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id",
            "auth_atom_id", "pdbx_PDB_model_num",
        ])
    for serial, (group, residue, number, atom, x, label, entity) in enumerate(records, 1):
        element = atom[0]
        if extension == "pdb":
            if serial == 5:
                lines.append("TER\n")
            lines.append(
                f"{group:6s}{serial:5d} {atom:^4s} {residue:3s} A{number:4d}    "
                f"{x:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00          {element:>2s}  \n"
            )
        else:
            lines.append(
                f"{group} {serial} {element} {atom} . {residue} {label} {entity} "
                f"{number if entity == 1 else '.'} ? {x} 0 0 1.0 20.0 "
                f"{number} {residue} A {atom} 1\n"
            )
    return "".join(lines) + ("END\n" if extension == "pdb" else "#\n")


@pytest.mark.parametrize("extension", ["pdb", "cif"])
def test_parser_retains_raw_residues_but_contacts_only_use_polymer(tmp_path, extension):
    path = tmp_path / f"1abc.{extension}"
    path.write_text(_structure_text(extension, RECORDS))
    (chain,) = pipeline.process_single_file(str(path))["atom_data"]
    residues = chain["residues"]

    assert [res["is_polymer"] for res in residues] == [True] * 4 + [False] * 6
    assert residues[1]["residue_name"] == "MET"
    assert residues[1]["original_residue_name"] == "MSE"
    assert uniprot_matcher.extract_sequence_from_parsed_data(chain) == "AMXG"
    assert uniprot_matcher.classify_molecule_type(chain) == "Protein"
    expected_coordinates = [[row[4], 0.0, 0.0] for row in RECORDS]
    assert [atom["coordinates"] for res in residues for atom in res["atoms"]] == expected_coordinates

    # Q01 uses polymer atoms for contacts without changing the source residues.
    distances.coords_cache.clear()
    try:
        points = distances.extract_coordinates(chain, "all_atoms")
        assert points.tolist() == expected_coordinates[:4]
        expected_contacts = sum(abs(row[4] - 1.0) <= 5.0 for row in RECORDS[:4])
        assert distances.count_nearby_atoms(
            distances.cKDTree(points), distances.cKDTree([[1.0, 0.0, 0.0]]), 5.0
        ) == expected_contacts
    finally:
        distances.coords_cache.clear()


@pytest.mark.parametrize(("polymer", "nonpolymer", "expected"), [
    (["ALA", "MSE"], ["HOH"] * 5 + ["GLC"], "Protein"),
    (["A", "U"], ["DA"] + ["HOH"] * 5, "RNA"),
    (["DA", "DT"], ["U"] + ["HOH"] * 5, "DNA"),
    (["DA", "U"], ["ALA"] * 5, "DNA/RNA"),
    (["ALA", "UNK", "UNK"], ["ALA"] * 5, "Unknown"),
    ([], ["ALA"], "Unknown"),
])
def test_molecule_type_ignores_nonpolymer_residues_including_subtype(polymer, nonpolymer, expected):
    chain = {"residues": (
        [{"residue_name": name, "is_polymer": True} for name in polymer]
        + [{"residue_name": name, "is_polymer": False} for name in nonpolymer]
    )}
    original = copy.deepcopy(chain)
    assert uniprot_matcher.classify_molecule_type(chain) == expected
    assert chain == original


def test_unknown_polymer_residues_remain_x_and_preserve_search_quality_filter(monkeypatch):
    chain = {
        "chain_id": "A", "unique_chain_id": "LOCAL:A", "molecule_type": "Unknown",
        "uniprot_id": None,
        "residues": [
            {"residue_name": "ALA", "is_polymer": True},
            {"residue_name": "UNK", "is_polymer": True},
            {"residue_name": "ALA", "is_polymer": False},
        ],
    }
    assert uniprot_matcher.extract_sequence_from_parsed_data(chain) == "AX"
    monkeypatch.setitem(uniprot_matcher.config, "blast", {"max_x_fraction": 0.2})
    monkeypatch.setitem(uniprot_matcher.config, "diamond", {"enabled": False})

    def no_queries(queries, **_kwargs):
        assert queries == [], "A genuine 50% X polymer sequence must still be excluded"
        return {}

    monkeypatch.setattr(uniprot_matcher, "_run_blastp_swissprot_batch", no_queries)
    uniprot_matcher.parallel_blast_search([{
        "pdb_id": "LOCAL", "file_path": "local.pdb", "atom_data": [chain],
    }], max_workers=1)
    assert chain["molecule_type"] == "Protein"
    assert chain["uniprot_id"] is None


def test_unmarked_internal_residue_records_keep_their_existing_meaning():
    chain = {"residues": [{"resname": name} for name in ["MSE", "UNK", "ALA"]]}
    assert uniprot_matcher.extract_sequence_from_parsed_data(chain) == "MXA"
    assert uniprot_matcher.classify_molecule_type(chain) == "Protein"


def test_nonpolymer_nucleotide_with_empty_query_does_not_require_blast(tmp_path, monkeypatch):
    source = tmp_path / "1abc.cif"
    source.write_text(_structure_text("cif", [RECORDS[6]]))
    assert pipeline.process_single_file(str(source))["atom_data"] == []
    # Keep the F02 empty-query guard covered independently of Q01 chain removal.
    parsed = file_parser.parse_structure_input(str(source))
    structure = data_processor.process_structure({**parsed, "file_path": str(source), "pdb_id": "1ABC"})
    (chain,) = structure["atom_data"]
    assert chain["residues"][0]["residue_name"] == "DA"
    assert chain["residues"][0]["is_polymer"] is False
    chain["molecule_type"] = "Unknown"  # No reference-based molecule assignment.
    assert uniprot_matcher.extract_sequence_from_parsed_data(chain) == ""
    assert uniprot_matcher.classify_molecule_type(chain) == "Unknown"
    monkeypatch.setitem(pipeline.config, "blastp_executable", str(tmp_path / "missing-blastp"))
    monkeypatch.setitem(pipeline.config, "diamond", {"enabled": False})

    assert pipeline._blast_fallback_needed([structure]) is False
    pipeline._validate_blast_ready_if_needed([structure])

    def unexpected_process(*_args, **_kwargs):
        raise AssertionError("An empty polymer sequence must not start a search process")

    monkeypatch.setattr(uniprot_matcher.subprocess, "run", unexpected_process)
    uniprot_matcher.parallel_blast_search([structure], max_workers=1)
    assert chain.get("uniprot_id") is None
