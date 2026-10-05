"""Small contact oracles: eligible sites, union of alternatives, nearest witness."""

import copy
import csv
import math
import random

import pytest

from pdb2net.contact_atoms import ContactView, count_contacts, unique_contacts
from pdb2net.detailed_results_exporter import DetailedInteractionBudget
from pdb2net.input_contract import InputValidationError
from pdb2net import distances, pipeline
from pdb2net.detailed_results_exporter import export_detailed_interactions


def residue(number, positions, *, name="ALA", polymer=True, kind="PeptideL", insertion=""):
    return {"residue_number": number, "residue_name": name, "is_polymer": polymer,
            "polymer_type": kind, "insertion_code": insertion,
            "atoms": [{"atom_name": atom, "coordinates": [x, 0., 0.], "altloc": alt,
                       "occupancy": occupancy} for atom, x, alt, occupancy in positions]}


def chain(*residues):
    return {"residues": list(residues)}


def test_selection_keeps_modification_atoms_and_unknown_occupancy_without_mutation():
    source = chain(
        residue(1, [("CA", 0, "", 1), ("SE", 1, "", None)], name="MET"),
        residue(2, [("CA", 2, "", 0), ("P", 3, "", .01)], name="SER"),
        residue(3, [("CA", 4, "", 1)], polymer=False),
        residue(4, [("O", 5, "", 1)], name="HOH", polymer=False),
        residue(5, [("CA", 6, "", 1)], name="GLC", kind="SaccharideD"),
        residue(6, [("P", 7, "", 1)], name="DA", kind="Dna"),
        residue(7, [("P", 8, "", 1)], name="U", kind="Rna"),
    )
    source["residues"][0]["original_residue_name"] = "MSE"
    source["residues"][1]["original_residue_name"] = "SEP"
    before = copy.deepcopy(source)
    assert [(p.residue, p.atom_name) for p in ContactView(source).positions] == [
        ("1:MET", "CA"), ("1:MET", "SE"), ("2:SER", "P"), ("6:DA", "P"), ("7:U", "P")]
    assert len(ContactView(source, "ca").positions) == 1
    assert source == before


@pytest.mark.parametrize("occupancy", [-1, math.nan, math.inf, "wrong"])
def test_invalid_contact_occupancy_is_a_content_error(occupancy):
    with pytest.raises(InputValidationError, match="INVALID_ATOM_OCCUPANCY"):
        ContactView(chain(residue(1, [("CA", 0, "", occupancy)])))
    # Irrelevant ligands do not invalidate a usable polymer.
    assert len(ContactView(chain(residue(1, [("CA", 0, "", occupancy)], polymer=False))).positions) == 0


def test_all_alternatives_count_once_and_csv_witness_is_nearest_independent_of_order():
    a = chain(residue(10, [("CA", 0, "A", .9), ("CA", 4, "B", .1)], insertion="A"))
    b = chain(residue(20, [("CA", 7, "A", .8), ("CA", 6, "B", .2)]))
    for reverse in (False, True):
        left, right = copy.deepcopy(a), copy.deepcopy(b)
        if reverse:
            left["residues"][0]["atoms"].reverse()
            right["residues"][0]["atoms"].reverse()
        av, bv = ContactView(left), ContactView(right)
        assert count_contacts(av, bv, 5) == 1
        (atom_a, atom_b, distance), = unique_contacts(av, bv, 5, budget=DetailedInteractionBudget(max_rows=1))
        assert (atom_a.residue, atom_a.altloc, atom_b.altloc, distance) == ("10A:ALA", "B", "B", 2)


def test_insertion_codes_and_atom_names_are_distinct_sites_but_microheterogeneity_is_not():
    a = chain(residue(1, [("CA", 0, "A", .5)], name="ALA"),
              residue(1, [("CA", 1, "B", .5), ("CB", 2, "B", .5)], name="SER"),
              residue(1, [("CA", 3, "", 1)], insertion="A"))
    b = chain(residue(2, [("CA", 4, "", 1)]))
    av, bv = ContactView(a), ContactView(b)
    assert count_contacts(av, bv, 5) == 3
    rows = list(unique_contacts(av, bv, 5))
    assert [(x.residue, x.atom_name, d) for x, _, d in rows] == [
        ("1:SER", "CA", 3), ("1:SER", "CB", 2), ("1A:ALA", "CA", 1)]


def test_tie_breaker_uses_stable_residue_and_altloc_identity():
    residues = [residue(1, [("CA", 1, "B", .5)], name="SER"),
                residue(1, [("CA", -1, "A", .5)], name="ALA")]
    b = ContactView(chain(residue(2, [("CA", 0, "", 1)])))
    for order in (residues, residues[::-1]):
        (a, _, distance), = unique_contacts(ContactView(chain(*order)), b, 5)
        assert (a.residue, a.altloc, distance) == ("1:ALA", "A", 1)


def test_cutoff_is_inclusive_and_applied_before_rounding():
    a = ContactView(chain(residue(1, [("CA", 0, "A", .5), ("CA", -1, "B", .5)])))
    b = ContactView(chain(residue(2, [("CA", 5, "", 1)]), residue(3, [("CA", 5.001, "", 1)])))
    assert count_contacts(a, b, 5) == 1
    assert [d for _, _, d in unique_contacts(a, b, 5)] == [5]


def test_no_altloc_count_uses_count_neighbors_without_enumerating_pairs():
    a = ContactView(chain(residue(1, [("CA", 0, "", 1)])))
    b = ContactView(chain(residue(2, [("CA", 1, "", 1)])))
    class CountOnly:
        def count_neighbors(self, other, radius):
            assert other is b.tree and radius == 5
            return 1
    a.tree = CountOnly()
    assert count_contacts(a, b, 5) == 1


def test_conflicting_positions_without_alternative_identity_fail():
    with pytest.raises(InputValidationError, match="AMBIGUOUS_ATOM_IDENTITY"):
        ContactView(chain(residue(1, [("CA", 0, "", 1), ("CA", 1, "", 1)])))
    # Exact duplicate records cannot manufacture extra contacts.
    view = ContactView(chain(residue(1, [("CA", 0, "", 1), ("CA", 0, "", 1)])))
    assert len(view.positions) == 1


@pytest.mark.parametrize("sites,expected", [(9, 0), (10, 1)])
def test_alternatives_cannot_manufacture_the_protein_ca_minimum(sites, expected, monkeypatch):
    monkeypatch.setitem(distances.config, "distance_thresholds", {"ca_radius": 12., "all_atoms_radius": 5.})
    monkeypatch.setitem(distances.config, "interaction_filters", {})
    chains = []
    for label, offset in [("A", 0), ("B", 1)]:
        item = chain(*(residue(i, [("CA", i * 30 + offset, "A", .5),
                                  ("CA", i * 30 + offset + .1, "B", .5)]) for i in range(sites)))
        item.update(chain_id=label, unique_chain_id=f"TST:{label}", molecule_type="Protein")
        chains.append(item)
    edges = distances.calculate_distances_with_ckdtree([{"file_path": "test", "atom_data": chains}])
    assert len(edges) == expected
    if edges:
        assert edges[0]["ca_neighbors"] == edges[0]["all_atoms_count"] == sites


def test_csv_and_counts_share_polymer_occupancy_and_altloc_selection(tmp_path):
    a = chain(residue(1, [("P", 0, "A", .7), ("P", 1, "B", .3)], name="DA", kind="Dna"),
              residue(2, [("P", 0, "", 0)], name="DA", kind="Dna"),
              residue(3, [("O", 0, "", 1)], name="HOH", polymer=False))
    b = chain(residue(1, [("P", 3, "A", .8), ("P", 2.004, "B", .2)], name="DA", kind="Dna"))
    for label, item in [("A", a), ("B", b)]:
        item.update(chain_id=label, unique_chain_id=f"TST:{label}", molecule_type="DNA")
    structure = {"pdb_id": "TST", "file_path": "test", "atom_data": [a, b]}
    edges = distances.calculate_distances_with_ckdtree([structure])
    assert edges[0]["all_atoms_count"] == 1
    budget = DetailedInteractionBudget(max_rows=1)
    export_detailed_interactions(structure, edges, str(tmp_path), budget=budget)
    with (tmp_path / "TST_detailed_interactions.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == budget.rows == 1
    assert rows[0]["Distance"] == "1.0" and rows[0]["Residue_A"] == "1:DA"


def test_parser_excludes_empty_contact_chain_but_retains_unknown_occupancy(tmp_path):
    source = tmp_path / "1abc.pdb"
    source.write_text("".join(
        f"ATOM  {i:5d}  CA  ALA {label}{1:4d}    {0.:8.3f}{0.:8.3f}{0.:8.3f}{occupancy} 20.00           C  \n"
        for i, label, occupancy in [(1, "A", "  0.00"), (2, "B", "      ")]))
    parsed = pipeline.process_single_file(str(source))
    assert [c["chain_id"] for c in parsed["atom_data"]] == ["B"]
    assert len(parsed["input_warnings"]) == 1
    assert parsed["input_warnings"][0]["code"] == "CHAIN_WITHOUT_CONTACT_ATOMS"
    assert "A" in parsed["input_warnings"][0]["message"]


@pytest.mark.parametrize("left_alternatives", [False, True])
def test_csv_budget_precheck_runs_before_neighbor_list_with_right_alternatives(left_alternatives):
    positions = [("CA", 0, "A", .5), ("CA", .1, "B", .5)] if left_alternatives else [("CA", 0, "", 1)]
    a = ContactView(chain(residue(1, positions), residue(2, positions)))
    b = ContactView(chain(residue(1, [("CA", 0, "A", .5), ("CA", 20, "B", .5)])))
    tree, calls = a.tree, []
    class ObservedTree:
        def query_ball_point(self, *args, **kwargs):
            calls.append("count" if kwargs.get("return_length") else "materialize")
            return tree.query_ball_point(*args, **kwargs)
    a.tree = ObservedTree()
    with pytest.raises(InputValidationError, match="DETAILED_INTERACTION_ROW_LIMIT_EXCEEDED"):
        list(unique_contacts(a, b, 5, budget=DetailedInteractionBudget(max_rows=1)))
    assert calls == ["count"]


def test_randomized_alternative_union_and_minima_match_exhaustive_oracle():
    rng = random.Random(20261003)
    for _ in range(20):
        left, right = [], []
        for residues in (left, right):
            for site in range(6):
                records = [("CA", rng.uniform(-10, 10), label, .5) for label in "AB"]
                residues.append(residue(site, records))
        for radius in (0., 1., 5.):
            expected = {}
            for a in left:
                for b in right:
                    minimum = min(math.dist(x["coordinates"], y["coordinates"])
                                  for x in a["atoms"] for y in b["atoms"])
                    if minimum <= radius:
                        expected[str(a["residue_number"]), str(b["residue_number"])] = minimum
            # Reverse input ordering and apply one rigid transformation.
            for transformed in (False, True):
                av, bv = copy.deepcopy(left), copy.deepcopy(right)
                if transformed:
                    for records in (av, bv):
                        records.reverse()
                        for res in records:
                            res["atoms"].reverse()
                            for atom in res["atoms"]:
                                x = atom["coordinates"][0]
                                atom["coordinates"] = [4., -x + 7., 2.]
                a, b = ContactView(chain(*av)), ContactView(chain(*bv))
                assert count_contacts(a, b, radius) == len(expected)
                observed = {(x.site[0], y.site[0]): d for x, y, d in unique_contacts(a, b, radius)}
                assert observed == pytest.approx(expected)
