"""Reproduce the approved Q01 rule alongside the historical decision examples.

Run from the repository root with:
  .venv/bin/python docs/examples/q01_contact_examples.py

Only tiny synthetic input files in a temporary directory are written. The
production parser/contact code is called directly, without annotation searches
or network export. Alternative counts use independently enumerated atom pairs.
"""

from collections import defaultdict
from dataclasses import dataclass, replace
import json
import math
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import gemmi

from pdb2net import config_loader, distances, pipeline
from pdb2net.uniprot_matcher import classify_molecule_type


@dataclass(frozen=True)
class Atom:
    chain: str
    number: int
    residue: str
    name: str
    x: float
    role: str = "polymer"
    alt: str = ""
    occ: float = 1.0

    @property
    def site(self):
        return self.chain, self.number, self.name


def pdb_text(rows, *, blank_occupancy=False):
    lines = ["HEADER".ljust(62) + "1ABC\n"]
    serial = 0
    for chain in sorted({a.chain for a in rows}):
        atoms = [a for a in rows if a.chain == chain]
        for polymer in (True, False):
            for atom in atoms:
                if (atom.role == "polymer") != polymer:
                    continue
                serial += 1
                record = "ATOM" if polymer and atom.residue != "MSE" else "HETATM"
                occupancy = "      " if blank_occupancy else f"{atom.occ:6.2f}"
                lines.append(
                    f"{record:6s}{serial:5d} {atom.name:^4s}{atom.alt or ' '}{atom.residue:3s} "
                    f"{chain}{atom.number:4d}    {atom.x:8.3f}{0.:8.3f}{0.:8.3f}"
                    f"{occupancy} 20.00          {atom.name[0]:>2s}  \n"
                )
            if polymer:
                lines.append("TER\n")
    return "".join(lines) + "END\n"


def cif_text(rows, *, unknown_occupancy=None):
    # Entity/asym metadata is written independently of the PDB TER convention.
    entities = {key: i + 1 for i, key in enumerate(dict.fromkeys((r.chain, r.role) for r in rows))}
    lines = ["data_1abc\n_entry.id 1abc\nloop_\n_entity.id\n_entity.type\n"]
    lines.extend(f"{value} {'non-polymer' if key[1] == 'ligand' else key[1]}\n" for key, value in entities.items())
    lines.append("loop_\n_entity_poly.entity_id\n_entity_poly.type\n")
    for (chain, role), entity in entities.items():
        if role == "polymer":
            kind = "polyribonucleotide" if any(a.chain == chain and a.residue == "U" for a in rows) else "polypeptide(L)"
            lines.append(f"{entity} '{kind}'\n")
    lines.append("loop_\n_struct_asym.id\n_struct_asym.entity_id\n")
    lines.extend(f"{chain}_{role} {entity}\n" for (chain, role), entity in entities.items())
    lines.append("loop_\n")
    fields = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
              "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code", "Cartn_x",
              "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv", "auth_seq_id", "auth_comp_id",
              "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"]
    lines.extend(f"_atom_site.{field}\n" for field in fields)
    for serial, a in enumerate(rows, 1):
        polymer = a.role == "polymer"
        group = "ATOM" if polymer and a.residue != "MSE" else "HETATM"
        lines.append(
            f"{group} {serial} {a.name[0]} {a.name} {a.alt or '.'} {a.residue} {a.chain}_{a.role} "
            f"{entities[a.chain, a.role]} {a.number if polymer else '.'} ? {a.x} 0 0 "
            f"{unknown_occupancy if unknown_occupancy is not None else a.occ} 20 "
            f"{a.number} {a.residue} {a.chain} {a.name} 1\n"
        )
    return "".join(lines)


def pairs(rows, radius, *, ca=False):
    return [(a, b) for a in rows if a.chain == "A" for b in rows if b.chain == "B"
            if (not ca or a.name == b.name == "CA") and abs(a.x - b.x) <= radius]


def representative(rows):
    """Illustrative residue-label selection, only for these controlled examples."""
    grouped = defaultdict(list)
    for atom in rows:
        grouped[atom.chain, atom.number].append(atom)
    selected = []
    for atoms in grouped.values():
        labels = {a.alt for a in atoms if a.alt}
        pick = min(labels, key=lambda label: (-sum(a.occ for a in atoms if a.alt == label)
                                             / sum(a.alt == label for a in atoms), label)) if labels else ""
        selected.extend(a for a in atoms if not a.alt or a.alt == pick)
    return selected


def count_result(rows, *, protein_pair, union=False):
    heavy, ca = pairs(rows, 5), pairs(rows, 12, ca=True)

    def count(ps):
        return len({(a.site, b.site) for a, b in ps}) if union else len(ps)

    heavy_count, ca_count = count(heavy), count(ca)
    return {"contacts_5A": heavy_count, "ca_pairs_12A": ca_count,
            "edge": heavy_count >= 1 and (not protein_pair or ca_count >= 10)}


def examples():
    far = Atom("A", 1, "ALA", "CA", 20)
    partner = Atom("B", 1, "U", "P", 3)
    cases = {
        "water": [replace(far, x=0), replace(partner, x=20), Atom("A", 101, "HOH", "O", 17, "water")],
        "ligand": [replace(far, x=0), replace(partner, x=20), Atom("A", 101, "GLC", "C1", 17, "ligand")],
        "zero_occupancy": [far, Atom("A", 1, "ALA", "CB", 0, occ=0), partner],
        "fractional_occupancy": [far, Atom("A", 1, "ALA", "CB", 0, occ=0.2), partner],
        "two_altloc_contacts": [replace(far, x=0, alt="A", occ=0.5), replace(far, x=1, alt="B", occ=0.5), partner],
        "minor_altloc_only": [replace(far, alt="A", occ=0.7), replace(far, x=0, alt="B", occ=0.3), partner],
        "modified_polymer": [far, Atom("A", 2, "MSE", "CA", 0), partner],
        "correlated_alternatives": [Atom(chain, 1, "U", "P", x + offset, alt=label, occ=0.5)
                                    for chain, offset in [("A", 0), ("B", 3)]
                                    for label, x in [("A", 0), ("B", 20)]],
    }
    protein_pair = [Atom(chain, i + 1, "ALA", "CA", 40.0 * i + offset)
                    for chain, offset in [("A", 0), ("B", 3)] for i in range(9)]
    protein_pair[0] = replace(protein_pair[0], alt="A", occ=0.5)
    protein_pair.insert(1, replace(protein_pair[0], x=1, alt="B"))
    cases["protein_threshold"] = protein_pair
    return cases


def main():
    root = Path(__file__).resolve().parents[2]
    cfg = json.loads((root / "pdb2net/configs/config.base.json").read_text())
    config_loader.activate_config(cfg)
    output = {"gemmi": gemmi.__version__, "cases": {}, "missing_occupancy": {}}
    # Hand-counted triples: raw location pairs / representative / unique sites.
    expected_counts = {
        "water": (1, 0, 0), "ligand": (1, 0, 0), "zero_occupancy": (1, 0, 0),
        "fractional_occupancy": (1, 1, 1), "two_altloc_contacts": (2, 1, 1),
        "minor_altloc_only": (1, 0, 1), "modified_polymer": (1, 1, 1),
        "correlated_alternatives": (2, 1, 1), "protein_threshold": (10, 9, 9),
    }
    with tempfile.TemporaryDirectory(prefix="pdb2net-q01-examples-") as temp:
        for name, rows in examples().items():
            protein_pair = name == "protein_threshold"
            positive_polymer = [a for a in rows if a.role == "polymer" and a.occ > 0]
            selected = representative(positive_polymer)
            alternatives = {
                "legacy_oracle": count_result(rows, protein_pair=protein_pair),
                "polymer_positive": count_result(positive_polymer, protein_pair=protein_pair),
                "representative": count_result(selected, protein_pair=protein_pair),
                "unique_site_union": count_result(positive_polymer, protein_pair=protein_pair, union=True),
                "occupancy_product_score": round(sum(a.occ * b.occ for a, b in pairs(positive_polymer, 5)), 6),
            }
            assert tuple(alternatives[key]["contacts_5A"] for key in
                         ("legacy_oracle", "representative", "unique_site_union")) == expected_counts[name]
            observed = {}
            for extension, text in [("pdb", pdb_text(rows)), ("cif", cif_text(rows))]:
                path = Path(temp) / f"{name}.{extension}"
                path.write_text(text)
                parsed = pipeline.process_single_file(str(path))
                # Isolate contact semantics from reference lookup and BLAST.
                for chain in parsed["atom_data"]:
                    chain["molecule_type"] = classify_molecule_type(chain)
                actual_atoms = [(c["chain_id"], r["residue_number"], a["atom_name"], a["coordinates"][0], r["is_polymer"])
                                for c in parsed["atom_data"] for r in c["residues"] for a in r["atoms"]]
                assert sorted(actual_atoms) == sorted((a.chain, a.number, a.name, a.x, a.role == "polymer") for a in rows)
                distances.coords_cache.clear()
                distances.tree_cache.clear()
                edges = distances.calculate_distances_with_ckdtree([parsed])
                expected = alternatives["unique_site_union"]
                assert bool(edges) == expected["edge"]
                if edges:
                    assert edges[0]["all_atoms_count"] == expected["contacts_5A"]
                    if protein_pair:
                        assert edges[0]["ca_neighbors"] == expected["ca_pairs_12A"]
                observed[extension] = {"edge": bool(edges), "edge_contacts": edges[0]["all_atoms_count"] if edges else None}
            output["cases"][name] = {**alternatives, "core": observed}

        one = [Atom("A", 1, "ALA", "CA", 0)]
        for label, text, reader in [
            ("pdb_blank", pdb_text(one, blank_occupancy=True), gemmi.read_pdb_string),
            ("pdb_zero", pdb_text([replace(one[0], occ=0)]), gemmi.read_pdb_string),
            ("cif_dot", cif_text(one, unknown_occupancy="."), lambda t: gemmi.make_structure_from_block(gemmi.cif.read_string(t)[0])),
            ("cif_question", cif_text(one, unknown_occupancy="?"), lambda t: gemmi.make_structure_from_block(gemmi.cif.read_string(t)[0])),
        ]:
            occ = reader(text)[0][0][0][0].occ
            output["missing_occupancy"][label] = occ if math.isfinite(occ) else "nonfinite"
    distances.coords_cache.clear()
    distances.tree_cache.clear()
    print(json.dumps(output, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
