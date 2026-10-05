"""Independent source/Euclidean oracle for the frozen real structures.

No PDB2Net imports, production atom selectors, KD-trees or exported snapshots.
Gemmi is used only to read source structure records. Small blocked NumPy
distance matrices count logical sites, taking the union of their alternatives.
Run manually with --output-dir outside expected/ to reproduce candidate facts.
"""

import argparse
import csv
import itertools
import json
from pathlib import Path

import gemmi
import numpy as np


DATASET = Path(__file__).parent / "fixtures/core_reference"
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL SEC PYL".split(),
              "ARNDCQEGHILKMFPSTWYVUO"))
PARENTS = {"MSE": "MET", "TPO": "THR", "PTR": "TYR", "SEP": "SER", "CSO": "CYS", "HYP": "PRO"}
POLYMERS = {"PeptideL": "Protein", "PeptideD": "Protein", "CyclicPseudoPeptide": "Protein",
            "Dna": "DNA", "Rna": "RNA", "DnaRnaHybrid": "DNA/RNA"}
STANDARD = set(AA) | {"A", "U", "G", "C", "I", "DA", "DT", "DG", "DC", "DI"}


def fasta_records(path):
    records = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                records.append([line[1:].strip(), ""])
            elif line.strip():
                records[-1][1] += line.strip()
    return records


def reference_annotations(dataset=DATASET):
    mappings = {}
    with (dataset / "references/sifts.tsv").open() as handle:
        for row in csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"):
            mappings.setdefault((row["PDB"], row["CHAIN"]), set()).add(row["SP_PRIMARY"])
    headers = {header.split()[0]: header for header, _ in fasta_records(dataset / "references/pdb_seqres.txt")}
    return mappings, headers


def read_model(path, model_index=0):
    """Return source-derived facts and eligible positions, before annotation."""
    structure = gemmi.read_structure(str(path))
    structure.merge_chain_parts()
    structure.setup_entities()
    kinds = {s: e.polymer_type.name for e in structure.entities for s in e.subchains}
    chains = {}
    for chain in structure[model_index]:
        residues = list(chain)
        polymer = [r for r in residues if r.entity_type == gemmi.EntityType.Polymer
                   and kinds.get(r.subchain) in POLYMERS
                   and any(a.element.name not in {"H", "D"} for a in r)]
        positions = []
        seen = set()
        for residue in polymer:
            for atom in residue:
                if atom.element.name in {"H", "D"} or atom.occ == 0:
                    continue
                site = (residue.seqid.num, residue.seqid.icode.strip(), atom.name)
                coord = tuple(atom.pos)
                identity = site, residue.name, atom.altloc, coord
                if identity in seen:
                    continue
                seen.add(identity)
                name = PARENTS.get(residue.name, residue.name)
                positions.append({"site": site, "coordinate": coord,
                                  "residue": f"{site[0]}{site[1]}:{name}", "atom": atom.name})
        if not positions:
            continue
        facts = {
            "raw_heavy_atoms": sum(a.element.name not in {"H", "D"} for r in residues for a in r),
            "polymer_residues": len(polymer),
            "contact_positions": len(positions),
            "contact_sites": len({p["site"] for p in positions}),
            "zero_occupancy_polymer_atoms": sum(a.element.name not in {"H", "D"} and a.occ == 0
                                                for r in polymer for a in r),
            "modified_residues": [[str(r.seqid), r.name, PARENTS.get(r.name, r.name)]
                                  for r in polymer if r.name not in STANDARD],
            "insertion_codes": [str(r.seqid) for r in polymer if r.seqid.icode.strip()],
        }
        kind = POLYMERS[kinds[polymer[0].subchain]]
        if kind == "Protein":
            facts["protein_sequence"] = "".join(AA.get(PARENTS.get(r.name, r.name), "X") for r in polymer)
        chains[chain.name] = {"facts": facts, "positions": positions, "molecule_type": kind}
    return len(structure), chains


def contact_pairs(left, right, radius, *, ca_only=False):
    """Logical site pairs -> minimum unrounded Euclidean distance.

    Blocking bounds memory. A dictionary unions all atom alternatives for the
    same site pair, independently of occupancy magnitude and input order.
    """
    left = [p for p in left if not ca_only or p["atom"] == "CA"]
    right = [p for p in right if not ca_only or p["atom"] == "CA"]
    if not left or not right:
        return {}
    b = np.asarray([p["coordinate"] for p in right], dtype=float)
    pairs = {}
    for start in range(0, len(left), 128):
        block = left[start:start + 128]
        a = np.asarray([p["coordinate"] for p in block], dtype=float)
        delta = a[:, None, :] - b[None, :, :]
        squared = np.sum(delta * delta, axis=2)
        for i, j in zip(*np.nonzero(squared <= radius * radius)):
            key = block[i]["site"], right[j]["site"]
            distance = float(np.sqrt(squared[i, j]))
            pairs[key] = min(distance, pairs.get(key, float("inf")))
    return pairs


def annotated_facts(path, dataset=DATASET):
    models, chains = read_model(path)
    mappings, headers = reference_annotations(dataset)
    is_af = path.name.startswith("AF-")
    pdbid = path.name.split(".")[0].lower()
    if path.name.startswith("pdb_00008aly"):
        pdbid = "8aly"
    display = "AF-P69905-F1" if is_af else pdbid.upper()
    for chain_id, chain in chains.items():
        header = headers.get(f"{pdbid}_{chain_id}", "")
        # PDB FASTA's explicit mixed-NA description is authoritative in these cases.
        if "mol:na" in header:
            chain["molecule_type"] = next((kind for kind in ["DNA/RNA", "DNA", "RNA"] if kind in header), "Nucleic Acid")
        accessions = sorted(mappings.get((pdbid, chain_id), []))
        source = "sifts" if accessions else None
        if is_af:
            accessions, source = ["P69905"], "direct_uniprot"
        elif path.name.startswith("pdb_00008aly"):
            accessions, source = ["Q9H2K2"], "embedded_sifts"
        elif len(accessions) > 1:
            source = "external_sifts"
        chain["facts"].update(molecule_type=chain["molecule_type"],
                              uniprot_id=accessions[0] if len(accessions) == 1 else None,
                              sifts_accessions=accessions, annotation_source=source)
    return {"display_id": display, "models": models,
            "chains": {name: chain["facts"] for name, chain in chains.items()}}, chains


def expected_contacts(chains, config):
    """All first-model pair counts and acceptance under the frozen policy."""
    rows = []
    for a, b in itertools.combinations(sorted(chains), 2):
        x, y = chains[a], chains[b]
        types = [x["molecule_type"], y["molecule_type"]]
        order = {"Protein": 0, "DNA": 1, "RNA": 2, "DNA/RNA": 3}
        interaction = "-".join(sorted(types, key=order.get))
        ca = len(contact_pairs(x["positions"], y["positions"], config["distance_thresholds"]["ca_radius"], ca_only=True))
        atoms = len(contact_pairs(x["positions"], y["positions"], config["distance_thresholds"]["all_atoms_radius"]))
        filters = config["interaction_filters"]
        if types == ["Protein", "Protein"]:
            accepted = ca >= filters["protein_protein_min_ca_neighbors"] and atoms >= filters["protein_protein_min_all_atom_contacts"]
        else:
            key = "protein_nucleic_acid_min_all_atom_contacts" if "Protein" in types else "nucleic_acid_min_all_atom_contacts"
            accepted = atoms >= filters[key]
        rows.append({"chains": [a, b], "ca_neighbors": ca, "all_atoms_count": atoms,
                     "interaction": interaction, "accepted": accepted})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.resolve() == (DATASET / "expected").resolve():
        parser.error("Write candidates to a separate directory; review before replacing expected data.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    config = json.loads((DATASET / "config.json").read_text())
    files = json.loads((DATASET / "provenance.json").read_text())["structures"]
    structures, networks = {}, {}
    for name in sorted(files):
        if name == "123456.cif":
            continue
        facts, chains = annotated_facts(DATASET / "structures" / name)
        structures[name] = facts
        networks[name] = expected_contacts(chains, config)
        print(name, len(chains), "chains", sum(r["accepted"] for r in networks[name]), "edges", flush=True)
    networks["1a1t.cif#all"] = {
        str(i + 1): expected_contacts(read_model(DATASET / "structures/1a1t.cif", i)[1], config)
        for i in range(structures["1a1t.cif"]["models"])
    }
    for name, value in [("structures.json", structures), ("networks.json", networks)]:
        (args.output_dir / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
