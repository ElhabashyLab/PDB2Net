from cytoscape import create_cytoscape_network
from config_loader import config

def create_protein_network(results, combined_data, run_output_path, network_config):
    """
    Creates a protein-level interaction network by combining chains
    based on their UniProt IDs.
    """
    if not network_config["protein_per_pdb"] and not network_config["combined_protein_network"]:
        print("❌ Protein network creation disabled.")
        return

    # 🔹 Map chain IDs to UniProt IDs
    chain_to_uniprot = {}
    for structure in combined_data:
        for chain in structure["atom_data"]:
            pdb_id = structure["pdb_id"]
            chain_id = chain["chain_id"]
            unique_chain_id = f"{pdb_id}:{chain_id}"
            uniprot_id = chain.get("uniprot_id")
            if uniprot_id:
                chain_to_uniprot[unique_chain_id] = uniprot_id

    # 🔹 Combine interactions by UniProt IDs
    protein_interactions = {}
    for entry in results:
        chain_a = entry["chain_a"]
        chain_b = entry["chain_b"]
        uniprot_a = chain_to_uniprot.get(chain_a)
        uniprot_b = chain_to_uniprot.get(chain_b)

        if not uniprot_a or not uniprot_b or uniprot_a == uniprot_b:
            print(f"⚠️ Skipping interaction: {chain_a} ({uniprot_a}) ↔ {chain_b} ({uniprot_b})")
            continue  # 🔹 Überspringt ungültige UniProt-Zuordnungen

        interaction_key = tuple(sorted([uniprot_a, uniprot_b]))
        if interaction_key not in protein_interactions:
            protein_interactions[interaction_key] = {
                "uniprot_a": uniprot_a,
                "uniprot_b": uniprot_b,
                "ca_nn_count": 0,
                "all_atoms_close_count": 0,
            }

        protein_interactions[interaction_key]["ca_nn_count"] += entry["ca_nn_count"]
        protein_interactions[interaction_key]["all_atoms_close_count"] += entry["all_atoms_close_count"]

    # 🔹 Falls aktiviert: Erstelle separate Protein-Netzwerke pro PDB
    if network_config["protein_per_pdb"]:
        print("\n🌐 Creating separate protein networks for each PDB file...")
        results_by_pdb = {}
        for interaction in protein_interactions.values():
            pdb_ids = [
                chain_id.split(":")[0] for chain_id, uni_id in chain_to_uniprot.items()
                if uni_id in [interaction["uniprot_a"], interaction["uniprot_b"]]
            ]
            pdb_ids = list(set(pdb_ids))
            for pdb_id in pdb_ids:
                results_by_pdb.setdefault(pdb_id, []).append(interaction)

        for pdb_id, pdb_results in results_by_pdb.items():
            network_title = f"Protein_Network_{pdb_id}"
            # 🛠 Sicherstellen, dass alle Keys existieren (verhindert KeyErrors)
            for interaction in pdb_results:
                interaction["chain_a"] = interaction.get("uniprot_a", "UNKNOWN_A")
                interaction["chain_b"] = interaction.get("uniprot_b", "UNKNOWN_B")
            create_cytoscape_network(pdb_results, network_title=network_title, run_output_path=run_output_path)

    # 🔹 Falls aktiviert: Erstelle ein kombiniertes Protein-Netzwerk
    if network_config["combined_protein_network"]:
        print("\n🌐 Creating a single combined protein network...")
        protein_results = []
        for interaction in protein_interactions.values():
            protein_results.append({
                "chain_a": interaction.get("uniprot_a", "UNKNOWN_A"),
                "chain_b": interaction.get("uniprot_b", "UNKNOWN_B"),
                "ca_nn_count": interaction["ca_nn_count"],
                "all_atoms_close_count": interaction["all_atoms_close_count"],
            })
        create_cytoscape_network(protein_results, network_title="Combined_Protein_Network", run_output_path=run_output_path)
