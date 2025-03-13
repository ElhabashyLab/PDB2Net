from cytoscape import create_cytoscape_network
from config_loader import config

def create_protein_network(results, combined_data, run_output_path, network_config):
    """
    Erstellt ein Protein-Level-Interaktionsnetzwerk basierend auf UniProt-IDs.
    Nur Interaktionen mit all_atoms_count > 0 werden berücksichtigt.
    """
    if not network_config["protein_per_pdb"] and not network_config["combined_protein_network"]:
        print("❌ Protein network creation disabled.")
        return

    # 🔹 Map chain IDs to UniProt IDs
    chain_to_uniprot = {}
    for structure in combined_data:
        pdb_id = structure["pdb_id"]
        for chain in structure["atom_data"]:
            chain_id = chain["chain_id"]
            unique_chain_id = f"{pdb_id}:{chain_id}"
            uniprot_id = chain.get("uniprot_id")
            if uniprot_id:
                chain_to_uniprot[unique_chain_id] = uniprot_id

    # 🔹 Combine interactions by UniProt IDs
    protein_interactions = set()  # Nutze ein Set, um doppelte Einträge zu vermeiden
    interaction_data = {}

    for entry in results:
        if entry.get("all_atoms_count", 0) > 0:  # 🔥 Nur mit all_atoms_count > 0
            chain_a = entry["chain_a"]
            chain_b = entry["chain_b"]
            uniprot_a = chain_to_uniprot.get(chain_a)
            uniprot_b = chain_to_uniprot.get(chain_b)

            if not uniprot_a or not uniprot_b or uniprot_a == uniprot_b:
                print(f"⚠️ Skipping interaction: {chain_a} ({uniprot_a}) ↔ {chain_b} ({uniprot_b})")
                continue  # 🔹 Ungültige UniProt-Zuordnungen überspringen

            # 🔹 Unterschiedliche Logik für separate und kombinierte Netzwerke
            if network_config["protein_per_pdb"]:
                pdb_id = chain_a.split(":")[0]  # PDB-ID aus der Kette extrahieren
                interaction_key = (pdb_id, tuple(sorted([uniprot_a, uniprot_b])))  # 🔥 Sortierte Key-Verwendung
            else:
                interaction_key = tuple(sorted([uniprot_a, uniprot_b]))  # 🔥 Kombiniertes Netzwerk ohne PDB-Trennung

            # 🔹 Verhindere doppelte Verbindungen mit einem Set
            if interaction_key not in protein_interactions:
                protein_interactions.add(interaction_key)  # Speichere die Interaktion im Set
                interaction_data[interaction_key] = {
                    "uniprot_a": uniprot_a,
                    "uniprot_b": uniprot_b,
                    "all_atoms_count": entry["all_atoms_count"],  # 🔥 Nur einmal speichern!
                }

    # 🔹 Falls aktiviert: Erstelle separate Protein-Netzwerke pro PDB
    if network_config["protein_per_pdb"]:
        print("\n🌐 Creating separate protein networks for each PDB file...")
        results_by_pdb = {}
        for (pdb_id, (uniprot_a, uniprot_b)), interaction in interaction_data.items():
            if pdb_id not in results_by_pdb:
                results_by_pdb[pdb_id] = []
            results_by_pdb[pdb_id].append(interaction)

        for pdb_id, pdb_results in results_by_pdb.items():
            network_title = f"Protein_Network_{pdb_id}"
            for interaction in pdb_results:
                interaction["chain_a"] = interaction.get("uniprot_a", "UNKNOWN_A")
                interaction["chain_b"] = interaction.get("uniprot_b", "UNKNOWN_B")
            create_cytoscape_network(pdb_results, network_title=network_title, run_output_path=run_output_path)

    # 🔹 Falls aktiviert: Erstelle ein kombiniertes Protein-Netzwerk
    if network_config["combined_protein_network"]:
        print("\n🌐 Creating a single combined protein network...")
        protein_results = []
        for interaction in interaction_data.values():
            protein_results.append({
                "chain_a": interaction.get("uniprot_a", "UNKNOWN_A"),
                "chain_b": interaction.get("uniprot_b", "UNKNOWN_B"),
                "all_atoms_count": interaction["all_atoms_count"],  # 🔥 Fix: Nur all_atoms_count nutzen
            })
        create_cytoscape_network(protein_results, network_title="Combined_Protein_Network", run_output_path=run_output_path)
