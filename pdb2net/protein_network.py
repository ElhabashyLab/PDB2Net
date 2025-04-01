from cytoscape import create_cytoscape_network


def create_protein_network(results, combined_data, run_output_path, network_config):
    """
    Builds and visualizes protein-level interaction networks using UniProt IDs.

    Depending on the configuration, the function creates:
      - separate networks per PDB file
      - a combined protein-level network

    Args:
        results (list): List of interaction dictionaries (chain pairs and distances).
        combined_data (list): Parsed structure data including atom chains and UniProt annotations.
        run_output_path (str): Output directory for saving Cytoscape files.
        network_config (dict): Configuration flags to control which networks to generate.
    """
    if not network_config["protein_per_pdb"] and not network_config["combined_protein_network"]:
        print("Protein network creation is disabled.")
        return

    chain_to_uniprot = {}
    uniprot_to_pdb_ids = {}
    uniprot_to_name = {}

    # Build mapping from chain IDs to UniProt and names
    for structure in combined_data:
        pdb_id = structure["pdb_id"]
        for chain in structure["atom_data"]:
            chain_id = chain["chain_id"]
            unique_chain_id = f"{pdb_id}:{chain_id}"
            uniprot_id = chain.get("uniprot_id")
            name = chain.get("molecule_name", "Unknown")
            if uniprot_id:
                chain_to_uniprot[unique_chain_id] = uniprot_id
                uniprot_to_pdb_ids.setdefault(uniprot_id, set()).add(pdb_id)
                uniprot_to_name[uniprot_id] = name

    protein_interactions = set()
    interaction_data = {}

    # Group interactions by UniProt-level pairs
    for entry in results:
        if entry.get("all_atoms_count", 0) > 0:
            chain_a = entry["chain_a"]
            chain_b = entry["chain_b"]
            uniprot_a = chain_to_uniprot.get(chain_a)
            uniprot_b = chain_to_uniprot.get(chain_b)

            if not uniprot_a or not uniprot_b or uniprot_a == uniprot_b:
                continue

            if network_config["protein_per_pdb"]:
                pdb_id = chain_a.split(":")[0]
                interaction_key = (pdb_id, tuple(sorted([uniprot_a, uniprot_b])))
            else:
                interaction_key = tuple(sorted([uniprot_a, uniprot_b]))

            if interaction_key not in protein_interactions:
                protein_interactions.add(interaction_key)
                interaction_data[interaction_key] = {
                    "uniprot_a": uniprot_a,
                    "uniprot_b": uniprot_b,
                    "all_atoms_count": entry["all_atoms_count"]
                }

    def get_color_group(uniprot_id):
        """
        Determines the color group label based on the number of PDB files a UniProt ID appears in.
        """
        pdbs = uniprot_to_pdb_ids.get(uniprot_id, set())
        if not pdbs:
            return "Multi"
        return "Multi" if len(pdbs) > 1 else list(pdbs)[0]

    def generate_nodes(interactions):
        """
        Generates Cytoscape nodes from interaction edges.

        Args:
            interactions (list): List of edges with UniProt IDs.

        Returns:
            list: Node dictionaries with ID, color group and name.
        """
        nodes = set()
        for inter in interactions:
            nodes.add(inter["chain_a"])
            nodes.add(inter["chain_b"])

        node_list = [{
            "id": node,
            "label": node,
            "color_group": get_color_group(node),
            "molecule_name": uniprot_to_name.get(node, node)
        } for node in nodes]

        return node_list

    # Generate networks per PDB if enabled
    if network_config["protein_per_pdb"]:
        print("\nCreating separate protein networks for each PDB file...")
        results_by_pdb = {}
        for (pdb_id, (a, b)), inter in interaction_data.items():
            results_by_pdb.setdefault(pdb_id, []).append({
                "chain_a": inter.get("uniprot_a", "UNKNOWN_A"),
                "chain_b": inter.get("uniprot_b", "UNKNOWN_B"),
                "all_atoms_count": inter["all_atoms_count"]
            })

        for pdb_id, pdb_results in results_by_pdb.items():
            nodes = generate_nodes(pdb_results)
            network_title = f"Protein_Network_{pdb_id}"
            create_cytoscape_network(pdb_results, network_title, run_output_path, nodes_data=nodes)

    # Generate one combined protein network if enabled
    if network_config["combined_protein_network"]:
        print("\nCreating a single combined protein network...")
        combined_results = []
        for inter in interaction_data.values():
            combined_results.append({
                "chain_a": inter.get("uniprot_a", "UNKNOWN_A"),
                "chain_b": inter.get("uniprot_b", "UNKNOWN_B"),
                "all_atoms_count": inter["all_atoms_count"]
            })

        nodes = generate_nodes(combined_results)
        create_cytoscape_network(combined_results, "Combined_Protein_Network", run_output_path, nodes_data=nodes)
