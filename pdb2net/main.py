from file_parser import read_files_from_csv
from data_processor import process_structure
from unknown_molecule_uniprot import process_molecule_info
from distances import calculate_distances_with_ckdtree
from cytoscape import create_cytoscape_network

# Define external file paths
PDB_FASTA_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\Neuer Ordner\\pdb_seqres.txt"
UNIPROT_FASTA_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\Neuer Ordner\\uniprot_sprot.fasta"

def main(csv_path):
    """
    Main function to load PDB structures, process molecular information, compute distances,
    and generate a network visualization in Cytoscape.

    Args:
        csv_path (str): Path to the CSV file containing PDB file paths.
    """
    try:
        # 🔹 Step 1: Load PDB files from CSV and parse structures
        print("\n📂 Loading PDB files from CSV...")
        structures = read_files_from_csv(csv_path)
        combined_data = []

        # Process structures and extract atom information
        for structure_data in structures:
            processed_data = process_structure(structure_data)
            combined_data.append(processed_data)

        # 🔹 Step 2: Determine molecule names and types using PDB FASTA and UniProt
        print("\n🔍 Determining molecule names, types, and sequences...")
        process_molecule_info(combined_data)

        # 🔹 Debugging Output: Show processed chains with molecule names and sequences
        print("\n📌 Debugging Output: Processed Chains with Name, Type, and Sequence:")
        for structure in combined_data:
            print(f"\n📄 File: {structure['file_path']} (PDB ID: {structure['pdb_id']})")
            for chain in structure["atom_data"]:
                print(f"  🔹 Chain: {chain['chain_id']}")
                print(f"     🏷 Name: {chain['molecule_name']}")
                print(f"     🔬 Type: {chain['molecule_type']}")
                print(f"     🧬 Sequence: {chain['sequence'][:50]}... (truncated)")

        # 🔹 Step 3: Compute distances between molecular chains
        print("\n📏 Computing atomic distances...")
        results = calculate_distances_with_ckdtree(combined_data)

        # 🔹 Output: Display distance calculation results
        print("\n📊 Distance Calculation Results:")
        for result in results:
            print(f"\n📄 File: {result['file_path']}")
            print(f"  🔗 Chain Pair: {result['chain_a']} - {result['chain_b']}")
            print(f"  ⚛ Cα/NN contacts (<15 Å): {result['ca_nn_count']}")
            print(f"  🔍 All atom contacts (<5 Å): {result['all_atoms_close_count']}")

        # 🔹 Step 4: Create the Cytoscape network visualization
        #print("\n🌐 Creating network in Cytoscape...")
        #create_cytoscape_network(results)

    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    # Define the default path for the CSV file containing PDB file paths
    csv_path = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\PDBFiles\\PathsCSV.csv"
    main(csv_path)
