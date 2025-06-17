# PDB2Net

PDB2Net is a tool extracts **Protein Interaction Networks (PINs)** from **PDB/mmCIF files** and visualizes them in **Cytoscape**.  
It utilizes **Biopython** for file parsing, **cKDTree** for efficient distance calculations, and **BLAST** for UniProt matching to also process files without PDB-ID.

## System Requirements & Setup

### 1️⃣ Install Python **3.11 or 3.12** 
- **Recommended Version:** Python **3.11**  
- [Download Python](https://www.python.org/downloads/)  
- Ensure that **pip** is installed:

>  python -m ensurepip --default-pip


### 2️⃣ Install Required Libraries  
> pip install -r requirements.txt


### 3️⃣ Install Cytoscape  
- Download **Cytoscape 3.10.3**:  
  [Cytoscape Download](https://cytoscape.org/download.html)  
- **Start Cytoscape manually** once before running the tool. After that it will start automatically when running the tool.
  ### 4️⃣ Download pdb_chain_uniprot.tsv
https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html

### 5️⃣ Download the PDB Single FASTA File (pdb_seqres.txt)
https://www.rcsb.org/downloads/fasta

### 6️⃣ Setting up BLAST for UniProt Matching

#### **Download & Install BLAST+**
1. **Go to the NCBI BLAST+ Download page:**  
   🔗 [NCBI BLAST+ Download](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)  
2. **Download the correct version for your OS:**
   - **Windows**: Download `ncbi-blast-*-win64.exe`
   - **Linux**: Download `ncbi-blast-*-x64-linux.tar.gz`
   - **MacOS**: Download `ncbi-blast-*-universal-macosx.tar.gz`
3. **Install BLAST+**:
   - **Windows**: Run the `.exe` file and follow the installation wizard.
   - **Linux/MacOS**: Extract the files and move them to `/usr/local/bin`:
     ```bash
     tar -xvzf ncbi-blast-*-x64-linux.tar.gz
     sudo mv ncbi-blast-* /usr/local/bin
     ```

## Run the Tool  
Once all dependencies are installed, you can run the tool with:
> python main.py

### ** Output Files**  
After running the tool, the following files will be generated:  
- **Interaction Data (CSV)**
  - Contains all atomic interactions between chains.
- **Cytoscape Networks**
  - **Chain Interaction Networks** → Networks where each chain is a node.
  - **Protein-Level Networks** → Networks based on UniProt IDs.
  - **Combined and Per-PDB Networks** (if enabled in `config.json`).


### **User input**  
A csv file (list) of PDB and CIF file paths in a structured table:

| file_path |
|--------|
| C:\Users\...\101m.pdb  |
| C:\Users\...\9jr2.cif  |
| ...  |






#### **Download the UniProt FASTA File**
The BLAST database will be built from a UniProt FASTA file.

1. **Download the latest UniProt Swiss-Prot database**
   - **Manual Download**: [UniProt Swiss-Prot](https://www.uniprot.org/uniprotkb?query=reviewed:true)

2. **Move the file to the BLAST database folder** (adjust the path if necessary):
   ```bash
   mkdir -p C:/blast_db   # Windows (Git Bash)
   mkdir -p ~/blast_db    # Linux/MacOS
   ```

#### **Create the BLAST Database**
Now, generate the BLAST database from the downloaded UniProt FASTA file.

1. Open a **terminal** (Linux/Mac) or **PowerShell/Git Bash** (Windows).
2. Run the following command:
   ```bash
   makeblastdb -in C:/blast_db/uniprot_sprot.fasta -dbtype prot -out C:/blast_db/uniprot_db
   ```
   Explanation:
   - `-in` → Input FASTA file.
   - `-dbtype prot` → Specifies a **protein** database.
   - `-out` → Output database name (`uniprot_db`).

3. Expected output:
   ```
   Building a new DB, current time: 03/16/2025 12:45:32
   New DB name:   C:/blast_db/uniprot_db
   Number of sequences: 570,000
   ```
   This confirms that BLAST has successfully created the database.


#### **Summary**
✔ **BLAST+ installed**  
✔ **UniProt FASTA downloaded**  
✔ **BLAST database created**  
✔ **Paths configured in `config.json`**  
✔ **BLAST verified with a test query**  

---

## ⚙️ Configuration (`config.json`)  
Before using the tool, adjust the following paths and parameters in `config.json`:

```json
{
    "input_csv_path": "C:/Path/to/CSV-file.csv",
    "pdb_fasta_path": "C:/Path/to/pdb_seqres.txt",
    "uniprot_fasta_path": "C:/Path/to/uniprot_sprot.fasta",
    "sifts_tsv_path": "C:/Path/to/pdb_chain_uniprot.tsv",
    "output_path": "C:/Path/to/output-folder",
    "cytoscape_path": "C:/Program Files/Cytoscape_v3.10.3/Cytoscape.exe",
    "blast_db_path": "C:/blast_db",
    "blastp_executable": "C:/Program Files/NCBI/blast-2.16.0+/bin/blastp.exe",

    "networks": {
        "chain_per_pdb": true,
        "combined_chain_network": true,
        "protein_per_pdb": true,
        "combined_protein_network": true
    },

    "keep_last_n_networks": 30,
    "export_detailed_interactions": true,

    "distance_thresholds": {
        "ca_radius": 15.0,
        "all_atoms_radius": 5.0
    }
}
```
⚠ **Important:** Adjust paths according to your system!





