# PDB2Net

**PDB2Net** automatically extracts **Protein Interaction Networks (PINs)** from **PDB/mmCIF files** and visualizes them as **Cytoscape networks**.  
It uses **Gemmi** for structure parsing, **SciPy cKDTree** for distance-based interaction detection, and **BLAST+** for **UniProt annotation** of unidentified chains.

## Features

- Automatic parsing of `.pdb`, `.cif`, and `.mmCIF` structures  
- Distance-based chain interaction detection  
- Protein-level and chain-level networks  
- Full UniProt annotation via SIFTS and BLAST+  
- Export of chain, protein, and combined networks (CX2 format)  

---

# System Requirements & Setup

### 1️⃣ Install Python **3.11 or 3.12** 
- **Recommended Version:** Python **3.11**  
- [Download Python](https://www.python.org/downloads/)  
- Ensure that **pip** is installed:

>  python3 -m ensurepip --default-pip


### 2️⃣ Install Required Libraries  
> python3 -m pip install -r requirements.txt

For local development notes and an environment pre-flight check, see [`docs/development.md`](docs/development.md).


### 3️⃣ Install Cytoscape  
- Download **Cytoscape 3.10.4** or newer:  
  [Cytoscape Download](https://cytoscape.org/download.html)  
- **Start once manually**, so it can auto-launch later via PDB2Net.
- On headless servers, Cytoscape is automatically disabled (`open_in_cytoscape = false`).

### 4️⃣ Reference Data (required)

| File                    | Source                                                                                                         | Purpose                          |
| ----------------------- | -------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| `pdb_seqres.txt`        | [https://www.rcsb.org/downloads/fasta](https://www.rcsb.org/downloads/fasta)                                   | PDB single-FASTA (chains)        |
| `pdb_chain_uniprot.tsv` | [https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)           | PDB ⇄ UniProt mapping (SIFTS)    |
| `uniprot_sprot.fasta`   | [https://www.uniprot.org/uniprotkb?query=reviewed:true](https://www.uniprot.org/uniprotkb?query=reviewed:true) | Swiss-Prot for building BLAST DB |


### 5️⃣ Setting up BLAST for UniProt Matching

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

### 6️⃣ **Create the BLAST Database**
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

---

# ⚙️ Configuration (Multi-Layer)

**PDB2Net** loads configuration in layers — later files override earlier ones:

1. `configs/config.base.json` — shared defaults  
2. `configs/config.{windows|linux|darwin}.json` — OS-specific overrides  
3. `configs/config.local.json` — user machine settings *(git-ignored)*  
4. **Environment variables** — highest priority  

> 🗂️ Paths support `~` and `$VARS` expansion.

### Core keys (by file)

`config.base.json`(defaults):
```json
{
  "networks": {
    "chain_per_pdb": true,
    "combined_chain_network": true,
    "protein_per_pdb": true,
    "combined_protein_network": true
  },
  "distance_thresholds": { "ca_radius": 15.0, "all_atoms_radius": 5.0 },
  "workers": { "parsing": "auto", "blast_threads": "auto" },
  "keep_last_n_networks": 46,
  "export_detailed_interactions": true
}
```
### OS examples (adjust to your system):
- `config.windows.json`
```json
{
  "input_folder_path": "%USERPROFILE%/pdb2net/pdb_inputs",
  "pdb_fasta_path": "%USERPROFILE%/pdb2net/reference/pdb_seqres.txt",
  "uniprot_fasta_path": "%USERPROFILE%/pdb2net/reference/uniprot_sprot.fasta",
  "sifts_tsv_path": "%USERPROFILE%/pdb2net/reference/pdb_chain_uniprot.tsv",
  "output_path": "%USERPROFILE%/pdb2net/outputs",
  "cytoscape_path": "C:/Program Files/Cytoscape_v3.10.4/Cytoscape.exe",
  "blast_db_path": "%USERPROFILE%/pdb2net/reference/blast_db",
  "blastp_executable": "blastp",
  "open_in_cytoscape": false
}
```
- `config.linux.json`
```json
{
  "input_folder_path": "$HOME/pdb2net/pdb_inputs",
  "pdb_fasta_path": "$HOME/pdb2net/reference/pdb_seqres.txt",
  "uniprot_fasta_path": "$HOME/pdb2net/reference/uniprot_sprot.fasta",
  "sifts_tsv_path": "$HOME/pdb2net/reference/pdb_chain_uniprot.tsv",
  "output_path": "$HOME/pdb2net/outputs",
  "blast_db_path": "$HOME/pdb2net/reference/blast_db",
  "blastp_executable": "blastp",
  "open_in_cytoscape": false
}
```
- `config.darwin.json` (macOS)
```json 
{
  "input_folder_path": "$HOME/pdb2net/pdb_inputs",
  "pdb_fasta_path": "$HOME/pdb2net/reference/pdb_seqres.txt",
  "uniprot_fasta_path": "$HOME/pdb2net/reference/uniprot_sprot.fasta",
  "sifts_tsv_path": "$HOME/pdb2net/reference/pdb_chain_uniprot.tsv",
  "output_path": "$HOME/pdb2net/outputs",
  "blast_db_path": "$HOME/pdb2net/reference/blast_db",
  "blastp_executable": "blastp",
  "open_in_cytoscape": true,
  "cytoscape_path": "/Applications/Cytoscape.app/Contents/MacOS/Cytoscape"
}
```
## Environment variable overrides
You can override individual settings via ENV:

| ENV var                     | Maps to config key                            |
| --------------------------- | --------------------------------------------- |
| `PDB2NET_INPUT`             | `input_folder_path`                           |
| `PDB2NET_OUTPUT`            | `output_path`                                 |
| `PDB2NET_PDB_FASTA`         | `pdb_fasta_path`                              |
| `PDB2NET_UNIPROT_FASTA`     | `uniprot_fasta_path`                          |
| `PDB2NET_SIFTS_TSV`         | `sifts_tsv_path`                              |
| `PDB2NET_CYTO_PATH`         | `cytoscape_path`                              |
| `PDB2NET_BLAST_DB`          | `blast_db_path`                               |
| `PDB2NET_BLAST_CACHE_PATH`  | `blast_cache_path` (optional SQLite cache)    |
| `PDB2NET_BLASTP`            | `blastp_executable`                           |
| `PDB2NET_DIAMOND_ENABLED`   | `diamond.enabled` (`true/false/1/0/yes/no`)   |
| `PDB2NET_DIAMOND`           | `diamond.executable`                          |
| `PDB2NET_DIAMOND_UNIREF90_DB` | `diamond.uniref90_db_path`                  |
| `PDB2NET_OPEN_IN_CYTOSCAPE` | `open_in_cytoscape` (`true/false/1/0/yes/no`) |
| `PDB2NET_EXPORT_DETAILED_INTERACTIONS` | `export_detailed_interactions` (`true/false/1/0/yes/no`) |
| `PDB2NET_WORKERS_PARSING`   | `workers.parsing` (`auto` or int)             |
| `PDB2NET_WORKERS_BLAST`     | `workers.blast_threads` (`auto` or int)       |
| `PDB2NET_CA_RADIUS`         | `distance_thresholds.ca_radius`               |
| `PDB2NET_ALL_ATOMS_RADIUS`  | `distance_thresholds.all_atoms_radius`        |
| `PDB2NET_PP_MIN_CA_NEIGHBORS` | `interaction_filters.protein_protein_min_ca_neighbors` |
| `PDB2NET_PP_MIN_ALL_ATOM_CONTACTS` | `interaction_filters.protein_protein_min_all_atom_contacts` |
| `PDB2NET_PNA_MIN_ALL_ATOM_CONTACTS` | `interaction_filters.protein_nucleic_acid_min_all_atom_contacts` |
| `PDB2NET_NA_MIN_ALL_ATOM_CONTACTS` | `interaction_filters.nucleic_acid_min_all_atom_contacts` |
| `PDB2NET_STRUCTURE_MODEL_POLICY` | `structure_model_policy` (`first` or `all`) |

## Examples:

**Windows PowerShell:**
```
setx PDB2NET_INPUT "E:\PDB_Files\Dataset"
setx PDB2NET_OUTPUT "E:\Networks"
setx PDB2NET_OPEN_IN_CYTOSCAPE "true"
```
**Linux/macOS:**
```
export PDB2NET_INPUT=~/pdb2net/pdb_inputs
export PDB2NET_OUTPUT=~/pdb2net/outputs
export PDB2NET_OPEN_IN_CYTOSCAPE=false
```

For server or read-only reference-data deployments, set `blast_cache_path` (or
`PDB2NET_BLAST_CACHE_PATH`) to a writable SQLite file outside the BLAST database
directory. If unset, PDB2Net keeps the previous default next to `blast_db_path`.

For automated webserver jobs, keep `open_in_cytoscape: false`, set
`blast_cache_path` to a writable cache directory, and leave
`structure_model_policy: "first"` unless you intentionally want all models from
multi-model structures represented as separate chain nodes.

### Optional DIAMOND/UniRef90 fallback

PDB2Net can optionally run DIAMOND against a local UniRef90 DIAMOND database
after the normal Swiss-Prot BLAST fallback fails to find an accepted hit. This
fallback is disabled by default and does not download or build UniRef90 itself.

When enabled, configure:

```json
{
  "diamond": {
    "enabled": true,
    "executable": "diamond",
    "uniref90_db_path": "/path/to/uniref90.dmnd",
    "max_target_seqs": 50,
    "assign_uniprot_id": "never"
  }
}
```

DIAMOND/UniRef90 hits are recorded with `annotation_source`,
`matched_database`, `matched_id`, `representative_accession`, and
`annotation_confidence`. By default they do not populate `uniprot_id`, keeping
cross-PDB UniProt identity bridges limited to higher-confidence annotation
sources unless you explicitly set `assign_uniprot_id` to `high_confidence` or
`always`.

## Run the Tool

Once all dependencies and reference files are configured, run PDB2Net headlessly
with explicit input and output folders:

```bash
python3 -m pdb2net run \
  --input-dir /path/to/input_structures \
  --output-dir /path/to/pdb2net_outputs \
  --headless
```

If the package is installed, the console command is also available:

```bash
pdb2net run \
  --input-dir /path/to/input_structures \
  --output-dir /path/to/pdb2net_outputs \
  --headless
```

The legacy config-driven entry point remains available:

```bash
python3 -m pdb2net.main
```

- Output goes to a **timestamped** subfolder in `output_path`, e.g.:
""/…/Networks/2025-10-20_18-32-45/"

For backend-style jobs, add `--web-output-dir` to collect stable
user-facing outputs:

```bash
python3 -m pdb2net run \
  --input-dir /path/to/job/inputs \
  --output-dir /path/to/job/work \
  --web-output-dir /path/to/job/outputs \
  --headless
```

See [`docs/server_backend_usage.md`](docs/server_backend_usage.md) for the
worker-facing output contract.

## **User input**  
Valid PDB/mmCIF files found in `input_folder_path`

## Outputs

| File/Folder                 | Description                                                                    |
| --------------------------- | ------------------------------------------------------------------------------ |
| `runtime_analysis.txt`      | Timing summary (parsing, classification, BLAST, interaction, exports)          |
| `manifest.json` / `run_summary.json` | Machine-readable run status, inputs, generated files, counts, config snapshot, contract/package versions, warnings, and errors |
| `outputs/summary.json`      | Stable web-output summary with copied outputs, counts, status, config snapshot, warnings, and errors when `--web-output-dir` is used |
| `*.cx2`                     | Cytoscape networks (Chain/Protein/Combined), portable CX2                      |
| `detailed_interactions.csv` | Per-atom residue/atom distance pairs (if `export_detailed_interactions: true`) |
| `error_in_batch_log/`                     | Batch/runtime logs                                                             |


## Network types

**PDB2Net** generates several network representations:

1. **Chain Interaction Network (per PDB)** — Nodes: chains; Edges: interactions  
2. **Combined Chain Network** — All chains across all PDBs  
3. **Protein Network (per PDB)** — Nodes: UniProt IDs; Edges aggregated over chains  
4. **Combined Protein Network** — UniProt nodes across all PDBs


## Cytoscape Behavior (important)

**Headless / Server** (`open_in_cytoscape: false`)  
→ Only **CX2** files are written (no `.cyjs`).  
→ Deterministic positions and visual mappings are embedded.

**Desktop** (`open_in_cytoscape: true`)  
→ Networks are created in **Cytoscape** via *py4cytoscape* and also exported as **CX2**.

#### **Download the UniProt FASTA File**
The BLAST database will be built from a UniProt FASTA file.

1. **Download the latest UniProt Swiss-Prot database**
   - **Manual Download**: [UniProt Swiss-Prot](https://www.uniprot.org/uniprotkb?query=reviewed:true)

2. **Move the file to the BLAST database folder** (adjust the path if necessary):
   ```bash
   mkdir -p C:/blast_db   # Windows (Git Bash)
   mkdir -p ~/blast_db    # Linux/MacOS
   ```
  
---

# Cite
Habitzreither, G., Gautam, Lupas, A., Elhabashy, H. PDB2Net: Automated extraction of biomolecular Interaction Networks from Three-Dimensional Structures. Manuscript in preparation.

# Authors
- Gregor Habitzreither
- Hadeer Elhabashy

# Contact
If you have any questions or inquiries, please feel free to contact Hadeer Elhabashy at (Elhabashylab [@] gmail.com))

# License
- The **PDB2NET code** in this repository is licensed under the [MIT License](./LICENSE).
