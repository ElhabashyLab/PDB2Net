# PDB2Net

**PDB2Net** automatically extracts **Protein Interaction Networks (PINs)** from **PDB/mmCIF files** and visualizes them as **Cytoscape networks**.  
It uses **Gemmi** for structure parsing, **SciPy cKDTree** for distance-based interaction detection, and **BLAST+** for **UniProt annotation** of unidentified chains.

## Features

- Automatic parsing of `.pdb`, `.cif`, and `.mmCIF` structures, optionally gzip-compressed
- Content-first structure identity with classic and Extended PDB IDs (`1abc` / `pdb_00001abc`)
- Distance-based chain interaction detection  
- Protein-level and chain-level networks  
- Full UniProt annotation via SIFTS and BLAST+, including embedded SIFTS in enriched NextGen mmCIF
- Configurable UniProt, Pfam, CATH, and SCOP2 segment annotations in existing CX2 tooltips
- Export of chain, protein, and combined networks (CX2 format)  
- Optional versioned per-PDB precompute store for fast repeated assembly
- Machine-readable, configuration-free capability contract for deployment compatibility checks

## Documentation

- [Scientific contact rules](docs/contact-semantics-decision.md) and
  [network projections](docs/network-semantics-decision.md)
- [Architecture](docs/architecture.md) and [development checks](docs/development.md)
- [Reference and search lifecycle](docs/runtime-reference-decisions.md)
- [Precomputed store](docs/precomputed_store.md) and [backend contract](docs/server_backend_usage.md)
- [Validation summary](docs/core-quality-audit.md) and [remaining work](docs/core-quality-plan.md)

## Installation

PDB2Net supports Python 3.11 and 3.12. A base installation is sufficient for
headless CX2 export and does not install Cytoscape integration:

```bash
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install .
```

From a source checkout, `python3 -m pip install -r requirements.txt` installs
the same base scientific dependencies. Desktop users who want live Cytoscape
integration can install the optional extra instead:

```bash
python3 -m pip install ".[cytoscape]"
```

Live mode also requires Cytoscape 3.10.4 or newer. Headless runs use
`--headless` or `open_in_cytoscape: false` and never import `py4cytoscape`.
For development setup, Ruff, tests, and the environment pre-flight check, see
[`docs/development.md`](docs/development.md).

## Reference data

| File                    | Source                                                                                                         | Purpose                          |
| ----------------------- | -------------------------------------------------------------------------------------------------------------- | -------------------------------- |
| `pdb_seqres.txt`        | [https://www.rcsb.org/downloads/fasta](https://www.rcsb.org/downloads/fasta)                                   | PDB single-FASTA (chains)        |
| `pdb_chain_uniprot.tsv` | [https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)           | PDB ⇄ UniProt mapping (SIFTS)    |
| `uniprot_sprot.fasta`   | [https://www.uniprot.org/uniprotkb?query=reviewed:true](https://www.uniprot.org/uniprotkb?query=reviewed:true) | Swiss-Prot for building BLAST DB |


### Setting up BLAST for UniProt matching

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

### Create the BLAST database
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

## Configuration

**PDB2Net** loads configuration in layers — later files override earlier ones:

1. `pdb2net/configs/config.base.json` — shared defaults
2. `pdb2net/configs/config.{windows|linux|darwin}.json` — OS defaults
3. `pdb2net/configs/config.local.json` — machine settings *(git-ignored)*
4. `PDB2NET_CONFIG_FILE` — an explicit JSON configuration file
5. individual `PDB2NET_*` environment overrides

> 🗂️ Paths support `~` and `$VARS` expansion.

Each CLI invocation builds a fresh merged configuration snapshot. `--config`
selects the explicit-file layer without changing `PDB2NET_CONFIG_FILE` for a
later invocation. Invalid types in validated configuration sections and invalid
environment values fail before the scientific pipeline starts, with a concise
configuration error.

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
  "network_annotations": {
    "use_embedded_sifts": true,
    "tooltip_fields": ["uniprot"],
    "max_tooltip_segments_per_database": 20
  },
  "distance_thresholds": { "ca_radius": 12.0, "all_atoms_radius": 5.0 },
  "workers": { "parsing": "auto", "blast_threads": "auto" },
  "keep_last_n_networks": 46,
  "export_detailed_interactions": false
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
| `PDB2NET_DIAMOND_THREADS`   | `diamond.threads`                             |
| `PDB2NET_DIAMOND_ITERATE`   | `diamond.iterate` (`true/false/1/0/yes/no`)   |
| `PDB2NET_DIAMOND_SENSITIVITY` | `diamond.sensitivity`                      |
| `PDB2NET_DIAMOND_BLOCK_SIZE` | `diamond.block_size`                        |
| `PDB2NET_DIAMOND_INDEX_CHUNKS` | `diamond.index_chunks`                    |
| `PDB2NET_DIAMOND_MAX_TARGET_SEQS` | `diamond.max_target_seqs`              |
| `PDB2NET_DIAMOND_BATCH_MAX_SEQUENCES` | `diamond.batch_max_sequences`      |
| `PDB2NET_DIAMOND_BATCH_MAX_FASTA_BYTES` | `diamond.batch_max_fasta_bytes`  |
| `PDB2NET_DIAMOND_TEMP_DIR`  | `diamond.temp_dir`                            |
| `PDB2NET_OPEN_IN_CYTOSCAPE` | `open_in_cytoscape` (`true/false/1/0/yes/no`) |
| `PDB2NET_EXPORT_DETAILED_INTERACTIONS` | `export_detailed_interactions` (`true/false/1/0/yes/no`) |
| `PDB2NET_WORKERS_PARSING`   | `workers.parsing` (`auto` or int)             |
| `PDB2NET_WORKERS_BLAST`     | `workers.blast_threads` (`auto` or int)       |
| `PDB2NET_MAX_INPUT_FILES`   | `resource_limits.max_input_files`             |
| `PDB2NET_MAX_TOTAL_INPUT_BYTES` | `resource_limits.max_total_input_bytes`   |
| `PDB2NET_MAX_SINGLE_INPUT_BYTES` | `resource_limits.max_single_input_bytes` |
| `PDB2NET_MAX_TOTAL_INPUT_EXPANDED_BYTES` | `resource_limits.max_total_input_expanded_bytes` |
| `PDB2NET_MAX_SINGLE_INPUT_EXPANDED_BYTES` | `resource_limits.max_single_input_expanded_bytes` |
| `PDB2NET_MAX_DETAILED_INTERACTION_ROWS` | `resource_limits.max_detailed_interaction_rows` |
| `PDB2NET_MAX_DETAILED_INTERACTION_BYTES` | `resource_limits.max_detailed_interaction_bytes` |
| `PDB2NET_MIN_OUTPUT_FREE_BYTES` | `resource_limits.min_output_free_bytes` |
| `PDB2NET_USE_EMBEDDED_SIFTS` | `network_annotations.use_embedded_sifts` |
| `PDB2NET_ANNOTATION_TOOLTIP_FIELDS` | comma-separated `network_annotations.tooltip_fields` |
| `PDB2NET_MAX_TOOLTIP_SEGMENTS_PER_DATABASE` | `network_annotations.max_tooltip_segments_per_database` |
| `PDB2NET_MAX_PROCESSING_BATCH_BYTES` | `resource_limits.max_processing_batch_bytes` |
| `PDB2NET_COMBINED_MAX_NODES` | `combined_graph_limits.max_nodes`            |
| `PDB2NET_COMBINED_MAX_EDGES` | `combined_graph_limits.max_edges`            |
| `PDB2NET_REFERENCE_MANIFEST_ID` | `reference_manifest_id`                   |
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

For `"all"`, disable both protein-network outputs; they require the first-model
policy. For example:

```json
{
  "structure_model_policy": "all",
  "networks": {
    "chain_per_pdb": true,
    "combined_chain_network": true,
    "protein_per_pdb": false,
    "combined_protein_network": false
  }
}
```

Detailed CSVs then use model-qualified chain IDs such as `1ABC:model2:A` in
`Chain_A` and `Chain_B`; the default `"first"` mode keeps author chain IDs.

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

In JSON configuration, `diamond.enabled` must be the Boolean `true` or `false`;
quoted strings such as `"false"` are rejected as configuration errors. Validation
uses the final merged value, so valid higher-priority file or environment
overrides retain their usual precedence.

For a small format-compatible integration check, the test suite includes an
original-style UniRef90 FASTA fixture at
`tests/fixtures/uniref90/mini_uniref90.fasta`. If DIAMOND is installed, this
test builds a temporary real `.dmnd` database and searches it end to end:

```bash
.venv/bin/python -m pytest -q \
  tests/test_uniprot_matching.py::test_real_diamond_search_against_original_format_mini_uniref90
```

The generated database stays in pytest's temporary directory. The fixture is
only for validating database construction, command invocation, output parsing,
and fallback behavior; it is not scientifically representative and must not be
used as a production UniRef90 database.

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

The legacy batch runner continues with later batches after a child failure or
timeout. It exits with status 1 for failed/timed-out batches or files skipped by
its size/stat checks. Known content exclusions inside an otherwise successful
child run are reported in that run's summaries and do not change its exit code.
Inspect `INPUT_FILE_SKIPPED` warnings when the full submitted dataset is required.

- Output goes to a **timestamped** subfolder in `output_path`, e.g.:
  `/…/Networks/2025-10-20_18-32-45/`. If that name is already reserved,
  PDB2Net creates a distinct suffixed directory such as
  `2025-10-20_18-32-45_2/`.

For backend-style jobs, add `--web-output-dir` to collect stable
user-facing outputs:

```bash
python3 -m pdb2net run \
  --input-dir /path/to/job/inputs \
  --output-dir /path/to/job/work \
  --web-output-dir /path/to/job/outputs \
  --headless
```

The `--web-output-dir` path must either not exist or be an empty regular
directory. PDB2Net reserves it before processing and rejects non-empty targets;
use a fresh path for every run.

See [`docs/server_backend_usage.md`](docs/server_backend_usage.md) for the
worker-facing output contract.

Uploaded/local enriched NextGen mmCIF files are ordinary supported inputs; no
separate mode is required. When `network_annotations.use_embedded_sifts` is
enabled, PDB2Net reads embedded SIFTS segments directly from the same mmCIF
document. `tooltip_fields` controls which of `uniprot`, `pfam`, `cath`, and
`scop2` are added to existing network-node tooltips. UniProt is enabled by
default; the other categories are opt-in. Ambiguous embedded UniProt mappings
remain explicitly unassigned instead of being guessed.

Valid structure files may have arbitrary filenames: an identifier in the file
content takes precedence. Renaming a file does not require editing its internal
IDs. Conflicting IDs inside the document are rejected, as are duplicate
structure identities supplied in the same run.

Single-chain AlphaFold DB structures are also supported as raw PDB/mmCIF input.
An AlphaFold-style filename or retained internal AlphaFold identity can provide
the UniProt accession directly. After renaming, mmCIF may retain that identity;
a PDB file without it may need BLAST fallback. Normal reference-data requirements
still apply. The [default real-file tests](docs/development.md#frozen-real-file-integration-tests)
exercise these cases alongside complete NextGen files.

Direct AlphaFold assignment respects the resolved content identity, existing
SIFTS assignments and explicit mapping ambiguity. A conflicting AlphaFold-style
filename cannot reclassify a nucleic-acid chain or replace another content
identity. Unresolved protein chains continue through sequence-search fallback.

Fallback sequence extraction and residue-based molecule typing use observed
polymer residues, using the roles already assigned by Gemmi. Waters and
non-polymer ligands are excluded from these annotation steps, even when they
share an author chain ID. Modified polymer residues such as MSE remain included
and normalized; unknown polymer residues remain `X` and still count against the
existing sequence-quality limit. This does not reconstruct missing residues.

Displayed `aa_len`/`nt_len` count the observed residue records by their polymer
role and type. Modified or unknown protein/DNA/RNA polymer residues count;
nonpolymers do not, even when named like an amino acid. These lengths do not
apply contact-occupancy filters or infer missing sequence regions. Atom
alternatives do not multiply the count. Callers without usable polymer metadata
retain the legacy name-based interpretation.

Contact calculations use only heavy atoms of protein/DNA/RNA polymer residues,
including the actual atoms of modified residues. Water, separate glycans, lipids,
cofactors, ions and other ligands do not contribute. Explicit occupancy zero is
excluded; positive or missing occupancy is retained without weighting or warnings.
All eligible alternative positions are considered, but each logical atom-site
pair counts only once. Detailed CSVs contain one row per such pair with the
minimum distance (rounded only for display). An insertion code distinguishes
residue sites; an altloc label or alternative residue name does not multiply a
site. These are possible geometric contacts, not binding probabilities or proof
that all contacts coexist in one conformation.

The existing cutoffs remain: protein pairs need at least ten Cα site pairs at
≤12 Å and one heavy-atom site pair at ≤5 Å; pairs involving DNA/RNA need at least
one heavy-atom site pair at ≤5 Å. The default uses the first structure model;
`structure_model_policy: "all"` processes models separately without cross-model
contacts or merging. No AlphaFold confidence filter is applied.

A supported chain without eligible heavy contact atoms is omitted with
`CHAIN_WITHOUT_CONTACT_ATOMS` in the run reports. A chain with atoms and zero
contacts remains valid. Invalid contact occupancy or ambiguous atom identity
is a file-content error; missing occupancy alone is not an error.

Raw-file runs continue past known file-content failures (including structures
with no supported chains). Each omitted file is reported as `INPUT_FILE_SKIPPED`
with its filename, cause and correction advice in the run summaries. The CLI
prints exclusions explicitly. Successful output metadata lists only analyzed
files; `resources.input.files` retains the original input count. A successful
run with exclusions exits with code 0, so automation must also inspect warnings
when it requires the entire submitted dataset. Combined networks contain only
the analyzed subset. If no input is usable, the run fails and retains all file
diagnostics without publishing results. Resource, integrity, tool and unexpected
program failures still abort the run. This adds no separate full parsing pass.

An empty or whitespace-only chain identifier in a selected model produces
`MISSING_CHAIN_ID` and excludes the whole file, even if it also contains named
chains. Correct the chain identifiers or export the structure again; Core does
not invent names. With `structure_model_policy: "first"`, later models do not
participate in this check; with `"all"`, every selected model does. Offline
`precompute` reports the same cause and does not publish a final store manifest
when any input fails.

### Sequence-search timeout and cancellation

`sequence_search_timeout_seconds` defaults to `null` (no additional Core limit).
Set positive whole seconds in the existing config layers or through
`PDB2NET_SEQUENCE_SEARCH_TIMEOUT_SECONDS`; the ENV value `null` disables it,
and an unset/empty ENV variable does not override config. Zero, negative values,
Booleans and fractional seconds are invalid. The limit applies to each BLASTP
or DIAMOND process, including a separate budget for a BLAST format-compatibility
retry. It is not a total-job limit and does not cover the separate makeblastdb
helper. Time to stop/reap a process is additional; OS-level process creation
cannot always be interrupted.

On timeout, pending search chunks do not start and active direct search processes
are stopped, given a shared grace period of at most five seconds, then killed
if necessary and reaped before query files are removed. Search processes remain
in the caller's process group. Arbitrary wrappers that detach their own process
trees are outside this guarantee. A timeout fails the run; it is neither a
negative search result nor a reason for DIAMOND fallback. Scientific cache/store
profiles are unchanged by this runtime setting.

Ctrl+C during sequence searching triggers the same cleanup. The library then
propagates `KeyboardInterrupt`; `run` and `precompute` CLI commands return 130.
Folder-based raw runs retain a failed summary with `RUN_INTERRUPTED` when the
reporting path is available. No success result is published. With `null`, a
standalone search can still run indefinitely unless interrupted; the existing
Web job time limit remains a separate outer bound.

### Reference generations and process lifetime

Use one unchanged, matching set of SIFTS, PDB FASTA, UniProt FASTA and search
databases for the lifetime of a Core process. Prepare updates in a new directory,
assign a new `reference_manifest_id`, and start a new Core process with those
paths. Keep the previous generation available until its running jobs finish.
For Python-library or notebook use, restart the process/kernel when changing
reference generations: `activate_config()` does not reload cached reference
tables. In-place replacement during a running process is unsupported and is
not detected by repeated file scans.

The caches serve different purposes: parsed reference lookups live in RAM;
sequence-search results use a database/policy-specific cache; precomputed entries
belong to a fixed scientific/reference profile. The Web workflow already starts
a new Core process per job and mounts references read-only. No automatic cache
invalidation, reference-file hashing per job or change in memory usage is implied.

### What the four network outputs mean

| Output | Nodes and connections |
| --- | --- |
| Chain per structure | One node per retained chain, with its measured geometric contacts. Keeps contacts between copies of the same protein and contacts of typed protein chains without a UniProt assignment. |
| Protein per structure | Merges chains with the same UniProt ID and sums contacts between different IDs. Includes DNA/RNA chains contacting assigned proteins. Omits contacts within one UniProt group, proteins without an assigned ID, and nucleic-acid-only edges. |
| Combined chain | Keeps separate structure/chain identities and their contacts. Connects components across files through shared UniProt IDs using red dotted `identity` edges. Unlinked structures remain available in per-structure outputs. |
| Combined protein | Summarizes those linked components into UniProt protein nodes, summing contacts between different IDs across the included chains/structures. Contains no DNA/RNA nodes or self-loops. |

For example, if two chains assigned to protein P contact Q with 5 and 7 atom
pairs, the protein edge P–Q has weight 12. A contact between those two P chains
remains visible in the chain network. Summed weights are contact totals; they
do not reconstruct individual chain contributions or measure binding affinity.
The chain network is the detailed view for contacts between copies of the same
protein and for chains without a UniProt assignment. A molecule typed as
`Unknown` is retained as a chain node but is not assigned contact edges.

Cross-file `identity` edges have `all_atoms_count: 0`: they represent a shared
UniProt assignment and carry no geometric contact weight. Their red dotted
style and technical layout attraction are independent of this value. When
analyzing contact topology, filter by `interaction` to exclude identity links;
a zero weight alone does not exclude them from unweighted graph statistics.

PDB2Net analyzes supplied coordinates. It does not generate assembly/symmetry
partners or determine whether a crystal contact is biologically relevant.
Already expanded coordinates can be analyzed when their chains have unique
identifiers within the selected model. `assemble` reads precomputed networks;
it does not build a biological assembly.

With `structure_model_policy: "all"`, per-structure chain networks are exported
separately for each model; no geometric contacts cross model boundaries.
Combined chain views preserve model-qualified chain IDs and may link them
across files by UniProt identity. Both protein outputs are rejected for `all`.
Model indices in outputs are sequential positions, not necessarily the source
file's original MODEL serial numbers. There is no merging of model contact counts.

### Optional precomputed PDB store

The normal laptop/desktop workflow above remains unchanged and does not require
a cache, database server, Docker, or the webserver repository. For repeated
analysis of a local PDB archive, PDB2Net can additionally precompute the compact
graph primitives for each structure:

```bash
export PDB2NET_REFERENCE_MANIFEST_ID="pdb-sifts-swissprot-2026-07"

pdb2net precompute \
  --input-dir /srv/pdb/archive \
  --store /srv/pdb2net/precomputed \
  --config /srv/pdb2net/config.json \
  --recursive \
  --headless
```

The store contains annotated chain metadata and filtered chain-pair interaction
edges for one versioned standard profile. It deliberately contains neither raw
coordinates nor detailed atom-pair CSV data. Assemble one or more cached PDB
entries into the usual CX2 outputs with:

```bash
pdb2net assemble \
  --store /srv/pdb2net/precomputed \
  --pdb-id pdb_00001abc \
  --pdb-id pdb_00002xyz \
  --output-dir ./outputs \
  --headless
```

`assemble` is read-only. An absent entry fails with
`PRECOMPUTED_ENTRY_MISSING`; add it by building a new store with the offline
`precompute` command and promote that store operationally.

Precomputed schema 3 uses canonical Extended PDB IDs and one fixed scientific
profile per store. It requires `structure_model_policy: "first"`, a non-empty
and operationally versioned `reference_manifest_id`, and
`export_detailed_interactions: false`. See
[`docs/precomputed_store.md`](docs/precomputed_store.md) for the artifact layout,
validation, and safe promotion procedure.

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

## Citation

Use the software citation in [`CITATION.cff`](CITATION.cff). Publication details
can be updated there when a corresponding manuscript or archived release is
available.

## Authors
- Gregor Habitzreither
- Hadeer Elhabashy

## Contact
If you have any questions or inquiries, please feel free to contact Hadeer Elhabashy at (Elhabashylab [@] gmail.com))

## License
- The **PDB2NET code** in this repository is licensed under the [MIT License](./LICENSE).
