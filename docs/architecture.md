# Architecture

PDB2Net is a batch pipeline for turning PDB/mmCIF structure files into protein and chain interaction networks.

## Pipeline

1. Discover valid `.pdb`, `.cif`, and `.mmcif` files and validate configured input limits.
2. Group files into optional raw-size-bounded processing batches.
3. Parse each batch with Gemmi and extract chain, residue, and atom data.
4. Annotate chains using SIFTS and PDB/UniProt FASTA reference files.
5. Run BLAST/DIAMOND fallback for unresolved protein chains.
6. Detect interactions with SciPy `cKDTree` distance searches and optionally write detailed tables.
7. Compact completed structures to chain metadata and polymer lengths, releasing atom coordinates.
8. Export per-PDB and combined chain/protein networks as CX2 from the accumulated summaries.
9. Optionally create networks in a running Cytoscape session.

## Main Modules

- `pdb2net/config_loader.py`: layered configuration, path expansion, environment overrides, and headless defaults.
- `pdb2net/file_parser.py`: valid input detection, PDB ID extraction, and Gemmi structure parsing.
- `pdb2net/reference_data.py`: cached loaders for PDB FASTA, SIFTS, and UniProt reference files.
- `pdb2net/data_processor.py`: conversion from parsed structures to chain and atom dictionaries used by later stages.
- `pdb2net/unknown_molecule_uniprot.py`: SIFTS, PDB FASTA, and UniProt FASTA annotation logic for known and unknown molecules.
- `pdb2net/uniprot_matcher.py`: BLAST database creation, BLAST lookup, diagnostics, and SQLite-backed BLAST result caching.
- `pdb2net/distances.py`: chain interaction detection using configurable CA and all-atom distance thresholds.
- `pdb2net/detailed_results_exporter.py`: optional atom-level distance CSV export.
- `pdb2net/visual_style.py`: shared visual profiles, color maps, and linked-identity border annotation.
- `pdb2net/cx2_export.py`: headless CX2 JSON generation.
- `pdb2net/cytoscape.py`: live Cytoscape/py4cytoscape integration and the public network export entry point.
- `pdb2net/protein_network.py`: aggregation from chain interactions to protein-level networks.
- `pdb2net/outputs.py`: per-run output folders and the versioned run/web output contract.
- `pdb2net/batching.py`: streamed batch sizing, timeout handling, and skipped-batch logging.
- `pdb2net/pipeline.py`: single-run orchestration, worker selection, parsing, annotation, interaction detection, and export sequencing.
- `pdb2net/precomputed_store.py`: optional profile-namespaced, validated per-PDB
  chain/edge artifacts plus targeted cache population and cached network assembly.
- `pdb2net/main.py`: stable public entry points for single-run and batch execution.

## External Inputs

The full pipeline depends on reference files that are intentionally not committed:

- `pdb_seqres.txt`: PDB chain FASTA reference.
- `pdb_chain_uniprot.tsv`: SIFTS PDB-to-UniProt mapping.
- `uniprot_sprot.fasta`: Swiss-Prot FASTA used to build the BLAST database.
- BLAST database files generated from the Swiss-Prot FASTA.

## Execution Modes

In headless mode, PDB2Net writes CX2 files directly and does not require Cytoscape. This is the preferred mode for automated checks and server runs.

In desktop mode, PDB2Net connects to Cytoscape through py4cytoscape, creates networks in the UI, applies styles, and also exports CX2.

The optional precomputed mode changes only the input adapter: `precompute`
persists compact annotated chains and standard-profile edges, and `assemble`
feeds validated entries into the same network exporters. It does not persist
raw coordinates or detailed atom-pair CSVs. Profile namespaces prevent silent
mixing across scientific/reference policies, while store ownership, daily mirror
promotion, obsolete-ID retention, and rollback remain deployment concerns.
The ordinary local `run --input-dir` path does not require this store.

## Configuration

Configuration is loaded in layers:

1. `pdb2net/configs/config.base.json`
2. `pdb2net/configs/config.{linux|windows|darwin}.json`
3. `pdb2net/configs/config.local.json`
4. `PDB2NET_CONFIG_FILE`
5. environment variable overrides

`config.local.json` is intended for machine-specific paths and is ignored by Git.

## Testing Strategy

Automated tests should start with small fixtures and headless behavior:

- config merge and environment override behavior
- file extension and PDB ID parsing
- distance calculations on small synthetic structures
- headless CX2 export shape and required fields

Full BLAST, full reference datasets, and Cytoscape UI behavior should be treated as integration checks rather than default test requirements.
