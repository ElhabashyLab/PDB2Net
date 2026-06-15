# Server Backend Usage

PDB2Net can be called by a future webserver worker without coupling the core
package to PHP, MySQL, Apache, or a specific queue system. The web layer should
own uploads, authentication, database rows, scheduling, and downloads. PDB2Net
should only receive filesystem paths and write outputs.

## Recommended Worker Command

```bash
pdb2net run \
  --input-dir /path/to/job/inputs \
  --output-dir /path/to/job/work \
  --web-output-dir /path/to/job/outputs \
  --headless
```

Without an installed console script, use:

```bash
python3 -m pdb2net run \
  --input-dir /path/to/job/inputs \
  --output-dir /path/to/job/work \
  --web-output-dir /path/to/job/outputs \
  --headless
```

## Required Configuration

Provide reference data through `pdb2net/configs/config.local.json`, an explicit
`--config` file, or environment variables:

- `PDB2NET_PDB_FASTA`: `pdb_seqres.txt`
- `PDB2NET_SIFTS_TSV`: `pdb_chain_uniprot.tsv`
- `PDB2NET_UNIPROT_FASTA`: Swiss-Prot FASTA
- `PDB2NET_BLAST_DB`: folder containing `uniprot_db.pin`, `uniprot_db.phr`,
  and `uniprot_db.psq`
- `PDB2NET_BLASTP`: `blastp` executable path or command name
- `PDB2NET_BLAST_CACHE_PATH`: writable SQLite cache path, recommended for
  server deployments

Use `--headless` or `PDB2NET_OPEN_IN_CYTOSCAPE=false` for web jobs. A running
Cytoscape desktop UI should not be required on servers.

## Output Contract

PDB2Net always keeps its normal timestamped internal output under
`--output-dir`, with folders such as:

- `combined/`
- `protein/`
- `chain/`
- `distances/`
- `runtime_analysis.txt`
- `manifest.json`
- `run_summary.json`

When `--web-output-dir` is provided, PDB2Net also creates a stable
user-facing structure:

```text
outputs/
├── summary.json
├── networks/
│   └── *.cx2
└── interactions/
    └── *.csv
```

`summary.json` lists copied network files, interaction tables, counts, runtime
analysis, warnings, and errors. It also includes `output_contract_version`
and `pdb2net_version` so workers can detect the summary structure they are
reading. Future webserver code can rely on this filesystem contract while
keeping database and job status logic outside the PDB2Net core.

## Notes For Workers

- Run each job in an isolated job directory.
- Do not write uploaded files into the PDB2Net repository.
- Do not require Cytoscape for automated jobs.
- Treat non-zero process exit as a failed job, then inspect
  `outputs/summary.json` or the internal `run_summary.json` when present.
- Keep large reference datasets and BLAST databases outside git.
