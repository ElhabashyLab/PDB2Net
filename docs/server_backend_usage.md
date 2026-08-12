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

For production worker images, install PDB2Net as a Python package and use the
`pdb2net` console script. Do not rely on mounting `../PDB2Net` into the worker
and running from the source tree.

The module form is useful for local development or transitional deployments
where the console script is not installed:

```bash
python3 -m pdb2net run \
  --input-dir /path/to/job/inputs \
  --output-dir /path/to/job/work \
  --web-output-dir /path/to/job/outputs \
  --headless
```

Use Python 3.11 or 3.12 in the worker image. The project is developed and
tested against those versions; avoid Python 3.13 unless the dependency stack is
explicitly verified there.

## Capability And Compatibility Probe

Before a real worker accepts jobs, it should run the configuration-free probe:

```bash
pdb2net --version
pdb2net capabilities --json
```

The JSON document is the stable handshake between Core and an external worker.
It reports the PDB2Net package version, capability schema, CLI contract, web
output contract, precomputed-artifact schema, web-config schema, available
commands, input formats and relevant feature allowlists. The probe does not
load machine-local configuration or inspect reference data.

Capability schema 2 also exposes the stable `server_interface` object. Its
`commands`, `contracts`, and `scientific_profiles` sections describe the exact
CLI options used by workers, generated config and output/artifact shapes, and
the scientific constraints of upload and precomputed execution. Precomputed
distance thresholds and interaction filters are content-addressed profile
inputs, not universal fixed values. Consumers should compare the complete
`server_interface` object and fail closed on missing, changed, or unknown
fields; generic feature lists outside that object remain additive.
The reported annotation-pipeline ID covers the SIFTS/FASTA/BLAST annotation
flow, including the optional DIAMOND/UniRef90 fallback and its guarded UniProt
assignment policy; changing any of those semantics requires a new ID.

An integration repository should keep its expected values in a tracked
compatibility file, install Core from one immutable full commit SHA, and refuse
to start a real worker when any reported value differs. A deployment-specific
environment file may contain paths and secrets, but must not override the
tracked compatibility keys. Changes to fields or semantics require a deliberate
contract/schema version change and corresponding consumer update.

## Package Installation Requirements

The installed package includes PDB2Net's shared default configuration files:

- `pdb2net/configs/config.base.json`
- `pdb2net/configs/config.linux.json`
- `pdb2net/configs/config.windows.json`
- `pdb2net/configs/config.darwin.json`
- `pdb2net/configs/config.local.example.json`

It intentionally does not include `pdb2net/configs/config.local.json`, because
that file is for machine-specific paths and may contain local reference data
locations. Server deployments should provide reference data paths through
environment variables or a per-job `--config` file.

Before switching a worker image to `PDB2NET_COMMAND=pdb2net`, verify the
installed package from outside the source tree:

```bash
python3 -m pip install .
cd /tmp
python3 -c 'import pdb2net, pathlib; root = pathlib.Path(pdb2net.__file__).parent; print((root / "configs/config.base.json").exists()); print((root / "configs/config.local.json").exists())'
pdb2net --help
pdb2net capabilities --json
```

Expected output for the config check is `True` for `config.base.json` and
`False` for `config.local.json`.

## Required Configuration

For web workers, provide reference data through an explicit `--config` file or
environment variables. `pdb2net/configs/config.local.json` is acceptable for a
local workstation checkout, but should not be used as the production deployment
mechanism for installed worker packages.

- `PDB2NET_PDB_FASTA`: `pdb_seqres.txt`
- `PDB2NET_SIFTS_TSV`: `pdb_chain_uniprot.tsv`
- `PDB2NET_UNIPROT_FASTA`: Swiss-Prot FASTA
- `PDB2NET_BLAST_DB`: folder containing `uniprot_db.pin`, `uniprot_db.phr`,
  and `uniprot_db.psq`
- `PDB2NET_BLASTP`: `blastp` executable path or command name
- `PDB2NET_BLAST_CACHE_PATH`: writable SQLite cache path, recommended for
  server deployments
- Optional DIAMOND/UniRef90 fallback:
  - `PDB2NET_DIAMOND_ENABLED=true`
  - `PDB2NET_DIAMOND`: `diamond` executable path or command name
  - `PDB2NET_DIAMOND_UNIREF90_DB`: local UniRef90 `.dmnd` path or prefix

Use `--headless` or `PDB2NET_OPEN_IN_CYTOSCAPE=false` for web jobs. A running
Cytoscape desktop UI should not be required on servers.

If DIAMOND/UniRef90 is enabled, the worker image or host must provision the
DIAMOND executable and database ahead of time. PDB2Net validates the configured
paths but does not download UniRef90 or build the `.dmnd` database.

For public workers, set hard limits in the server-controlled JSON configuration.
The package defaults are `null` (unlimited) for backward compatibility:

```json
{
  "workers": {"parsing": 2},
  "resource_limits": {
    "max_input_files": 50,
    "max_total_input_bytes": 2500000000,
    "max_single_input_bytes": 600000000,
    "max_processing_batch_bytes": 600000000,
    "max_total_input_expanded_bytes": 2500000000,
    "max_single_input_expanded_bytes": 600000000
  },
  "network_annotations": {
    "use_embedded_sifts": true,
    "tooltip_fields": ["uniprot"],
    "max_tooltip_segments_per_database": 20
  },
  "combined_graph_limits": {
    "max_nodes": 50000,
    "max_edges": 500000
  }
}
```

Compressed upload limits and expanded limits are independent. Gzip input is
streamed to validate decompression, CRC, nested-gzip rejection, and expansion
limits before scientific processing. The processing-batch limit bounds the sum
of expanded file sizes parsed, annotated,
and distance-checked at once. After each batch, atom coordinates are discarded;
only compact chain metadata and interaction summaries remain for final network
exports. It is a scheduling bound, not a guaranteed process-RSS ceiling, because
parsed in-memory expansion depends on structure complexity.

Run one scientific job per 8-CPU/32-GB worker until target-host measurements
show that parallel jobs are safe. Reference lookup tables, imported scientific
libraries, accumulated interaction summaries, and combined-network preparation
create a non-zero memory floor outside the raw-file batch limit. A representative
RSS/load test is therefore a deployment gate for the 2.5-GB public profile.

## Optional PDB-ID Store Mode

The optional `precompute`/`assemble` path is documented in
[`precomputed_store.md`](precomputed_store.md). Schema 2 caches compact geometry
and standard-profile edges separately from refreshable chain annotations, not
coordinates or detailed atom-pair CSVs. It therefore rejects detailed
interaction export and enforces fixed per-entry and per-request artifact
ceilings.

For a shared deployment, run offline precompute and the worker under the same
UID/GID or a common store-owning group. Use setgid directories (normally mode
`2770`) so nested profile/shard directories inherit the group. Published Core
artifacts use `0640` and lock files use `0660`; directory write/search access is
required for lazy population and atomic replacement.

Set `PDB2NET_REFERENCE_MANIFEST_ID` to an immutable manifest identity covering
reference checksums, the Core/worker build, and exact BLAST+/DIAMOND versions.
Do not reuse it after a reference or toolchain change. Mirror refreshes should
flow through staging, incremental precompute, scientific/contract validation,
and an external atomic store promotion with the prior generation retained for
rollback. Core never automatically removes cached IDs that disappear from a
new mirror; obsolete-ID policy belongs to the operator's staging catalog and
retention workflow.

For a managed server mirror, keep exactly one current Core mmCIF gzip per entry
at `entries/<shard>/<extended-id>/structures/<extended-id>.cif.gz`; for example
`entries/ab/pdb_00001abc/structures/pdb_00001abc.cif.gz`. The targeted resolver
also accepts documented legacy layouts for local/backward compatibility and
fails closed when more than one candidate matches. This server fast path is
additive: ordinary local `pdb2net run --input-dir ...`, including enriched
NextGen uploads, remains independent of a precomputed store and all web
infrastructure.

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

Contract `1.2` keeps the existing fields and adds `identities` and
`structure_inputs` alongside annotation counts, reference identity/metadata,
resource observations, and `skipped_outputs`. Combined graphs
over the configured node or edge cap are listed in `skipped_outputs` and
`warnings`; the job still succeeds and retains its other outputs. `summary.json`
also includes `output_contract_version` and `pdb2net_version` so workers can
detect the structure they are reading. The `config` object mirrors the relevant
runtime settings from the internal `run_summary.json`, which lets a web UI show
the exact defaults or advanced settings used for a job.

## Web UI Settings Model

The recommended webserver model is:

- keep the default submission form simple and do not expose scientific tuning
  unless the user opens advanced settings
- validate advanced settings in the webserver worker
- write validated values into a per-job JSON config
- pass that config with `--config <job>/pdb2net_config.json`

For upload jobs, expose the embedded NextGen tooltip fields as a validated
allowlist: UniProt on by default, with Pfam, CATH, and SCOP2 opt-in. Hide that
selector for server PDB-ID jobs because the Core-only archive is intentionally
not an enriched NextGen mirror. Other useful advanced settings are
distance/contact thresholds, selected network types, detailed interaction
export, and the optional DIAMOND/UniRef90 fallback.
Do not pass arbitrary user-provided shell strings to PDB2Net. Keep reference
data paths, DIAMOND database paths, and executable paths server-controlled.

## Notes For Workers

- Run each job in an isolated job directory.
- Do not write uploaded files into the PDB2Net repository.
- Do not require Cytoscape for automated jobs.
- Treat non-zero process exit as a failed job, then inspect
  `outputs/summary.json` or the internal `run_summary.json` when present.
- Keep large reference datasets, BLAST databases, and DIAMOND databases outside
  git.
