# PDB2Net clean-release implementation plan

## Status and authority

This document is the binding hand-off for the next implementation chat. It
captures the decisions already made for cleaning up and testing the existing
PDB2Net Core and Webserver repositories. It is not a request to redesign the
system.

The plan is **closed world**: behavior not required here is out of scope. If an
implementation detail conflicts with an old hardening-branch document or test,
this plan defines the intended target unless changing it would alter the
scientific results of the established raw-file pipeline. In that case the
implementer must stop and ask.

Only these existing repositories are in scope:

- Core: `/home/gregor/projects/PDB2Net`
- Webserver: `/home/gregor/projects/pdb2net-webserver`

No clone, copied repository, “clean repo”, or Git worktree may be created.
Working branches inside the existing repositories are allowed as defined
below. The Core remains standalone; it must not import, call, or depend on PHP,
Apache, MariaDB, Docker, or the Webserver repository.

## 1. Goal and definition of done

The goal is a maintainable, publication-ready PDB2Net Core and a small,
deployable Webserver integration that another person can understand, test,
operate, and modify without learning a custom release framework. The first
server release must support both structure-file uploads and PDB-ID submissions.
PDB-ID jobs use a store produced completely offline; runtime lazy population is
not part of the product.

The implementation is done only when all of the following are true:

1. Core supports local/headless raw-file processing, offline precomputation,
   read-only assembly, a small capability document, and the stable web-output
   layout described below.
2. Core's established scientific raw-file behavior is preserved by regression
   tests; the same fixed-profile input produces scientifically equivalent raw
   and precomputed results.
3. Core is installable outside the source tree, contains the required packaged
   configs but not `config.local.json`, has current publication metadata and
   documentation, and passes its supported Python test/quality matrix.
4. Web exposes exactly two submission modes: `upload` and `pdb_ids`.
   Directory/server-folder input is gone.
5. Web uses possession of a random UUIDv4 job ID as the unlisted result link.
   There is no separate job access token, fragment-token exchange, private-job
   session table, job listing, or searchable job index.
6. Email is required at submission and continues to use the existing external
   `mailtoken` helper integration. The helper and its secrets are not copied
   into either repository.
7. The production worker mounts reference data and the selected precomputed
   store read-only. It has no PDB archive mount and no write permission to the
   store. Exactly one administrative offline process may build a store.
8. The Webserver is rebuilt cleanly from its local `master` history inside the
   same repository, carrying forward only required product, security,
   deployment, and test behavior. The release-readiness ledger, receipts,
   attestation machinery, migration rehearsals, and custom rollback framework
   are absent.
9. The complete Core suite, Web suite, package-install checks, Compose/admin
   checks, and real cross-repository upload and PDB-ID smoke tests pass on tiny
   fixtures.
10. Neither repository contains generated jobs, reference datasets, database
    dumps, BLAST/DIAMOND databases, secrets, machine-local paths, or deployment
    state. No push or server deployment is performed by this plan.

## 2. Starting state

The implementation chat must re-check these facts before changing anything.
They describe the state when this plan was written.

### Core

- Path: `/home/gregor/projects/PDB2Net`
- Branch: `codex/release-hardening`
- HEAD: `38ababe` (`Harden detailed interaction export contract`)
- Working tree: clean before this plan file was added
- Branch position: 13 commits ahead of local `main`, 9 commits ahead of
  `origin/codex/release-hardening`
- Package version: `0.2.0rc4`
- Test inventory: 257 tests collected with `.venv/bin/python -m pytest
  --collect-only -q`
- Large hardening additions include a 1,734-line `precomputed_store.py`, a
  723-line `server_interface.py`, capability schema 2, CLI contract 1, web
  output contract 1.2, and precomputed schema 2.
- Lazy population is currently exposed through `assemble --source-dir
  --populate-missing`; both flags and all supporting runtime writer/lock/archive
  logic are target removals.

### Webserver

- Path: `/home/gregor/projects/pdb2net-webserver`
- Branch: `codex/contract-1.2-hardening`
- HEAD: `38493bd` (`Align webserver with current PDB2Net core`)
- Working tree: `docs/deployment.md` has an uncommitted user change of 64
  inserted and 3 deleted lines. It must not be discarded, stashed, committed,
  or moved without explicit user approval.
- Branch position: 43 commits and 116 changed files relative to local `master`,
  with 37,209 insertions and 934 deletions
- Test inventory: 730 tests collected with the Core virtualenv
- The hardening branch contains required ideas and tests but also the
  overengineered surfaces being removed: the 6,038-line `pdb2netctl.py`, image
  receipts/activation state, append-only release-readiness governance, exact
  Core semantic snapshots, legacy migration audit/rehearsal machinery,
  directory mode, lazy population, and token/session-based job access.
- Local `master` is the clean rebuild base and already contains the basic
  public pages, upload path, worker, MariaDB schema, Docker Compose, direct
  downloads, and `mailtoken` integration.

The collected test counts are inventory baselines, not permission to preserve
tests for behavior that this plan deliberately removes. Retain or rewrite tests
that protect required behavior; delete tests whose only purpose is enforcing a
removed mechanism.

## 3. Binding product and architecture decisions

### Repository boundary

- Core is the scientific product and future publication artifact.
- Web is a separate adapter around the installed Core package and CLI.
- Production Web images install Core through `PDB2NET_PACKAGE_SPEC` (a package
  release or immutable VCS revision). They never depend on a sibling source
  checkout. The existing local Compose override may use the sibling checkout
  only for development.
- Most future scientific or CLI changes belong in Core. Web changes only when
  the small compatibility contract, job parameters, or output consumption must
  change.
- Core owns parsing and semantic validation of PDB/mmCIF, CX2, CSV, and the
  precomputed artifact. Web owns HTTP input validation, queueing, filesystem
  confinement, artifact delivery, and operational configuration.

### Target contract and package versions

The cleanup intentionally breaks unreleased RC4 server-facing interfaces and
re-versions them cleanly:

| Contract | Target |
| --- | --- |
| Core package during cleanup | `0.2.0rc5` |
| Capabilities schema | `3` |
| CLI contract | `2` |
| Web output contract | `2.0` |
| Precomputed artifact schema | `3` |

There is no compatibility adapter for capability schema 2, CLI contract 1,
output contract 1.2, precomputed schema 2, or the removed RC4
`server_interface` tree. No released public consumer depends on those RC
interfaces.

### Scientific behavior

- Preserve the existing config layering and established raw-file scientific
  behavior unless this plan explicitly changes a server interface.
- Keep Gemmi parsing, SIFTS/FASTA annotation, Swiss-Prot BLAST fallback,
  cKDTree distance calculations, network construction, and headless CX2 export.
- Keep DIAMOND/UniRef90 as an optional annotation fallback. It is disabled by
  default, requires an explicitly configured local `.dmnd` database, and
  remains labelled lower-confidence. Nothing downloads or builds UniRef90.
- Do not add benchmarking, paper experiments, new scientific methods, or
  generalized experiment infrastructure. Those are later projects.
- Keep live Cytoscape support, but move `py4cytoscape` from mandatory runtime
  dependencies into a `cytoscape` optional extra. Headless use must work with
  the base installation.

## 4. Public Core contracts

### 4.1 CLI contract 2

The installed console command is `pdb2net`; `python3 -m pdb2net` exposes the
same parser. The legacy `python3 -m pdb2net.main` entry point remains functional
for existing local users but is not used by Web automation.

The top-level public commands remain exactly:

| Command | Required interface and behavior |
| --- | --- |
| `run` | Requires `--input-dir` and `--output-dir`; accepts optional `--web-output-dir`, `--config`, `--headless`, current reference/search overrides, current distance/network/model/export options, and optional DIAMOND options. Processes local `.pdb`, `.ent`, `.cif`, and `.mmcif` inputs, including supported gzip forms. |
| `precompute` | Requires `--input-dir` and `--store`; accepts `--config`, `--recursive`, and `--headless`. It is the only command that writes precomputed entries and is an offline administrative operation. |
| `assemble` | Requires `--store`, one or more repeated `--pdb-id`, and `--output-dir`; accepts `--web-output-dir`, `--config`, and `--headless`. It only reads and validates an already-published store. |
| `capabilities` | Accepts `--json`; produces the configuration-free document below without loading local configs or scientific reference data. |

Global `--version` remains. `assemble --source-dir` and
`assemble --populate-missing` are removed from the parser, implementation,
documentation, and tests. There is no hidden environment equivalent and no
runtime archive resolver.

Exit status is zero only on complete success and non-zero on validation,
processing, or output failure. Structured error codes remain in run/output
diagnostics; no server-specific logging or database behavior enters Core.

### 4.2 Capability schema 3

`pdb2net capabilities --json` has exactly these required top-level fields and
no `server_interface` subtree:

```json
{
  "capabilities_schema_version": "3",
  "pdb2net_version": "0.2.0rc5",
  "cli_contract_version": "2",
  "output_contract_version": "2.0",
  "precomputed_artifact_schema_version": "3",
  "commands": ["run", "precompute", "assemble", "capabilities"],
  "input_formats": [".pdb", ".cif", ".mmcif", ".pdb.gz", ".cif.gz", ".mmcif.gz"],
  "network_outputs": [
    "chain_per_pdb",
    "protein_per_pdb",
    "combined_chain_network",
    "combined_protein_network"
  ],
  "structure_model_policies": ["first", "all"],
  "features": [
    "headless_cx2",
    "detailed_interactions",
    "embedded_sifts",
    "diamond_uniref90",
    "precomputed_store"
  ]
}
```

The runtime value of `pdb2net_version` comes from the installed package. Field
order is not semantically relevant. A future schema may add top-level fields;
Web schema-3 validation checks required scalar versions and that required list
members are present. It does not require byte-for-byte or whole-document
equality and does not carry a duplicate 28-KiB semantic snapshot.

Web keeps one small checked-in compatibility requirement document as the only
source consumed by worker startup, Doctor, and tests:

```json
{
  "capabilities_schema_version": "3",
  "cli_contract_version": "2",
  "output_contract_version": "2.0",
  "precomputed_artifact_schema_version": "3",
  "required_commands": ["run", "precompute", "assemble", "capabilities"],
  "required_input_formats": [".pdb", ".cif", ".mmcif", ".pdb.gz", ".cif.gz", ".mmcif.gz"],
  "required_network_outputs": ["chain_per_pdb", "protein_per_pdb", "combined_chain_network", "combined_protein_network"],
  "required_structure_model_policies": ["first", "all"],
  "required_features": ["headless_cx2", "precomputed_store"]
}
```

It records no exact Core commit and no duplicate CLI/config semantics. The
installed package revision is pinned by the image/package build mechanism and
reported for diagnostics, not encoded into compatibility logic.

### 4.3 Web output contract 2.0

Both `run --web-output-dir` and `assemble --web-output-dir` publish the same
stable layout:

```text
outputs/
├── summary.json
├── networks/
│   └── *.cx2
└── interactions/
    └── *.csv
```

`summary.json` is valid UTF-8 JSON, contains no non-finite numbers or absolute
server paths, and has these top-level fields:

```text
output_contract_version  string, exactly "2.0"
pdb2net_version          non-empty string
status                   "success" or "failed"
started_at               ISO-8601 string or null
finished_at              ISO-8601 string or null
input_files              array of portable labels
identities               array of canonical structure identities
structure_inputs         array of portable input metadata
networks                 array of relative networks/*.cx2 paths
interactions             array of relative interactions/*.csv paths
artifacts                object with networks[] and interactions[]
runtime_analysis         relative path or null
counts                   object with networks, interactions, structures,
                         chains, and skipped_outputs integers
annotations              JSON object
references               object containing non-empty manifest_id on success;
                         precomputed jobs also identify schema/profile
resources                JSON object
skipped_outputs          array
config                   public scientific runtime snapshot only
errors                   array of {code, message[, output]}
warnings                 array of {code, message[, output]}
```

Each network artifact record contains `path`, `size_bytes`, `sha256`, `nodes`,
and `edges`. Each interaction artifact record contains `path`, `size_bytes`,
`sha256`, `rows`, and `columns`. On `success`, `errors` is empty and every
listed regular file exists beneath the output root. On `failed`, no result
artifact is published and diagnostics are bounded and path-redacted.

Core is the single semantic validator and producer for CX2 and CSV. Web checks
the summary version/status/types, relative-path confinement, regular-file and
symlink rules, configured count/size ceilings, allowed extensions, and that
file sizes/hashes match the Core records. Web does not implement a second CX2
aspect parser or CSV row/column parser.

## 5. Precomputed store schema 3

### Purpose and ownership

The store accelerates PDB-ID jobs for one fixed scientific profile. It is not a
general cache and not a replacement for the authoritative source archive.

- One store directory contains exactly one full scientific profile.
- One administrative `precompute` process is the only writer.
- Runtime `assemble` and every Web worker are read-only consumers.
- The PDB archive is mounted only into the administrative precompute context,
  never into the runtime worker.
- Missing entries fail with `PRECOMPUTED_ENTRY_MISSING`; they are never fetched,
  parsed, or populated at request time.
- Store promotion is an operator action: build a new versioned directory,
  validate it, update the configured host path/symlink outside Core, and
  recreate the worker with that directory mounted read-only. Core implements no
  release pointer, receipt, rollback state machine, or cleanup policy.

### Fixed PDB-ID scientific profile

The full profile fingerprint covers all of these values:

- artifact schema 3;
- current scientific, parser, interaction, and annotation semantic versions;
- concrete Gemmi version;
- source scope `asymmetric_unit`;
- `structure_model_policy: first`;
- C-alpha radius `12.0` Å;
- all-atom radius `5.0` Å;
- protein/protein minimum C-alpha neighbours `10`;
- protein/protein minimum all-atom contacts `1`;
- protein/nucleic-acid minimum all-atom contacts `1`;
- nucleic-acid minimum all-atom contacts `1`;
- embedded SIFTS enabled;
- Swiss-Prot BLAST search policy and its thresholds;
- optional DIAMOND policy, disabled by default and included when enabled;
- a required immutable `reference_manifest_id` identifying the exact PDB FASTA,
  SIFTS, Swiss-Prot, optional UniRef90, and external search-tool versions.

The manifest's `profile` object has this exact field structure; the named
semantic-version and search-policy values come directly from their existing
Core owners rather than being duplicated as new constants:

```text
artifact_schema_version
scientific_pipeline_version
interaction_pipeline_version
annotation_pipeline_version
source_scope
structure_model_policy
parser: {semantics, gemmi_version}
distance_thresholds: {ca_radius, all_atoms_radius}
interaction_filters: {
  protein_protein_min_ca_neighbors,
  protein_protein_min_all_atom_contacts,
  protein_nucleic_acid_min_all_atom_contacts,
  nucleic_acid_min_all_atom_contacts
}
annotations: {
  reference_manifest_id,
  use_embedded_sifts,
  search_policy
}
```

Network export selections and tooltip display choices are not part of the
profile because networks are assembled from the same stored chains and edges.
Detailed atom/residue interaction export is not supported by `assemble`.

### Filesystem layout

There is no `profiles/<hash>/` nesting:

```text
store/
├── manifest.json
└── entries/
    └── ab/
        └── pdb_00001abc.json.gz
```

IDs are deduplicated canonical lowercase Extended PDB IDs. The two-character
shard follows the existing identity/sharding helper. Entries are deterministic
gzip-compressed JSON and are written atomically. Symlinks and paths escaping
the store root are rejected.

`manifest.json` is written atomically **last**, only after every discovered
source in the requested offline build has either produced or reused a valid
entry and the build has zero failures. A failed build leaves entries available
for a retry but no published manifest; `assemble` rejects such a directory. A
retry may reuse entries only when their source hash and complete profile ID
match. Building against an already-published store fails; updates use a new
operator-selected directory.

The manifest has exactly:

```json
{
  "artifact_schema_version": "3",
  "created_at": "<UTC ISO-8601>",
  "producer": {"name": "pdb2net", "version": "0.2.0rc5"},
  "profile_id": "<sha256 of canonical profile JSON>",
  "profile": {"...": "complete fixed scientific and annotation profile"},
  "entry_count": 1
}
```

Each entry has exactly these top-level fields:

```text
artifact_schema_version  "3"
created_at               UTC ISO-8601 string
producer                 {name: "pdb2net", version: <package version>}
profile_id               manifest profile SHA-256
pdb_id                   canonical lowercase Extended PDB ID
structure_identity       canonical existing StructureIdentity mapping
source                   {sha256, size_bytes, scope}
geometry                 {structure, interactions}
annotations              {references, chains}
counts                   {chains, interactions}
```

The nested compact chain, interaction, structure-identity, annotation segment,
and count semantics remain scientifically equivalent to the validated schema-2
representation. Geometry and annotations remain separate sections, but there
is one combined store profile and no refreshable annotation sub-profile. The
store contains no atom coordinates, cKDTree state, detailed atom/residue pairs,
absolute source paths, source-machine filenames, or rendered CX2/CSV outputs.

Retain the fail-closed validation ceilings from the current implementation: 64
MiB compressed and 128 MiB expanded per entry, 50,000 chains and 500,000 edges
per entry, and 256 MiB expanded JSON, 200,000 chains, and 2,000,000 edges across
one assemble request. The configured total-input limit may be lower. These are
security ceilings, not capacity promises; changing them requires a separate
user decision.

## 6. Webserver contract

### Submission modes

The only accepted `input_mode` values are:

1. `upload`: one to 50 `.pdb`, `.cif`, or `.mmcif` files, including supported
   gzip forms, processed with `pdb2net run`.
2. `pdb_ids`: one to 50 unique PDB IDs. Classic four-character aliases are
   normalized at submission to canonical lowercase Extended IDs and processed
   with `pdb2net assemble`.

The worker constructs argv arrays directly and never invokes a shell. The two
real-mode command shapes are exactly:

```text
pdb2net run
  --input-dir /jobs/<jobid>/inputs
  --output-dir /jobs/<jobid>/work
  --web-output-dir /jobs/<jobid>/outputs
  --config /jobs/<jobid>/pdb2net_config.json
  --headless

pdb2net assemble
  --store /precomputed
  --pdb-id <canonical-id> [--pdb-id <canonical-id> ...]
  --output-dir /jobs/<jobid>/work
  --web-output-dir /jobs/<jobid>/outputs
  --config /jobs/<jobid>/pdb2net_config.json
  --headless
```

Directory input and all `server_dir`, `/server-input`, allowed-directory-base,
directory-mode feature flags, frontend controls, worker branches, Compose
mounts, documentation, and tests are removed.

Both supported modes are product behavior, not deployment feature flags.
Remove `PDB2NET_ENABLE_PDB_ID_MODE`, `PDB2NET_ENABLE_DIRECTORY_MODE`,
`PDB2NET_ALLOWED_SERVER_INPUT_BASE`, `PDB2NET_PDB_ARCHIVE_HOST_DIR`, and
`PDB2NET_PRECOMPUTED_POPULATE_MISSING`. Keep one required
`PDB2NET_PRECOMPUTED_HOST_DIR` for the published read-only store.

Upload jobs keep the existing public scientific choices and bounds from
`contracts/web-job-contract.json`: cutoffs, four network selections,
interaction filters, `first|all` model policy, detailed-interaction export, and
annotation tooltip settings. Its schema version remains `1`, and that file
remains the single source of UI and server parameter defaults/bounds.

PDB-ID jobs always use the fixed store profile. They require first-model
semantics and standard cutoffs/filters, prohibit detailed interactions, and may
choose only which of the four networks are rendered. The frontend must explain
this and must not submit incompatible hidden values.

Keep these currently hardened upper bounds. Changing one requires a separate
user decision:

- 50 uploaded files;
- 600,000,000 bytes per compressed/uploaded file;
- 2,500,000,000 uploaded bytes per job;
- the same current per-file and aggregate expanded-byte ceilings;
- 50 PDB IDs;
- 20 active jobs globally, three active jobs per email, and ten submissions per
  email per hour by default, with the existing bounded environment overrides;
- 5,000,000 detailed-interaction rows and 500,000,000 detailed-interaction
  bytes per job;
- a 6,000,000,000-byte runtime free-space reserve;
- a 1,000,000-byte summary, at most 1,000 artifacts, 512 MiB per artifact, and
  2 GiB across published artifacts;
- the current filename, gzip/archive, symlink, and path-confinement guards.

### Job identity and access

- Submission creates a lowercase RFC-4122 UUIDv4 using a cryptographically
  secure system source.
- Possession of that UUID is the only result-access capability. Result links use
  the existing `my_jobs.php?jobid=<uuid>` shape; status and download endpoints
  accept the same validated UUID.
- Remove `access_token_hash`, fragment `#token=...`, token exchange,
  `pdb2net_job_sessions`, browser authorization tables/cookies, token TTL
  settings, and token recovery UI.
- Do not expose an endpoint or page that lists jobs. Do not link job pages from
  an index or sitemap; result pages use `noindex`/`no-store` responses.
- This is unlisted-link access, not authenticated private storage. User-facing
  text and documentation must say so accurately.
- Every status/download lookup is an exact UUID lookup. Download paths are
  derived from the configured jobs root and allowed artifact records, never
  trusted from user-controlled `result_path` or filenames.

### Email

- A syntactically valid normalized email address is required for both modes.
- Production configuration requires the external `mailtoken` helper path and
  public base URL. The helper mount is read-only and available only to the Web
  container; OAuth credentials and helper implementation stay outside the repo.
- Submission emails the UUID result link. Mail failure is logged without
  secrets and the same unlisted link is shown once in the submission response
  as a fallback; failure does not invent a second credential scheme.
- Keep queue/admission limits by normalized email. Do not expose email in job
  URLs, summaries, logs, or downloaded artifacts.

### Fresh MariaDB schema

There is no migration requirement for existing jobs or databases. Replace the
historical migration chain with one current fresh-install schema and a minimal
ordered migration mechanism for future releases. Do not keep legacy fixtures,
rehearsal databases, migration provenance hashes, database identity/backup
attestation tables, or compatibility transforms.

The current `pdb2net_jobs` table contains only fields needed by the two modes,
queue lifecycle, email notification, cleanup, and diagnostics:

```text
id                 BIGINT UNSIGNED AUTO_INCREMENT primary key
jobid              CHAR(36) ASCII/BINARY, unique, non-null UUIDv4
email              VARCHAR(320), non-null
input_mode         ENUM('upload','pdb_ids'), non-null
input_files        JSON nullable; upload labels only
pdb_ids            JSON nullable; canonical IDs only
input_file_count   INT non-null default 0
input_bytes        BIGINT UNSIGNED non-null default 0
params             JSON non-null
status             ENUM('pending','in_progress','completed','failed')
status_message     TEXT nullable
notification_status ENUM('pending','sent','failed') default 'pending'
heartbeat_at       TIMESTAMP nullable
failure_code       VARCHAR(64) nullable
core_version       VARCHAR(64) nullable
contract_version   VARCHAR(32) nullable
created_at         TIMESTAMP non-null default current timestamp
started_at         TIMESTAMP nullable
finished_at        TIMESTAMP nullable
expires_at         TIMESTAMP nullable; job-retention cleanup only
```

Remove `job_name`, `description`, `server_dir`, `access_token_hash`, and
`result_path`. There is no `pdb2net_job_sessions` table. Filesystem paths are
deterministic from `jobs/<jobid>`. Replace the removed token-TTL setting with
one directly named cleanup setting, `PDB2NET_JOB_RETENTION_DAYS`, defaulting to
30 days; submission sets `expires_at` from it and cleanup uses that timestamp.

The initial schema creates the current table directly. `db migrate` remains a
small forward-only mechanism for migrations introduced after this clean
baseline; it does not replay the discarded hardening migration history.

### Runtime and Compose

Use standard Docker Compose with three runtime services:

- Apache/PHP Web, serving only `app/public`;
- Python worker with the installed Core package;
- MariaDB on the internal Compose network.

Runtime mounts are minimal:

- job data: writable where submission/worker responsibilities require it;
- scientific references and BLAST/optional DIAMOND data: worker read-only;
- precomputed store: worker read-only;
- `mailtoken` helper: Web read-only;
- database storage: MariaDB only.

There is no worker PDB archive mount, no writable precomputed mount, no
`server-input` mount, and no PHP access to scientific reference/store data.
Production PDB-ID support is not a lazy-population feature toggle: Doctor and
startup fail clearly when the configured published store is absent or
incompatible.

## 7. Lean operations interface

Replace the 6,038-line controller and related release-state modules with a thin
Python standard-library CLI plus small responsibility-focused helpers. It may
wrap Docker Compose, MariaDB tools, filesystem setup, Core precompute, and HTTP
smokes; it must not become a release orchestrator or persistence framework.

The supported commands are exactly:

```text
./pdb2netctl init
./pdb2netctl doctor
./pdb2netctl build
./pdb2netctl up
./pdb2netctl down
./pdb2netctl status
./pdb2netctl logs
./pdb2netctl db init
./pdb2netctl db migrate
./pdb2netctl db backup
./pdb2netctl db restore
./pdb2netctl precompute
./pdb2netctl cleanup [--apply]
./pdb2netctl smoke
./pdb2netctl render-apache
```

Behavior:

- `init` creates/validates the required runtime directories and an environment
  file from the documented template without touching existing data.
- `doctor` performs direct configuration, permissions, Compose, Core
  capability, store-manifest, database, mail-helper, and optional online health
  checks. It reports pass/fail; it emits no receipt.
- `build`, `up`, `down`, `status`, and `logs` are transparent Compose wrappers.
- `db init|migrate|backup|restore` use the checked-in current schema and normal
  MariaDB dump/restore tools. Restore requires an explicit target and normal
  confirmation; there is no custom provenance ledger.
- `precompute` invokes the installed Core offline against explicit archive,
  reference/config, and target-store paths. It is not available inside the
  runtime worker.
- `cleanup` is a dry run unless `--apply` is supplied and only targets expired
  deterministic job directories and matching database rows.
- `smoke` exercises health, one tiny upload, and one tiny PDB-ID job through the
  public route using explicit test inputs.
- `render-apache` renders the small native Apache route/template; it does not
  install or attest the host configuration.

Remove image receipts, image activation state, deployment records, immutable
evidence capture, JSON release trackers, release finding IDs, automatic
rollback state, migration rehearsals, backup identity catalogs, and validators
that exist solely for those mechanisms. Replace them with CI, ordinary tests,
and one human-readable deployment checklist.

## 8. Implementation phases

### Phase 0 — Safety, branches, and baseline

1. Re-read both `AGENTS.md` files and this plan; inspect both full statuses and
   diffs.
2. In Core, create `codex/clean-release` from the current
   `codex/release-hardening` HEAD inside the existing checkout. Add this plan as
   the first docs-only commit on that branch.
3. In Web, create branch reference `codex/clean-rebuild` from local `master`.
   Before switching the existing checkout, stop and ask the user how the dirty
   hardening-branch `docs/deployment.md` change should be preserved. Do not
   choose commit/stash/patch/discard on the user's behalf. After the user has
   made that custody decision, switch the same repository checkout to the new
   branch.
4. Record baseline test commands/results. Do not modify existing hardening
   branches, remote branches, tags, or user data.

### Phase 1 — Core public-contract reduction

1. Apply the target version constants and capability schema 3.
2. Remove `server_interface.py`, its exact semantic tree, sync/snapshot
   consumers, and tests that only enforce that tree. Preserve reusable
   scientific constants in the modules that actually own them.
3. Keep the four commands, remove lazy CLI options, and update CLI/help/docs and
   focused contract tests.
4. Make web output 2.0 the single stable Web adapter and keep internal
   timestamped output backward-compatible for local users.
5. Run focused CLI, config, output, and regression tests before proceeding.

### Phase 2 — Core precompute simplification

1. Replace the monolithic module with the minimum package split already agreed:
   `precomputed/schema.py` for profile/schema validation,
   `precomputed/io.py` for bounded atomic reads/writes,
   `precomputed/build.py` for the sole offline writer, and
   `precomputed/assemble.py` for the read-only consumer. `__init__.py` exposes
   only the small public functions needed by CLI/tests.
2. Implement schema 3, the one-profile layout, manifest-last publication, and
   retry/reuse rules.
3. Delete archive lookup, lazy population locks, runtime mutation, profile
   namespaces, annotation refresh-in-place, and their tests/docs.
4. Reuse the normal parser/annotation/distance pipeline rather than copying its
   scientific logic. Keep geometry and annotations separately represented.
5. Prove raw/precomputed equivalence and fail-closed behavior for missing,
   partial, corrupt, oversized, symlinked, stale-profile, and duplicate-ID
   artifacts.

### Phase 3 — Core maintainability and publication readiness

1. Streamline pipeline/import boundaries only where required to remove the
   server-interface/precompute coupling; do not rewrite working scientific
   algorithms for style.
2. Keep BLAST and optional DIAMOND behavior isolated behind their existing
   annotation responsibilities; DIAMOND remains off by default and
   lower-confidence.
3. Make Cytoscape an optional extra; verify base headless installation.
4. Add Ruff to the development quality checks without mass stylistic churn.
5. Complete `pyproject.toml`, README, CHANGELOG, license/CITATION, architecture,
   development, precomputed-store, and biologist-friendly usage documentation
   needed for a future public package/code release.
6. Test Python 3.11 and 3.12 in CI, build/install outside the source tree, and
   verify packaged configs and console commands.

### Phase 4 — Clean Web rebuild from local master

1. Start from the small `master` implementation, not the 43-commit hardening
   branch. Use the latter only as read-only reference for tested security and
   PDB-ID/mail behavior; do not merge or cherry-pick the governance stack.
2. Retain the public site, upload flow, installed-Core boundary, direct
   downloads, MariaDB queue, Compose base, and external `mailtoken` integration.
3. Port only required upload hardening, canonical PDB-ID parsing, admission
   limits, worker lifecycle, path confinement, safe output publication, health,
   and the two-mode UI.
4. Remove directory mode and implement PDB-ID assembly against the read-only
   schema-3 store.
5. Replace exact Core snapshots with the small compatibility requirement
   document. Simplify output collection to metadata/security validation only.
6. Keep modules responsibility-focused and direct. Do not recreate the removed
   release framework under new names.

### Phase 5 — Web database, UUID access, mail, and operations

1. Install the fresh schema and delete token/session/legacy migration code and
   tests.
2. Implement UUID-only unlisted links, exact UUID lookups, noindex/no-store
   responses, confined downloads, and no job listing.
3. Require email and wire the UUID link through existing `mailtoken`, with the
   one-time response fallback on delivery failure.
4. Reduce Compose/runtime mounts to the stated boundary and make the
   precomputed store read-only.
5. Implement the lean `pdb2netctl` command set and replace release governance
   with concise architecture, deployment, mailtoken, precompute, and operator
   documentation.
6. Port tests only for target behavior. Delete tests tied exclusively to
   removed receipts, trackers, sessions, directory input, lazy population, and
   legacy migrations.

### Phase 6 — Cross-repository integration and final review

1. Build/install Core as a package rather than importing the source checkout.
2. Validate the small capability contract from Web Doctor and worker startup.
3. Build a tiny schema-3 store offline from repository fixtures.
4. Run one real headless upload job and one real PDB-ID job through the complete
   queue/worker/output/download path.
5. Verify an absent PDB ID fails clearly without writing to the store or
   attempting archive access.
6. Review both diffs against this plan, remove unjustified files/options/layers,
   run all final checks, and produce a concise deployment-readiness handoff. Do
   not push or deploy.

## 9. Test and acceptance matrix

### Core required checks

- `.venv/bin/python -m pytest` passes.
- Ruff check passes on the agreed maintained source/test scope.
- Config layering order and environment overrides remain tested.
- `pdb2net --version` and help for all four commands work from an installed
  package.
- Capability schema 3 exact required shape and configuration-free import are
  tested.
- `run` on tiny PDB/mmCIF fixtures works headlessly and writes valid internal
  plus web-output 2.0 artifacts.
- Offline precompute/assemble equivalence covers chain data, annotations,
  interaction edges, network nodes/edges, identities, and warnings relevant to
  scientific output.
- Store tests cover fixed profile hashing, manifest-last publication, retry of
  an unpublished build, read-only assembly, canonical/duplicate IDs, missing
  entry, corrupt JSON/gzip, wrong schema/profile/identity, size/count ceilings,
  symlinks/path escape, and no coordinate/detail data.
- BLAST behavior remains covered with tiny fixtures/mocks; optional DIAMOND has
  focused disabled/default/lower-confidence tests and its existing tiny local
  end-to-end fixture where available.
- Python 3.11/3.12 CI builds the package and verifies from outside the source
  tree that distributable configs exist and `config.local.json` does not.

### Web required checks

- The Web Python/PHP test suite passes from the existing repository.
- PHP syntax checks pass for every PHP source.
- Submission tests cover upload and PDB-ID success/failure, required email,
  canonical ID deduplication, exact limits, gzip expansion, filenames,
  symlinks, partial uploads, admission limits, and removal of directory mode.
- DB tests start from an empty MariaDB and cover schema creation, atomic job
  claim, heartbeat/stale recovery, completion/failure, notification status,
  cleanup, and no token/session tables or legacy migrations.
- UUID access tests cover valid UUIDv4 links, invalid IDs, unknown IDs, no
  listing, noindex/no-store, path traversal, symlinks, allowlisted artifacts,
  and absence of tokens/cookies/session authorization.
- Mail tests use a stub helper and cover required recipient, correct UUID link,
  successful delivery, safe failure logging, and response fallback without
  leaking email or adding a credential.
- Worker tests assert exact argv arrays for `run` and `assemble`, never a shell,
  never lazy flags, and never a PDB archive.
- Output collector tests cover contract/version/status, bounded summary,
  relative paths, regular files, extensions, sizes, hashes, and no duplicate
  semantic CX2/CSV parsing.
- Compose/config tests prove Apache serves only `app/public`, MariaDB is
  internal, reference/store mounts are worker read-only, mailtoken is Web
  read-only, and forbidden archive/server-input/writable-store mounts are
  absent.
- Every lean `pdb2netctl` command has focused argument/error tests; destructive
  cleanup/restore behavior requires its explicit flag/target.

### Cross-repository required checks

- Install the candidate Core package into the real worker image/environment.
- `pdb2net capabilities --json` satisfies the small Web requirements without
  whole-document equality.
- Real upload E2E: submit tiny structure plus email, queue, run Core headlessly,
  finish, view status by UUID, and download every summary-listed artifact.
- Real PDB-ID E2E: build a tiny store offline, mount it read-only, submit classic
  and Extended aliases, assemble, finish, and download outputs.
- Negative PDB-ID E2E: missing/corrupt/profile-mismatched entry fails with a
  useful public error and leaves the store byte-for-byte unchanged.
- Core/Web output hashes, counts, versions, config snapshot, status, warnings,
  and error handling agree.
- `docker compose config`, image builds, service health, and the lean smoke
  command pass without a live Cytoscape UI or external dataset download.

## 10. Explicit non-goals

- No lazy population, runtime PDB download, runtime archive resolution, or
  worker write access to precomputed data.
- No full PDB archive download, full-store build, UniRef90 download/build, large
  batch run, performance benchmark, scientific benchmark, paper experiment, or
  benchmark framework in this task.
- No change to scientific thresholds/algorithms except expressing the already
  agreed fixed PDB-ID profile.
- No new annotation database, network type, frontend framework, API framework,
  ORM, queue broker, object store, Kubernetes layer, cloud service, telemetry
  platform, or authentication system.
- No Directory input mode.
- No private-token/session access model and no job listing.
- No compatibility with unreleased RC4 server contracts or existing Web
  databases/jobs.
- No custom release attestation, receipts, evidence ledger, rollout state
  machine, migration rehearsal framework, or automated rollback framework.
- No server deployment, DNS/TLS/Apache installation, real mail delivery,
  production database restore, public launch, tag, release upload, or Git push.
- No decision about the paper text or benchmarking methodology beyond keeping
  the code and test seams clean enough to support that later work.

## 11. Simplicity guardrails and stop conditions

These rules are acceptance criteria, not suggestions:

1. The plan is closed world: if a feature is not required above, omit it.
2. A possible future use case is not a current requirement.
3. Add no dependency, configuration layer, compatibility layer, schema surface,
   metadata format, service, infrastructure component, or abstraction without a
   direct numbered plan requirement.
4. Build no generic framework for a single concrete use case.
5. Remove superseded mechanisms; do not preserve them behind wrappers, feature
   flags, aliases, deprecation layers, or “legacy” modules.
6. Keep one source of truth for each contract. Tests consume that source; they
   do not duplicate it as a large snapshot.
7. Prefer a direct function/module and standard tools over registries,
   factories, plugin systems, custom state machines, and orchestration layers.
8. Do not use line count alone as a quality target, but every newly added file,
   module, dependency, config key, schema field, or command must have a direct
   plan justification at its phase gate.
9. Preserve scientific behavior and security bounds; “simpler” does not mean
   silently weakening validation, resource limits, path confinement, or tests.
10. Use no subagents or parallel code edits unless the user explicitly requests
    them. One lead agent owns all implementation decisions and cross-repo
    compatibility.
11. Never discard, reset, stash, commit, or overwrite user work without explicit
    authorization. In particular, stop before handling the current dirty Web
    `docs/deployment.md` needed for the branch switch.
12. Stop and ask before changing this plan, introducing a new architecture
    decision, altering scientific output, expanding scope, performing a large
    download/run, pushing, or deploying.

At each phase boundary, the implementation agent must re-read this document,
run the relevant checks, review the diff for unnecessary moving parts, and
report: what was simplified/removed, what was added and why, test results, and
the next bounded phase.
