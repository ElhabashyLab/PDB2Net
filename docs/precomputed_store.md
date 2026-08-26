# Precomputed PDB store

The schema-3 store accelerates PDB-ID jobs for one fixed scientific profile.
It is optional and does not change the standalone raw-file pipeline.

## Ownership

`pdb2net precompute` is the only writer. It runs offline with access to the
source archive and reference data. `pdb2net assemble` is a read-only consumer:
it never resolves an archive source, downloads a structure, creates a missing
entry, or changes the store.

Build each generation in a new operator-selected directory. After validation,
promote that directory outside PDB2Net and mount it read-only for workers. A
missing ID fails with `PRECOMPUTED_ENTRY_MISSING`.

## Fixed profile

One store contains exactly one complete profile. The profile fingerprint covers:

- artifact schema 3 and the current scientific, parser, interaction, and
  annotation semantic versions;
- the concrete Gemmi version, asymmetric-unit scope, and first-model policy;
- C-alpha radius 12.0 Å and all-atom radius 5.0 Å;
- minimum contact counts 10/1 for protein pairs and 1 for protein/nucleic-acid
  and nucleic-acid pairs;
- embedded SIFTS enabled;
- the complete Swiss-Prot BLAST and optional DIAMOND search policy;
- a required immutable `reference_manifest_id` for the exact references and
  external search-tool versions.

Nonstandard thresholds, filters, model policy, or disabled embedded SIFTS are
rejected. Network selections and tooltip display fields are not profile inputs;
they affect rendering from the same stored chains and edges. Detailed atom or
residue interaction export is unavailable during assembly.

## Layout and publication

```text
store/
├── manifest.json
└── entries/
    └── ab/
        └── pdb_00001abc.json.gz
```

IDs are canonical lowercase Extended PDB IDs. Entries are canonical JSON in
deterministic gzip files and are written atomically. Symlinked roots, manifests,
entry paths, and files are rejected.

`manifest.json` is written atomically last, only after every discovered source
has produced or reused a valid entry and the build has no failures. A failed
build leaves completed entries but no manifest, so assembly rejects it. Retrying
the same unpublished build may reuse an entry only when its complete profile ID
and source hash match. Once a manifest exists, the store is immutable; updates
use a new directory.

The manifest has exactly:

```json
{
  "artifact_schema_version": "3",
  "created_at": "<UTC ISO-8601>",
  "producer": {"name": "pdb2net", "version": "0.2.0rc5"},
  "profile_id": "<sha256 of canonical profile JSON>",
  "profile": {"...": "complete fixed profile"},
  "entry_count": 1
}
```

Each entry has exactly these top-level fields:

```text
artifact_schema_version
created_at
producer
profile_id
pdb_id
structure_identity
source
geometry
annotations
counts
```

`geometry` contains compact chain structure metadata and filtered interaction
edges. `annotations` contains reference identity and chain annotations. Entries
contain no atom coordinates, residue/atom pairs, source-machine filename or
absolute path, cKDTree state, or rendered CX2/CSV output.

## Commands

Build a new store offline:

```bash
pdb2net precompute \
  --input-dir /srv/pdb/archive-generation \
  --store /srv/pdb2net/store-staging-v1 \
  --config /srv/pdb2net/precompute.json \
  --recursive \
  --headless
```

Assemble one or more published entries:

```bash
pdb2net assemble \
  --store /srv/pdb2net/store-current \
  --pdb-id pdb_00001abc \
  --pdb-id pdb_00002xyz \
  --output-dir /srv/jobs/example/work \
  --web-output-dir /srv/jobs/example/outputs \
  --config /srv/jobs/example/pdb2net_config.json \
  --headless
```

## Validation ceilings

Assembly fails closed at:

- 64 MiB compressed and 128 MiB expanded per entry;
- 50,000 chains and 500,000 edges per entry;
- 256 MiB expanded JSON, 200,000 chains, and 2,000,000 edges per request.

The configured total-input limit may be lower. Corrupt JSON/gzip, wrong schema,
profile or identity, inconsistent counts/endpoints, non-finite values, unsafe
paths, and unpublished stores are rejected before network assembly.
