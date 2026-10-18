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

The current profile uses scientific pipeline v4, parser v6, interaction v2,
annotation v6 and Swiss-Prot selection policy v3; the public schema remains 3.
Scientific pipeline v4 records polymer-role/type-based lengths (F06/F07):
modified nucleotides contribute to `nt_len`, while nonpolymer amino acids no
longer inflate `aa_len`.

Parser/annotation v6 additionally retain fully modified supported polymers and
use consistent polymer types when residue names leave a chain unknown (F08).
Swiss-Prot policy v3 corrects alignment-span coverage when `qcovs` is unavailable
(F09). The existing profile check rejects older stores; build a new store
offline if one is needed. The search-policy signature also separates older
positive and negative cache rows, without deleting them. F10 only aligns
readiness checks with the existing query-quality rule.

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

The parser fingerprint includes merging noncontiguous author-chain parts within
each model. Entries from a different parser fingerprint are incompatible and
are not reused, including during retries of unpublished builds.

Parser fingerprint v5 introduced Q01 together with interaction fingerprint v2:
protein/DNA/RNA polymer heavy atoms only, explicit zero occupancy excluded,
unknown occupancy retained, and one contact per atom-site pair across all
eligible alternatives. Stores built with parser v4 or interaction v1 are rejected
and require an offline rebuild in a new directory. Q01 itself kept schema 3 and
annotation v5 unchanged and did not invalidate sequence-search caches; the
current F08/F09 profile changes described above additionally apply.
Raw altloc/occupancy metadata is not added to compact entries.

Q02 changes only the reported count of cross-file identity links to zero while
preserving their layout attraction. Those links are generated at export and
are absent from stored geometry. Q02 alone did not require rebuilding
Q01-compatible entries. Entries must still match the current complete profile
described above. Existing
CX2 exports need regeneration to receive the corrected identity count.

Offline build reports and CLI stderr include input warnings such as
`CHAIN_WITHOUT_CONTACT_ATOMS`. Keep the build log with the store's scientific
provenance; these source diagnostics are not stored in compact entries and are
not recreated during `assemble`. A source with no usable contact chains fails
its entry and prevents publication. A valid chain without contacts is allowed.

The annotation fingerprint also distinguishes the correction that prevents
AlphaFold-style filenames from overriding resolved structure identity,
nucleic-acid type or explicit SIFTS ambiguity. Stores using the earlier
annotation fingerprint are rejected and must be rebuilt offline in a new
directory; the schema remains 3. Direct AlphaFold assignments are not stored in
the sequence-search cache, so this correction does not invalidate BLAST hits.

Annotation fingerprint v5 additionally records polymer-only fallback sequences
and residue-based molecule typing. Earlier profiles, including v4, must likewise
be rebuilt offline. Residue-role metadata is internal to raw processing and is
not added to the compact store schema. The sequence-search policy is unchanged:
a corrected query gets a different existing sequence cache key; an identical
query against the same database/policy may still reuse its cached result.

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

Before processing, `precompute` inventories each canonical source once for its
content identity, file signature, compressed and expanded sizes, and initial
SHA-256. Pending batches reuse that inventory without separate identity or
fingerprint discovery passes; parsing still revalidates the expected identity
and content. Failures during per-source transport inspection, identity
preflight, or parsing are isolated; other valid sources may still produce
entries for a later retry, but any failure still prevents manifest publication.
Count and aggregate byte limits apply to the complete canonical input set,
independent of cache hits, and the source fingerprint is checked again
immediately before writing an entry to detect input mutation.

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
