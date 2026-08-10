# Precomputed PDB graph store

The precomputed store is an optional acceleration layer. It does not replace
the normal local `pdb2net run --input-dir ...` pipeline and has no dependency on
PHP, Apache, MariaDB, Docker, or a particular webserver.

## Stored scientific data

For each PDB entry, schema 2 stores two explicit sections:

- reference-independent compact chain geometry metadata and filtered chain-pair
  interaction edges for one geometry profile;
- refreshable chain annotations plus compact embedded UniProt, Pfam, CATH, and
  SCOP2 residue segments from enriched mmCIF input;
- source SHA-256, source size, source scope (`asymmetric_unit`), package and
  artifact versions, and the canonical Extended PDB identity;
- the geometry profile fingerprint and a separate annotation-profile ID tied to
  the reference manifest and search policy.

It does not store atom coordinates, residue-level atom pairs, detailed
interaction CSVs, or any absolute source-machine path. Consequently,
`assemble` rejects `export_detailed_interactions: true`; use the normal raw-file
mode when those tables are required. Keep the authoritative raw mmCIF archive:
the compact cache cannot reconstruct coordinates, re-run changed geometry
semantics, or recover mmCIF categories that were not selected for the cache.

Schema 2 applies fail-closed artifact limits before assembly:

- at most 64 MiB compressed and 128 MiB expanded JSON per entry;
- at most 50,000 chains and 500,000 compact interaction edges per entry;
- at most 256 MiB expanded JSON, 200,000 chains, and 2,000,000 edges across one
  assemble request.

The configured `resource_limits.max_total_input_bytes` may impose a lower limit
on the sum of requested compressed entries. These are validation ceilings, not
capacity targets or guarantees about peak process memory.

Artifacts are gzip-compressed JSON and deterministically sharded by PDB ID:

```text
store/
└── profiles/
    └── <profile-sha256>/
        ├── manifest.json
        └── entries/
            └── ab/
                └── pdb_00001abc.json.gz
```

Different geometry profiles coexist under different SHA-256 namespaces. They
are never silently mixed or overwritten. Annotation/reference changes update
the validated annotation section within the matching geometry namespace and
reuse geometry only when the source and chain/molecule-type geometry remain
compatible.

## Runtime ownership and permissions

Precompute and worker processes must use the same runtime UID/GID or belong to a
shared group that owns the store. Prepare each store generation and its parent
with the setgid bit so newly created profile, shard, and lock directories inherit
that group. For example, with a deployment-specific group named `pdb2net`:

```bash
install -d -o pdb2net-precompute -g pdb2net -m 2770 /srv/pdb2net/store-staging
```

Run the offline precompute process with that group, and give the worker a
compatible primary or supplementary GID. Core atomically publishes manifests
and entries with mode `0640`; persistent per-entry lock files use `0660`.
Therefore group members need read/write/search permission on the store
directories even though published entries themselves are group-readable rather
than group-writable. Do not populate as `root` and later run the worker under an
unrelated UID/GID. Verify effective ownership and permissions from the actual
worker runtime before enabling lazy population.

## Commands

Populate a store from a directory. Parsing, annotation, BLAST/DIAMOND fallback,
and distance detection use the same bounded batching as a normal run:

```bash
pdb2net precompute \
  --input-dir /srv/pdb/archive \
  --store /srv/pdb2net/precomputed \
  --config /srv/pdb2net/config.json \
  --recursive \
  --headless
```

Assemble cached graph entries into the normal timestamped output layout and,
optionally, the stable web output contract:

```bash
pdb2net assemble \
  --store /srv/pdb2net/precomputed \
  --pdb-id pdb_00001abc \
  --pdb-id pdb_00002xyz \
  --output-dir /srv/jobs/123/work \
  --web-output-dir /srv/jobs/123/outputs \
  --config /srv/pdb2net/config.json \
  --headless
```

On a cache miss, `--source-dir ... --populate-missing` performs a targeted
lookup and precomputes only the requested entry. It does not recursively scan
the complete archive for every request.

## Archive layout and compatibility

The recommended managed mirror has one authoritative representation: the
current Core archive mmCIF gzip under the Extended PDB ID. With the archive
generation itself as `--source-dir`, the exact layout is:

```text
archive/current/
└── entries/
    └── ab/
        └── pdb_00001abc/
            └── structures/
                └── pdb_00001abc.cif.gz
```

Do not mirror a second legacy PDB representation, uncompressed duplicate, or
enriched NextGen copy for the PDB-ID service. Enriched NextGen mmCIF remains a
supported local/upload input; it is not the server's canonical archive.

For local and migration compatibility, the targeted resolver also recognizes
plain or gzip `.cif`/`.mmcif` Extended-ID names at the archive root or shard,
classic-ID `.cif`/`.mmcif`/`.pdb`/`.ent` names, and the established wwPDB
`divided/{mmCIF,pdb}/<shard>` variants. The two-character shard is derived from
the final four legacy-compatible characters (`ab` for `pdb_00001abc`).

If no exact candidate exists, lookup fails. If several supported formats or
locations match, lookup fails with `PDB_ARCHIVE_ENTRY_AMBIGUOUS` instead of
choosing an arbitrary representation or revision. The bulk `precompute` command
can discover nested layouts with `--recursive`; multiple sources resolving to
one PDB ID are likewise rejected.

## Profile invalidation and promotion

The geometry-profile fingerprint covers artifact/scientific-pipeline semantics,
Gemmi/parser semantics, distance cutoffs, contact filters, and first-model
policy. Source-content changes invalidate the entire entry. A separate
annotation-profile fingerprint covers embedded-SIFTS use, BLAST/DIAMOND search
policy, and `reference_manifest_id`. Schema 2 requires a non-empty reference
manifest. The value must identify more than the reference filenames: make it a
content address or immutable release-manifest key that records the exact PDB
FASTA, SIFTS, Swiss-Prot and optional UniRef90 inputs, the Core/worker build, and
the concrete `blastp` and optional `diamond` versions. For example,
`pdb-sifts-swissprot-2026-07-worker-a1b2c3` may refer to an immutable external
manifest containing their checksums and versions. Core cannot infer external
binary versions from an operator-chosen label, so never reuse a manifest ID
after any of those inputs changes.

For a reference, policy, or scientific-core update:

1. provision the new references and assign a new `reference_manifest_id`;
2. stage a copy/snapshot of the store; a changed geometry profile creates a new
   namespace, while an annotation-only change refreshes annotation sections and
   can retain compatible geometry;
3. precompute a small representative set and compare normalized chain,
   annotation, and edge data with raw runs;
4. complete the desired background population and inspect the report for
   failures;
5. atomically switch the worker's store mount/config to the validated store
   (read-only for cache-only service, writable only when lazy population is
   enabled) and keep the previous generation available for rollback.

Do not edit profile manifests or move entry files between geometry-profile namespaces.
Cache entries are treated as untrusted input and are rejected on schema,
profile, ID, endpoint, count, gzip/JSON-size, or symlink violations.

## Daily mirror and promotion workflow

Treat the mirror and store as versioned generations; do not update the live
worker view in place:

1. synchronize the upstream wwPDB mirror into a staging snapshot and finish the
   sync before starting precompute;
2. select the immutable reference/toolchain/worker manifest ID described above;
3. prepare a setgid staging-store generation with runtime-compatible UID/GID;
4. when the scientific profile is unchanged, seed staging from the prior
   generation and run `pdb2net precompute --recursive`; unchanged source
   SHA-256 values become cache hits, while new or changed sources are published
   atomically;
5. when only the reference/annotation profile changes, refresh annotations in
   staging and verify geometry reuse counts; when the geometry profile changes,
   populate its new namespace rather than moving entries into it;
6. require a zero-failure precompute report, validate representative raw versus
   cached chain/edge results, run single- and multi-ID assemble smokes, and
   verify permissions from the worker identity;
7. compare the new mirror inventory with the served PDB-ID catalog. IDs removed
   upstream are obsolete, but Core does not automatically delete their cached
   entries. Apply an explicit operator retention/removal or catalog-deny policy
   to staging; never infer deletion merely from one failed mirror sync;
8. promote the validated generation atomically outside Core, for example by
   switching the host release pointer or bind mount and recreating the worker so
   its in-container store path is canonical;
9. retain the previous mirror/store/reference/worker generation for rollback,
   and reverse the external pointer if post-promotion checks fail.

If offline precompute and lazy worker population can overlap, they must still
share compatible permissions. Per-entry locks serialize population of the same
profile/PDB artifact; they do not make an in-place mirror/reference replacement
safe or provide an obsolete-entry cleanup policy.

The store remains entirely optional. A local user can continue to run
`pdb2net run --input-dir ...` with their own PDB/mmCIF files without a store,
PDB mirror, webserver, PHP, database, or server-specific filesystem layout.
