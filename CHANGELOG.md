# Changelog

## v0.2.0-rc4 - Unreleased

- Add `pdb2net --version` and a configuration-free
  `pdb2net capabilities --json` handshake for server integrations.
- Centralize CLI, web-output, web-config, capability, and precomputed-artifact
  contract constants so producers and compatibility probes cannot drift through
  duplicated literals.

## v0.2.0-rc3 - 2026-07-15

- Add an optional, portable precomputed graph store for PDB-ID based runs while
  preserving the standalone local `pdb2net run --input-dir` workflow.
- Add `pdb2net precompute` to cache compact annotated chain metadata and
  standard-profile chain-pair edges without retaining atom coordinates or
  detailed atom-pair tables.
- Add `pdb2net assemble` to create per-PDB and combined CX2 networks from
  validated cache entries, including targeted lazy population from a local PDB
  archive and normal output-contract 1.1 success/failure summaries.
- Namespace artifacts by a scientific profile fingerprint so reference,
  annotation-policy, cutoff, filter, and schema changes never silently reuse or
  overwrite incompatible entries.
- Support canonical plain and gzip-compressed PDB/mmCIF archive names.
- Document shared-group/setgid store permissions, fail-closed cache ceilings,
  versioned daily mirror promotion and rollback, explicit obsolete-ID handling,
  and immutable reference/toolchain manifest requirements.

## v0.2.0-rc2 - 2026-07-15

- Add bounded large-job processing controls and additive output-contract 1.1
  metadata for resource use, annotation sources, warnings, and skipped outputs.
- Batch optional DIAMOND queries against small or externally provisioned
  UniRef90-compatible databases while preserving Swiss-Prot-first annotation.
- Add a real, tiny UniRef90-format DIAMOND fixture for deterministic integration
  tests without downloading or committing large reference datasets.
- Keep operational BLAST/DIAMOND failures out of the negative cache, invalidate
  cached results when search policies change, and reject non-UniProtKB cluster
  representatives for canonical `uniprot_id` assignment.
- Reduce the per-job reference-memory floor by reusing cached mappings and
  retaining PDB SEQRES headers only for structures in the current batch.
- Add Python 3.11/3.12 CI and installed-package verification.

## v0.2.0-rc1 - 2026-06-15

- Backend-ready CLI and optional web-output adapter with stable `summary.json`,
  `networks/`, and `interactions/` outputs.
- Headless CX2 export improvements for automated runs without Cytoscape UI.
- Release hardening for versioned output contracts, failed-run summaries,
  CLI error exit codes, reference-data preflight checks, and safer KDTree
  caching across multiple structures processed in one run.
- Optional DIAMOND/UniRef90 fallback scaffolding for broader annotation without
  bundling or downloading large external databases.
