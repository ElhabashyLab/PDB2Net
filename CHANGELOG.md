# Changelog

The RC2–RC4 entries below record unreleased development history. Their older
server-interface snapshots, output contracts, profile namespaces, and rollout
ideas are superseded by the RC5 clean-release baseline and must not be treated
as current compatibility requirements.

## v0.2.0-rc5 - Unreleased

- Retain fully modified supported polymer chains and use unambiguous polymer
  types when residue-name inference is unknown (F08). Preserve determinate
  name-based types and reference precedence; parser/annotation fingerprints v6
  distinguish corrected stores.
- Correct BLAST alignment-span coverage when `qcovs` is unavailable (F09);
  Swiss-Prot selection policy v3 separates prior search-cache and store results.
  Apply the same X-fraction eligibility in preflight and actual search (F10).
- Correct the independent combined-network test oracle to accept identity
  spanning trees, and add negative search-threshold and compatibility cases.
- Run installed-wheel `run`, `precompute` and read-only `assemble` in CI,
  outside the checkout with runtime dependencies and without Cytoscape, checking
  fixed real-fixture contacts, annotations, lengths and protein aggregation.
- Include 23 frozen real structures and authentic mini-references in the
  default test suite, with independent scientific expectations and strict
  input hashes.
- Correct F06/F07 lengths using existing residue polymer roles and types:
  include modified nucleotides and exclude nonpolymer amino acids. Keep legacy
  name-based calls, contact rules, sequence annotation and public schemas;
  scientific store profile v4 distinguishes corrected lengths. The three
  real-file length regressions now pass normally.
- Replace the unreleased capability schema 2 semantic snapshot with the small,
  configuration-free capability schema 3 document.
- Publish CLI contract 2 and web-output contract 2.0.
- Remove runtime lazy-population options; `assemble` is a read-only store
  consumer and `precompute` is the only writer.
- Publish the schema-3 one-profile store with manifest-last offline builds and
  bounded read-only assembly.
- Make live Cytoscape support an optional package extra and add Ruff plus
  installed-wheel checks to the Python 3.11/3.12 quality matrix.
- Remove the unused experimental external Java/Prefuse layout backend; headless
  exports continue to use the deterministic built-in Python layout.
- Build and validate a fresh configuration snapshot per CLI invocation, apply
  path normalization after all override layers, and propagate the snapshot to
  spawned parser workers without mutating `PDB2NET_CONFIG_FILE`.
- Reuse one bounded source inventory during offline precompute and isolate
  invalid sources so completed entries remain available for unpublished-store
  retries while global input limits and manifest-last publication stay strict.
- Preserve authoritative structure identities, annotation ambiguity, and
  nucleic-acid types when considering direct AlphaFold accession assignment.
- Derive fallback sequences and residue-based molecule types from observed
  polymer residues; retain modified residues.
  Annotation fingerprint v5 rejects older precomputed annotations.
- Preserve residue insertion codes in the existing detailed CSV labels, and
  reject non-Boolean `diamond.enabled` values before invoking search tools.
- Strengthen scientific CX2 comparison and standard-suite coverage of distance
  boundaries, model separation, and the checked-in real NextGen fixture.
- Continue raw runs past explicitly recognized per-file content failures,
  recording every exclusion and retaining correct subset counts. Failed runs
  retain file diagnostics without publishing results; global/integrity failures
  remain fatal. No additional full parsing pass or public contract version change.
- Restrict contacts to protein/DNA/RNA polymer heavy atoms, including modified
  residues; omit explicit zero occupancy and retain unknown occupancy silently.
  Count each atom-site pair once across eligible alternatives and export its
  minimum distance once in CSV. Thresholds and model separation are unchanged.
  Diagnose chains without contact atoms and invalid contact metadata through
  existing reports. Parser fingerprint v5 / interaction fingerprint v2 require
  offline store rebuilds; public schemas and annotation/search policy are unchanged.
- Report cross-file identity links with zero geometric contact count while
  preserving their topology, style and `python_fast` layout strength. Clarify
  the four network projections, aggregation limits and supplied-coordinate scope;
  preserve model separation and the existing restrictions on all-model outputs.
- Reject invalid or empty CX2 comparison roots so the regression tool cannot
  report PASS without comparing files. Preserve valid empty-graph comparisons.
- Document immutable reference generations per process and test fresh-process
  annotation and search-cache separation without changing reference loading.
- Add an optional per-process sequence-search timeout (disabled by default),
  coordinated cancellation and cleanup for native BLASTP/DIAMOND processes,
  and CLI exit 130 with failed raw-run reporting on user interruption.
- Report empty chain identifiers in selected models as `MISSING_CHAIN_ID`
  file exclusions, allowing valid files in a mixed raw run to finish. Keep
  whole-file rejection, model selection and strict offline store publication;
  do not invent chain names or change scientific/store fingerprints.

## v0.2.0-rc4 - 2026-08-26

- Add `pdb2net --version` and a configuration-free
  `pdb2net capabilities --json` handshake for server integrations.
- Centralize CLI, web-output, web-config, capability, and precomputed-artifact
  contract constants so producers and compatibility probes cannot drift through
  duplicated literals.
- Stream detailed atom-interaction CSVs with cumulative row, byte, and
  free-space limits so web jobs cannot exhaust memory or staging storage before
  artifact validation.
- Publish the existing single-block mmCIF and content-identity conflict policy
  through the stable server interface.

## v0.2.0-rc3 - 2026-07-15

- Add an optional, portable precomputed graph store for PDB-ID based runs while
  preserving the standalone local `pdb2net run --input-dir` workflow.
- Add `pdb2net precompute` to cache compact annotated chain metadata and
  standard-profile chain-pair edges without retaining atom coordinates or
  detailed atom-pair tables.
- Add `pdb2net assemble` to create per-PDB and combined CX2 networks from
  validated cache entries with normal output-contract 1.1 success/failure
  summaries.
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
