# Reference generations and search lifecycle

## One immutable reference generation per process

SIFTS, PDB FASTA, UniProt FASTA and search databases form a coherent reference
generation. They must not be replaced or redirected during a Core process.
Loaders retain path-keyed reference mappings in memory; `activate_config` and a
new `reference_manifest_id` do not promise to reload those mappings.

Prepare updated references in a new directory with a new reference identifier,
then start new Core processes against that generation. Existing processes may
finish with their original generation. Notebook users must restart the kernel
when switching generations. This is a usage contract, not a per-access content
hash check or an automatic invalidation mechanism.

RAM reference mappings, SQLite sequence-search results and prepared structure facts
are separate layers. Search signatures distinguish relevant database/policy
changes; geometry profiles distinguish structure rules independently of reference generations.
Deleting a search cache does not refresh a process's reference mappings.

## Optional search timeout

`sequence_search_timeout_seconds` is shared by BLASTP and DIAMOND. Its default
is `null`, leaving each search without a Core time limit. A positive integer
limits each external process invocation. The value follows normal configuration
layering and is included in the runtime configuration snapshot.

The limit is not per file or per whole job. A BLAST output-format compatibility
retry is a separate invocation with its own limit. A caller may impose an
additional overall job deadline. No fixed timeout suits every database,
hardware configuration and query batch.

A timeout is an operational failure, not a biological non-hit. It must not
create a negative cache entry or trigger DIAMOND as a consequence of a BLAST
failure. Raw processing fails; prepared analysis also fails without publishing successful results. Already completed valid cache results may remain.

## Cancellation and boundaries

On timeout or interruption, pending chunks do not start. Running direct native
search processes are terminated and awaited before temporary query files are
released. Concurrent workers share cancellation state; cancellation does not
wait for result delivery in input order. CLI interruption returns 130 and does
not produce a successful run report.

Search children retain their caller's process-group membership. The Core
guarantee covers direct native BLASTP/DIAMOND programs and their threads, not
arbitrary wrappers that detach new sessions or independent process trees.
An external adapter owns its own job supervision and cleanup policy.

Timeout configuration does not alter the scientific search/geometry profile.
Changing it alone requires neither a store rebuild nor a search-cache purge.
Scientific policy changes are handled separately by the existing fingerprints.

## Verification

`tests/test_search_runtime.py` uses bounded fake processes to check timeouts,
launch/cancellation races, interruption, cleanup and error propagation.
`tests/test_reference_data.py` checks reference-generation behavior. The example
`docs/examples/reference_generation_probe.py` reproduces the difference between
replacing a file in one process and reading a new generation in a fresh process.

Signal/process-group verification is limited to the tested Linux environment.
Windows/macOS behavior requires its own platform acceptance; see
[current status](core-quality-plan.md).
