# Current verification status and remaining work

The current development tree passes 664 tests under Python 3.12/Linux, including
44 frozen real-input integration cases. Ruff and the isolated installed-wheel
pipeline check pass. These results describe the local development tree, not a
published release or a completed remote CI run. See the
[validation summary](core-quality-audit.md) for scope and reproduction.

## Completed

- Scientific regression uses frozen structure/reference fixtures and independent
  expected sequences, lengths, contacts and network projections.
- Polymer-role/type lengths exclude nonpolymers and include modified residues.
  Supported fully modified chains are retained; consistent polymer types fill
  otherwise unknown name-based classifications.
- BLAST alignment-span coverage is corrected, and readiness checks share the
  actual query-quality rule. Negative selection boundaries are tested.
- Identity-network verification accepts the intended spanning tree and rejects
  missing, extra or nonzero identity edges.
- CI includes installed-wheel `run`, `precompute` and read-only `assemble`
  outside the checkout without the Cytoscape extra.

## Remaining work

| Item | Scope |
| --- | --- |
| Windows source labels | Protein `source_files` currently splits only on `/`; native Windows paths can retain a full path instead of a basename. Contact weights are unaffected. |
| Platform acceptance | Verify supported Windows/macOS installation, paths, external tools and cancellation. Linux checks do not establish those results. |
| Remote CI and publication | Confirm CI for the eventual submitted revision. Release, tag, push and deployment are separate authorized actions. |
| Optional test maintenance | Reduce repeated full CLI runs for already-covered polymer lengths where this preserves focused regression diagnostics. |

There is no active archive-processing, store-migration or deployment task.
This status file records work; it does not grant permission to perform it.
