# Current verification status and remaining work

The published prepared-input revision `61888ec` passed 653 tests under Python
3.11 and 3.12/Linux, plus Ruff and the isolated installed-wheel check in
[CI](https://github.com/ElhabashyLab/PDB2Net/actions/runs/39248747446).
This includes the Windows source-label correction, exercised through Windows
path semantics on Linux rather than a complete Windows platform acceptance.
See the [validation summary](core-quality-audit.md).

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
- CI includes installed-wheel `run` and the prepared-input API outside the
  checkout without the Cytoscape extra.
- Protein source labels use native path rules, including Windows drive and UNC
  paths, while preserving POSIX filenames and prepared-input source labels.

## Remaining work

The prepared-input transition is implemented: Core 0.2.0rc6 provides
capabilities 4, CLI 3, outputs 3.0 and facts 1. The former Core store commands
are removed; adapters own persistence while the shared API prepares, validates,
annotates and exports. The final Core suite passed 653 tests on Python 3.12;
Ruff and the installed-package raw/prepared check passed without Cytoscape.
Independent step reviews found and resolved the pair-iteration memory issue;
iteration now retains memory proportional to chains. Scientific fixtures and
independent expected values were retained when replacing old store tests.

Existing real fixtures also passed a bounded preparation/streaming measurement,
including the 20-chain full 8ALY structure and its enriched counterpart. Those
results do not establish whole-archive capacity. The final independent review
found no remaining blocker. The working branch is published and remote CI passed;
it has not been merged into `main` or released.

| Item | Scope |
| --- | --- |
| Platform acceptance | Verify supported Windows/macOS installation, paths, external tools and cancellation. Linux checks do not establish those results. |
| Optional test maintenance | Reduce repeated full CLI runs for already-covered polymer lengths where this preserves focused regression diagnostics. |

No full archive run or production deployment is part of this local transition.
This status file records work; it does not grant permission to perform it.
