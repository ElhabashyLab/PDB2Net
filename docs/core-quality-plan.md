# Current verification status and remaining work

## Precompute audit corrections

Four subsequent audit findings are corrected locally:

- Structural sequences retain SEC as U and PYL as O. Existing search masking
  of O remains unchanged. Prepared sequence/type policy is now
  `polymer-residues-v2`.
- Contact decisions at numerical cutoff boundaries use exact decimal squared
  distances, with bounded KD-tree candidate selection. Contact policy is now
  `pdb2net-distance-classification-v3`; just-outside pairs remain excluded.
- Same-accession, equal-bitscore HSP reduction prefers the lower E-value.
  Swiss-Prot search selection policy 4 invalidates earlier cached decisions.
- Combined protein metadata uses only chains in its linked component, including
  source labels, names, embedded annotation evidence and tooltips.

Local Python 3.12 verification passed 697 tests (26 additional cases), Ruff,
and the isolated installed-wheel raw/prepared pipeline check without Cytoscape.
The suite reports three existing `datetime.utcnow()` deprecation warnings.
Independent code and documentation reviews found no remaining blockers;
the independent numerical review additionally checked 640 boundary cases against
direct exact arithmetic, including alternative positions and large translations.
Frozen scientific fixtures and their expected results were not changed.

The profile changes deliberately reject all earlier-profile prepared facts,
including internally consistent old hashes. Adapters need a new store generation
prepared from sources; preserve the previous Core/store pair until transition
validation completes. The facts schema and output contracts remain unchanged.
These local corrections are not yet published or deployed; the archive build
must wait for the reviewed Core/store transition. No automated migration is added.

## Earlier verification baseline

The published prepared-input revision `61888ec` passed 653 tests under Python
3.11 and 3.12/Linux, plus Ruff and the isolated installed-wheel check in
[CI](https://github.com/ElhabashyLab/PDB2Net/actions/runs/39248747446).
This includes the Windows source-label correction, exercised through Windows
path semantics on Linux rather than a complete Windows platform acceptance.
See the [validation summary](core-quality-audit.md).

## Completed

- Polymer-to-author chain mapping no longer consumes the 200,000-row embedded
  SIFTS segment budget. Structural mapping uses column views and a separate
  20-million-row ceiling (the existing atom-mapping bound); true segment limits
  and mapping semantics remain unchanged. Regression covers large polymer
  tables with and without embedded annotations. Parser/store and annotation
  fingerprints remain unchanged because previously successful inputs retain
  their facts and annotation results; previously rejected inputs have no
  successful prepared revision to invalidate. Local Python 3.12 verification
  passed 671 tests, Ruff and the isolated installed-wheel raw/prepared check.
  Independent review found no remaining blockers, including preservation of
  CIF token interpretation and existing tag-suffix case handling.
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
