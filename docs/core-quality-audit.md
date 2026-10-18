# Validation summary

The scientific pipeline is verified through small independent cases, frozen
real structures and installed-package integration. These checks complement one
another: raw/store agreement alone cannot reveal a shared scientific error.

## Verified development baseline

The current development tree completed 664 tests under Python 3.12/Linux with
no skipped or expected-failure cases. Three occurrences of the known
`datetime.utcnow` deprecation warning remain. Ruff passed. The isolated wheel
check passed with its runtime dependencies and without the Cytoscape extra.
These are local results; remote CI and other platform acceptance are listed in
[current status](core-quality-plan.md).

| Area | Evidence |
| --- | --- |
| Input handling | PDB/mmCIF, gzip, content identity, merged chain parts, insertion codes, model separation and typed file exclusions |
| Annotation | Embedded/external SIFTS precedence and ambiguity, direct AlphaFold safeguards, polymer-only sequences and real mini-database fallback |
| Contact calculations | Independent Euclidean counts, inclusive cutoffs, rigid transformations, occupancy, alternative-site unions and minimum-distance witnesses |
| Detailed CSV | All expected atom pairs in the 1BNA case, residue identity, rounding and cumulative row/byte limits |
| Lengths and compaction | Modified nucleotides and nonpolymer amino-acid counterexamples, raw/compact export and store regression |
| Networks | Four declared projections, hand-counted aggregation and independent identity-tree invariants |
| Search lifecycle | Positive hits and controlled misses, quality boundaries, format fallback, cache separation, timeout and cancellation |
| Store and package | Profile rejection, manifest-last publication, immutable assembly and real installed-console calculations |

## Regression evidence

The fixed dataset contains 23 structure files and authentic mini-references.
Hashes and provenance are in the [fixture documentation](../tests/fixtures/core_reference/README.md).
Expected results are independently derived rather than adopted from existing
Core output. The full suite includes 44 real-input integration cases.

Recent counterexamples protect fully modified RNA/DNA chains, BLAST coverage
without `qcovs`, and candidate sequences excluded by the X-fraction rule.
Removing the BLAST coverage threshold in a controlled in-memory mutation makes
six focused tests fail. The network oracle accepts a valid two-structure,
four-chain identity tree and rejects missing, extra and nonzero bridges.

The installed-wheel check processes 4HHB and 1BNA through the real console
command, checks fixed annotations, lengths, chain contact counts, protein sums
and raw CSV row counts, then assembles after removing the source/reference
copies. The store is read-only and its file hashes remain unchanged.

## Reproduction and limits

Use [Development](development.md) for environment setup and exact commands:

```bash
.venv/bin/python -m pytest -q -ra
.venv/bin/ruff check pdb2net tests scripts
```

The same guide describes `scripts/check_installed_pipeline.py`, which must run
with the separately installed wheel outside the checkout.

Passing tests do not prove every branch or every possible structure correct.
Live-Cytoscape behavior is partly represented by stubs and fallback tests;
this is not a complete live acceptance. Large archive workloads, Windows/macOS
and external deployment state are outside this local validation.
