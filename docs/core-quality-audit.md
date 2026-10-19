# Validation summary

The scientific pipeline is verified through small independent cases, frozen
real structures and installed-package integration. These checks complement one
another: raw/prepared agreement alone cannot reveal a shared scientific error.

## Verified development baseline

Development baseline `5ca056c` completed 664 tests under Python 3.11 and
3.12/Linux in [CI](https://github.com/ElhabashyLab/PDB2Net/actions/runs/39139887748),
with no skipped or expected-failure cases. Python 3.12 reported three occurrences
of the known `datetime.utcnow` deprecation warning. Ruff and the isolated wheel
check passed with runtime dependencies and without the Cytoscape extra.

The subsequent source-label correction passed 46 focused local protein/export,
headless and raw/store cases. Seven new cases exercise native Windows and POSIX
path rules through actual CX2 output, including tooltips and nucleic-acid neighbor
metadata. Three Windows cases failed before correction. This is not a full
Windows installation or process-lifecycle acceptance; see [current status](core-quality-plan.md).

| Area | Evidence |
| --- | --- |
| Input handling | PDB/mmCIF, gzip, content identity, merged chain parts, insertion codes, model separation and typed file exclusions |
| Annotation | Embedded/external SIFTS precedence and ambiguity, direct AlphaFold safeguards, polymer-only sequences and real mini-database fallback |
| Contact calculations | Independent Euclidean counts, inclusive cutoffs, rigid transformations, occupancy, alternative-site unions and minimum-distance witnesses |
| Detailed CSV | All expected atom pairs in the 1BNA case, residue identity, rounding and cumulative row/byte limits |
| Lengths and compaction | Modified nucleotides and nonpolymer amino-acid counterexamples, raw/compact export and store regression |
| Networks | Four declared projections, hand-counted aggregation and independent identity-tree invariants |
| Search lifecycle | Positive hits and controlled misses, quality boundaries, format fallback, cache separation, timeout and cancellation |
| Prepared inputs and package | Profile rejection, reference-independent facts, current annotation, atom-contact validation and installed-console/API calculations |

## Prepared-input transition verification

The local 0.2.0rc6 transition passed **653 tests on Python 3.12/Linux**, Ruff,
and the installed-package check outside the source checkout without Cytoscape.
The package contains the prepared API and no former `pdb2net.precomputed` module.
Independent reviews covered raw/prepared scientific equivalence, fixed expected
values, source provenance, contact validation, reference changes and package
contracts. This is local evidence; the published baseline above remains the
last remote CI result recorded here.

The reduction from the earlier test count reflects removal of Core-owned store
publication tests and porting their scientific cases to the new API. External
persistence tests belong to their adapters. No fixed scientific oracle was
replaced with the implementation's own output.

## Regression evidence

The fixed dataset contains 23 structure files and authentic mini-references.
Hashes and provenance are in the [fixture documentation](../tests/fixtures/core_reference/README.md).
Expected results are independently derived rather than adopted from existing
Core output. The published baseline included 44 real-input integration cases; the prepared
transition retains their independent expectations and adds boundary checks.

Recent counterexamples protect fully modified RNA/DNA chains, BLAST coverage
without `qcovs`, and candidate sequences excluded by the X-fraction rule.
Removing the BLAST coverage threshold in a controlled in-memory mutation makes
six focused tests fail. The network oracle accepts a valid two-structure,
four-chain identity tree and rejects missing, extra and nonzero bridges.

The installed-wheel check processes 4HHB and 1BNA through the real console
command, checks fixed annotations, lengths, chain contact counts, protein sums
and raw CSV row counts, then analyzes prepared facts after removing coordinate
source copies. Reference data remains available for fresh annotation. Facts
are unchanged and both paths produce byte-equal CSVs.

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
