# Current verification status and remaining work

The current Core 0.2.0rc6 contract is capabilities 4, CLI 3, output 3.1 and facts 1,
including `mixed_structure_analysis`. The latest full local Linux verification below
passed 746 tests; native Windows passed 732 with 14 explicit platform/privilege
skips. Earlier counts and publication records describe their named
milestones; they are not evidence of publication or deployment of the current pair.

## Native Windows file I/O corrections

Windows path-stat and handle-stat calls can disagree on `ctime` for unchanged
files. Only that cross-API Windows comparison omits `ctime`; full path-to-path
checks, file identity/size/mtime checks, existing digests and POSIX behavior
remain intact. Public summary publication retains file fsync and replacement
on both platforms; only Windows omits unsupported parent-directory fsync.
Scientific rules, fingerprints and public contracts are unchanged.

The symlink tests skip only on Windows error 1314, and home-directory/path-label
expectations now account for Windows conventions. Seventeen new cases cover
handle identity mismatches and summary synchronization ordering/error propagation.
Python 3.12.3 on Linux passed 746 tests without skips. Native Windows x64 with
Python 3.12.10 passed 732 tests, with nine existing POSIX-only skips and five
missing-symlink-privilege skips. Both suites retain three existing `utcnow()`
deprecation warnings. Isolated installed-wheel raw/prepared checks passed on both
platforms without the Cytoscape extra. Ruff also passed. Native Python 3.11 and
macOS were not rerun for this correction.
Independent review repeated 159 focused tests, checked the native test reports
and installed product files, and confirmed narrow privilege skips and preserved
POSIX directory-error propagation without open findings.

## Bounded negative-occupancy acceptance

Finite negative contact-atom positions are now omitted when they make up at most
2% of the otherwise eligible nonzero/unknown positions in each selected model.
The gate shares the existing contact-eligibility pass; counts and CSV share its
validated selection. Nonfinite/nonnumeric values, excessive negative fractions,
zero-only structures and unsupported polymer types retain their exclusions.
A bounded `NEGATIVE_OCCUPANCY_ATOMS_SKIPPED` warning records the omission.

The acceptance extension preserves the successful-input geometry profile and
facts schema. Twenty existing successful real fixtures produced byte-identical
facts and identical contact streams before/after the change. All 105 audited
sources were checked locally: 29 negative-occupancy cases now produce validated
facts/contact bundles; the other 76 retain `NO_VALID_CHAINS`. Original sources
and frozen scientific expectations were not changed.

Python 3.12 verification passed 729 tests, including 17 new boundary/reporting/
raw-prepared parity cases, with three existing `datetime.utcnow()` warnings and
no skips. Independent review repeated the full suite and verified all 29 local
bundles. Ruff and the isolated installed-wheel raw/prepared check also passed.
This is local verification, not publication or deployment; update runtime and
offline preparation together before publishing newly accepted entries.

## Earlier profile-changing audit corrections

Four audit findings established the current geometry/annotation baseline:

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

Verification at that milestone passed 697 tests on Python 3.12 (26 additional cases), Ruff,
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
This requirement applies when upgrading from a profile predating those corrections.
The later mixed-input and bounded negative-occupancy changes preserve this profile;
compatible successful facts do not need rebuilding. No automated migration is added.

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

The prepared-input transition and mixed-input orchestration are implemented.
Adapters own persistence while the shared API prepares, validates, annotates and
exports. The earlier 653-test milestone resolved the pair-iteration memory issue;
iteration retains memory proportional to chains. Scientific fixtures and independent
expected values were preserved. The current 729-test verification is recorded above.

Existing real fixtures also passed a bounded preparation/streaming measurement,
including the 20-chain full 8ALY structure and its enriched counterpart. Those
results do not establish whole-archive capacity. Publication and deployment of the
current paired adapter/Core revisions require their own verification; the earlier
published baseline and its CI result are documented separately above.

| Item | Scope |
| --- | --- |
| Platform acceptance | Verify supported Windows/macOS installation, paths, external tools and cancellation. Linux checks do not establish those results. |
| Optional test maintenance | Reduce repeated full CLI runs for already-covered polymer lengths where this preserves focused regression diagnostics. |

No full archive run or production deployment is part of this local transition.
This status file records work; it does not grant permission to perform it.
