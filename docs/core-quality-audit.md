# Validation of the scientific baseline

The baseline completed 617 tests under Python 3.12/Linux with no skipped or
expected-failure cases and two occurrences of the known datetime.utcnow
warning. Ruff passed. These are local results, not remote CI or platform-wide
acceptance.

The suite includes 44 real-input integration cases. The frozen dataset contains
23 structures and authentic mini-references with independent expected sequences,
lengths, contacts and network projections. Sources and checksums are documented
in the [fixture README](../tests/fixtures/core_reference/README.md).

Small cases independently check inclusive distance cutoffs, rigid transformations,
alternative-site unions, occupancy and insertion codes. The 1BNA CSV is checked
against every expected atom pair. Raw/store regressions check compact lengths,
profile compatibility and read-only assembly. Real mini-database searches and
bounded fake processes cover fallback and search lifecycle behavior.

Passing checks do not prove every intermediate decision correct. Remaining
modified-polymer, BLAST-selection and network-oracle findings are tracked in
[current status](core-quality-plan.md). Scientific installed-wheel calculations,
Windows/macOS, live Cytoscape and large archive workloads still need their
respective acceptance checks.

Use [Development](development.md) for the reproducible test commands. Raw/store
agreement alone is not an independent scientific oracle.
