# Development Setup

This project can be edited without the full reference datasets, but meaningful runtime checks need local Python dependencies and, for the full annotation pipeline, BLAST reference data.

## Local Python Environment

Create a virtual environment and install the package in editable mode:

```bash
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install --upgrade pip setuptools wheel
python3 -m pip install -e ".[dev]"
```

After the venv exists, run development commands through the venv interpreter
unless you have explicitly activated it:

```bash
.venv/bin/python -m pytest
.venv/bin/ruff check pdb2net tests scripts
.venv/bin/python scripts/check_environment.py
```

This matters because the system `python3` may point to a different Python
version without development dependencies. For example, a shell can have
`python3` from Conda while this project uses `.venv/bin/python` with Python
3.11 or 3.12 and `pytest` installed.

If editable install is not needed, installing the pinned runtime dependencies is enough:

```bash
python3 -m pip install -r requirements.txt
```

If editable install needs to be refreshed in an offline or sandboxed
environment, avoid build isolation so pip uses the build tooling already
installed in the venv:

```bash
.venv/bin/python -m pip install --no-build-isolation --no-deps -e .
```

## Small real fixture in the default suite

The normal suite runs the checked-in, provenance-documented 8ALY NextGen excerpt
through the real headless CLI, offline precompute, and read-only assembly. It
asserts known author-chain mapping, UniProt/Pfam segment ranges and CX2 contents.
Tiny synthetic FASTA/TSV lookups are created under pytest's temporary directory;
this test needs neither an external dataset nor a BLAST database or Cytoscape.
It runs in the existing Python 3.11/3.12 CI jobs. The excerpt has two coordinate
residues, so it complements rather than replaces the full-structure suite below.

## Frozen real-file integration tests

The default suite now includes the real-file checks using the 23 full structures
and three frozen mini-reference files in `tests/fixtures/core_reference/`:

```bash
.venv/bin/python -m pytest -q -ra
```

CI also installs the built wheel and its runtime dependencies into a separate
environment without the Cytoscape extra. From outside the checkout,
`scripts/check_installed_pipeline.py` exercises the installed console command
through `run`, `precompute`, and `assemble` using frozen 4HHB/1BNA inputs.
Both paths are checked against fixed chain identities, annotations, polymer
lengths, contact counts and protein aggregation. Raw CSV row counts are also
checked. Sources and reference copies are removed before assembly; store
permissions and file hashes enforce the read-only check.

To repeat with an already installed wheel (use a new work directory):

```bash
cd /tmp
/path/to/package-venv/bin/python -I /path/to/PDB2Net/scripts/check_installed_pipeline.py \
  --fixtures /path/to/PDB2Net/tests/fixtures/core_reference \
  --work-dir /tmp/pdb2net-installed-pipeline
```

The checker rejects source-tree imports and environments containing
`py4cytoscape`. This check needs no BLAST database or reference download.

Run only the real cases with `.venv/bin/python -m pytest tests/integration -q -ra`.
For an explicitly smaller selection, use `-m "not integration"`. Collection and
non-integration tests do not require BLAST or load the real dataset. Actual real
runs require local `blastp` and `makeblastdb`; missing tools fail clearly rather
than silently skipping these checks. CI installs both, plus DIAMOND for the
existing separate mini-UniRef90 check, under Python 3.11 and 3.12.

The [fixture README](../tests/fixtures/core_reference/README.md) lists every
file, its purpose, reference selection and scientific expectations. The dataset
contains the previous 19 reviewed inputs plus pure DNA (1BNA), modified RNA
(1EHZ), phosphorylated protein (2ERK), and ambiguous SIFTS mapping (1AO7).
These are about 26 MB of structures and 42 KB of reference text. No external
Datensatz22 directory, production references or runtime downloads are needed.
Input hashes are checked against `provenance.json`; Git preserves payload bytes
across checkouts. Older unverified source dates are explicitly identified.

`config.json` in that fixture directory freezes scientific settings. Test CLI
processes isolate local configuration and PDB2NET environment overrides. Tests
build a small BLAST database once per session, recording reference hashes and
tool versions. Databases, caches, configurations, stores and outputs remain in
pytest's temporary directory. Use `--junitxml /path/to/results.xml` to retain a
summary; CLI cases also retain command, stdout and stderr logs.

Existing external path options remain available for identical copies of this
same pinned dataset, including all four new files and the updated mini-references:

```bash
.venv/bin/python -m pytest tests/integration --run-integration \
  --integration-data /path/to/frozen/structures \
  --integration-references /path/to/frozen/references -q -ra
```

Both paths and `--run-integration` are required when overriding locations.
`--run-integration` alone selects the same built-in default data. Missing files
and mismatching hashes fail; these options do not accept arbitrary new reference
snapshots under old scientific expectations.

The existing real suite checks parsing, annotation, all four network views,
CSV, model separation, batching and raw/read-only-store output. Fixed structure
and network expectations are independently derived from source records and
blocked Euclidean calculations in `tests/reference_oracle.py`, without Core's
atom selectors or KD-trees. Reproduce candidates outside the accepted directory:

```bash
.venv/bin/python tests/reference_oracle.py --output-dir /tmp/pdb2net-reference-candidates
```

Review differences before changing accepted expectations. Raw/store agreement
alone is not independent scientific validation. Tests also check all 652 atom
pairs in the small 1BNA CSV, a real competing BLAST hit, and a controlled miss.
The F06/F07 length regressions run normally: 1EHZ has 76 nt, both 170D chains
have 12 nt, and 6M17 chains A/C have 605 aa. A dedicated raw/precompute/assemble
case verifies these lengths in compact entries and read-only output as well.
`tests/test_polymer_lengths.py` checks raw/compact export consistency, modified
and unknown polymers, excluded nonpolymers and legacy name-only calls. Length
counting does not apply atom-occupancy or contact filters.

The independent pairwise distance oracle runs in the normal test suite without
external fixtures. It checks inclusive cutoffs, atom ordering, translation, and
rotation. A small overlapping-model fixture also verifies that geometrically
close atoms from alternative models never form cross-model contacts.

The normal suite also checks polymer-only sequence/type annotation with small,
independently written PDB/mmCIF fixtures containing waters, ligands, MSE and an
unknown polymer residue. Atom retention and the existing X-fraction filter are
checked separately. The frozen real-file suite forces actual BLAST fallback
for hydrated 1A8O in both formats by removing its direct SIFTS/FASTA lookup from
temporary test references. It requires the known 70-position sequence and
P12497 assignment while retaining all 644 raw heavy atoms; Q01 contact selection
uses only eligible polymer positions from those raw records.

Q01 regressions cover missing versus zero occupancy in PDB/mmCIF, stable
metadata identity, modified polymer atoms, excluded nonpolymers, all-alternative
site unions, insertion codes and residue alternatives, minimum-distance CSV
witnesses, inclusive thresholds, the Cα minimum and row-budget prechecks.
`tests/test_contact_atoms.py` contains the small contact oracles;
`tests/test_atom_metadata.py` tests source metadata independently of selection.
The precomputed end-to-end test checks raw/store agreement with alternatives
and water, while profile tests reject the old parser/interaction versions.
The nine paired-format examples in `docs/examples/q01_contact_examples.py`
also compare the implemented union rule with hand-counted alternatives.

Q02 tests exercise real validated CX2 for all four network views, with manually
specified aggregation inputs: same/different/missing UniProt IDs, DNA/RNA,
summed weights, linked components and identity links with zero contact count.
They preserve the independent tests of geometry and annotation rather than
substituting for them. Layout regression compares identity links against the
historical raw weight 1000 through the ordinary weighted layout path. Existing
model tests remain authoritative; small PDB/mmCIF assembly-metadata cases verify
that recognizing operators does not generate additional coordinate chains.

Raw-run exclusion regressions are in `tests/test_input_exclusions.py`. They
compare a mixed-input run with a valid-only control, check complete reports and
counts, and exercise real CLI exits. Explicit counterexamples keep resource,
integrity and unexpected parser/annotation failures fatal. These tests use tiny
local fixtures and run in the standard suite without external search databases.
Missing-chain cases cover real PDB/mmCIF parsing, complete file omission despite
remaining named chains, offline manifest rejection and real CLI worker error
propagation. `tests/test_data_processor.py` also checks whitespace-only names,
the unchanged rejection before residue selection, and first/all model scope.

## Local Configuration

Search lifecycle tests use finite fake executables rather than large databases.
They check actual process exit, query cleanup, cancelled launch races, a timeout
behind an earlier `executor.map` input, SIGTERM escalation and parent-only SIGINT.
Configuration accepts null or positive integer search seconds; runtime limits
must not affect store/search fingerprints. Real signal/shebang probes are POSIX
tests; no Windows/macOS cancellation validation is implied by a Linux run.

Reference-loader caches follow a one-generation-per-process contract. New
reference generations require new paths/manifest IDs and a fresh process (or
notebook kernel); changing config alone does not invalidate RAM lookups.
`tests/test_reference_data.py` exercises fresh-process annotation, repeatability,
search-cache separation and changed store profiles with tiny local references.
This work intentionally does not change reference loading or add memory profiling.

Copy the example local config and adjust paths:

```bash
cp pdb2net/configs/config.local.example.json pdb2net/configs/config.local.json
```

`config.local.json` is ignored by Git. Keep machine-specific input folders, output folders, reference files, and BLAST database paths there.

For automated or Codex-driven checks, prefer:

```json
{
  "open_in_cytoscape": false,
  "workers": {
    "parsing": 1,
    "blast_threads": 1
  }
}
```

This avoids opening Cytoscape and keeps fixture runs small and predictable.

## Environment Check

Run:

```bash
python3 scripts/check_environment.py
```

The check reports:

- Python version
- required Python packages
- optional external commands such as `blastp`, `makeblastdb`, `diamond`, and
  `cytoscape`
- whether configured local inputs, required references, and output storage are
  usable

The check uses the same configuration loader and layer order as the CLI. It
returns a non-zero status for invalid configuration, unsupported Python,
missing or broken required imports, unusable configured inputs or output
storage, missing required reference files, and unavailable explicitly enabled
live Cytoscape support. Missing BLAST/DIAMOND fallback tools or databases and a
missing `makeblastdb` builder are reported as non-fatal readiness warnings
because inputs resolved through direct or SIFTS annotation do not use those
fallbacks.

Missing Cytoscape is acceptable when `open_in_cytoscape` is false. Missing
required references is always fatal for a normal raw-file run; missing BLAST
tools only blocks inputs that actually require BLAST-backed annotation.

`py4cytoscape` is intentionally absent from the base package. Install
`.[cytoscape]` only when testing live Cytoscape behavior; the headless suite and
base package must pass without it.

If Matplotlib warns that its default config directory is not writable, set a temporary cache directory before running imports or scripts:

```bash
export MPLCONFIGDIR=/tmp/pdb2net-matplotlib
```

## Packaging Checks

PDB2Net is intended to run as an installed package in worker containers. The
wheel must include the shared config defaults and OS-specific config files, but
must not include `pdb2net/configs/config.local.json`, which is reserved for
machine-local paths and is ignored by Git.

Check this from outside the source tree so Python imports the installed copy,
not the working directory:

```bash
PDB2NET_PY="$PWD/.venv/bin/python"
.venv/bin/python -m pip install --no-build-isolation --no-deps \
  --target /tmp/pdb2net-package-check .
cd /tmp
PYTHONPATH=/tmp/pdb2net-package-check "$PDB2NET_PY" -c 'import pdb2net, pathlib; root = pathlib.Path(pdb2net.__file__).parent; print(root); print((root / "configs/config.base.json").exists()); print((root / "configs/config.local.json").exists())'
```

Expected result:

- `configs/config.base.json` exists in the installed copy
- `configs/config.local.json` does not exist in the installed copy
- `pdb2net --version` and help for `run`, `precompute`, `assemble`, and
  `capabilities` work after installation

## Full Pipeline Requirements

The full pipeline needs:

- a folder with `.pdb`, `.cif`, or `.mmcif` inputs, optionally gzip-compressed
- `pdb_seqres.txt`
- `pdb_chain_uniprot.tsv`
- `uniprot_sprot.fasta`
- a BLAST database and `blastp` when unresolved inputs require sequence-search
  fallback; `makeblastdb` is needed to build that database

Run the current script entry point with:

```bash
python3 -m pdb2net.main
```

The legacy runner continues after individual batch failures and timeouts so it
can process the remaining inputs. Its final exit status is non-zero whenever
any batch fails, times out, or an eligible input is skipped.

Do not run full batch jobs on large datasets during routine development. Use tiny fixtures or a small local input folder first.

## Detailed interaction CSVs

With `structure_model_policy: "first"`, `Chain_A` and `Chain_B` contain the
author chain IDs. With `"all"`, they contain the same model-qualified IDs used
by the chain networks, for example `1ABC:model2:A`. This distinguishes contacts
from different models without changing the CSV columns.

`Residue_A` and `Residue_B` retain the author residue number and any insertion
code before the residue name: for example `10:ALA` and `10A:ALA` identify
different residues. Negative and zero author numbers are retained. Rows without
an insertion code keep their existing labels, and the CSV columns are unchanged.
This preserves residue identity without changing atom selection or distances.

## Goldstandard CX2 Regression Check

Use the accepted CX2 output directory as a semantic goldstandard. The check
compares graph semantics, annotations, visual style semantics, and layout
metrics without relying on byte-for-byte CX2 equality.

Compare an existing output directory:

```bash
python3 scripts/run_goldstandard_check.py \
  --actual /path/to/pdb2net_outputs/2026-06-13_14-44-59 \
  --expected /path/to/goldstandard/expected
```

Run the configured headless pipeline first, then compare the newest output:

```bash
python3 scripts/run_goldstandard_check.py \
  --expected /path/to/goldstandard/expected
```

Reports are written to `/tmp/pdb2net-goldstandard` by default. A `FAIL` status
means graph semantics, annotations, required style semantics, or expected files
changed. A `WARN` status is acceptable for intentional layout-only differences
unless `--fail-on-warn` is used.

Comparison inputs must be existing directories, each containing at least one
regular `.cx2` file (nested directories are searched). Individual files, missing
paths and directories without CX2 files are rejected with a diagnostic and a
nonzero exit; the direct comparison CLI uses exit 2 and writes no new report.
The Python directory API raises `ValueError` for these invalid inputs. An actual
CX2 file containing an empty graph can still be compared. Check the reported
number of compared files against the intended dataset as part of acceptance.

The comparison checks explicit biological node attributes, including polymer
lengths, molecule names/types, UniProt assignments, annotation source/confidence,
embedded mapping status, `annotation_*` fields, and structure/chain provenance.
Changed or one-sided missing values fail, even when a tooltip still contains the
old value. Contact counts and other supported scientific edge values are also
compared. Duplicate numeric or semantic identities and missing edge endpoints
fail on either side rather than disappearing during graph normalization.

`source_file` and `source_files` are intentionally excluded because raw-file and
assembled networks use different file labels. Structure IDs and `source_chains`
remain significant. Export IDs, record order, tooltip formatting and timestamps
are not scientific differences; layout differences remain warnings. This tool
compares the selected semantics and graph integrity, not every CX2 schema rule.
