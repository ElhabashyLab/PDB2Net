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

## Optional real-file integration tests

The normal test command runs the repository's small tests and explicitly skips
`tests/integration/`. This is also the default CI behavior. Enable the separate
32-case real-file suite with:

```bash
.venv/bin/python -m pytest tests/integration --run-integration \
  --integration-data /path/to/structure-fixtures \
  --integration-references /path/to/reference-fixtures -q -ra
```

To run both suites together, replace `tests/integration` with `tests` in that
command. Include the test path so pytest discovers these project-specific options.
Explicit integration selection without `--run-integration`, missing options,
missing/unreadable fixtures, or unavailable `blastp`/`makeblastdb` fails with a
clear error. The tests never download structures or reference datasets.

The structure directory must contain these fixed fixtures. Files may also be
in its `review-ergaenzungen/` subdirectory, so the existing reviewed dataset can
be reused directly:

| Group | Required files |
| --- | --- |
| Existing structures | `4hhb.pdb`, `9rat.pdb`, `1TUP.pdb`, `1tup.cif`, `1a34.cif`, `124d.cif`, `170d.cif`, `5rlu.pdb`, `6w41.cif`, `6m17.cif` |
| Identity conflict fixture | `123456.cif`: the deliberately modified 9JR2 mmCIF with `data_123456` and `_entry.id 123456`, while its PDB database reference remains `9JR2` |
| Multiple models and modified residues | `1a1t.pdb`, `1a1t.cif`, `1a8o.pdb`, `1a8o.cif` |
| Regular and full NextGen | `8aly.cif`, `pdb_00008aly_xyz-enrich.cif.gz` |
| AlphaFold DB | `AF-P69905-F1-model_v6.pdb`, `AF-P69905-F1-model_v6.cif` |

Use matching PDB/mmCIF snapshots and retain the files' provenance and checksums.
This fixed suite uses individual structures up to 12 MB. It does not scan the
rest of the supplied directory or run a full archive batch.

The reference directory contains the three small, matching reference subsets:
`pdb_seqres.txt`, `sifts.tsv`, and `swissprot.fasta` (at most 5 MB each). Reuse
the prepared review subsets; the full reference corpus is not needed here.
Old review `config.json`, `inventory.json`, output folders, and prebuilt BLAST
databases are not required. The fixture builds its own small BLAST database
with `makeblastdb`, records reference hashes/tool versions, and creates fresh
headless configuration from the distributable defaults.

The fixture's generated database, BLAST caches, configurations, and run outputs
stay under pytest's temporary directory. Individual CLI cases retain `command.json`,
`stdout.log`, and `stderr.log`; `--junitxml /path/to/results.xml` can save a test
summary. Keep these artifacts and the external datasets out of Git. A future
CI job for this suite also needs these fixtures provisioned separately.

The suite covers real parser/model handling, NextGen and AlphaFold annotation,
free filenames, failure cases, CSV/CX2 export, batching, and raw versus read-only
store assembly. CX2 comparisons check scientific and style semantics, allowing
the input modes' different source-file provenance. They do not assert byte
equality or replace an independently accepted scientific goldstandard.

## Local Configuration

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
