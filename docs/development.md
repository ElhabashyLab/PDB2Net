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
- whether configured local reference/input/output paths exist

Missing Cytoscape is acceptable for headless checks. Missing BLAST tools or reference files only block the BLAST-backed annotation path, not pure code editing.

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
- `pdb2net --help` and `pdb2net run --help` work after installation

## Full Pipeline Requirements

The full pipeline needs:

- a folder with `.pdb`, `.cif`, or `.mmcif` inputs
- `pdb_seqres.txt`
- `pdb_chain_uniprot.tsv`
- `uniprot_sprot.fasta`
- a BLAST database created from `uniprot_sprot.fasta`
- `blastp` and `makeblastdb` available in `PATH` or configured explicitly

Run the current script entry point with:

```bash
python3 -m pdb2net.main
```

Do not run full batch jobs on large datasets during routine development. Use tiny fixtures or a small local input folder first.

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
