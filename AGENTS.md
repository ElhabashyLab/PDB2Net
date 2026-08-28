# AGENTS.md

## Project

PDB2Net extracts protein interaction networks from PDB/mmCIF structure files and exports Cytoscape-compatible CX2 networks. The pipeline combines Gemmi parsing, SIFTS/FASTA-based annotation, optional BLAST fallback, cKDTree distance searches, and headless or live Cytoscape export.

## Current Baseline And Authority

- The clean-release implementation is complete.
  `docs/history/cleanup-plan.md` is retained
  as historical implementation and decision context, not as a phase checklist
  to run again.
- Current work is governed by this file and the maintained documents in
  `docs/`. Old hardening branches and their documents are historical reference,
  not the current architecture.
- The current public contracts are package `0.2.0rc5`, capability schema `3`,
  CLI contract `2`, Web output contract `2.0`, and precomputed artifact schema
  `3`.
- Preserve the standalone raw-file pipeline, the configuration layering, and
  established scientific behavior. Do not change scientific semantics,
  contract versions, dependencies, schemas, or the Core/Web boundary without
  explicit user approval.
- `precompute` is the sole offline store writer. `assemble` is read-only and
  must never resolve an archive, download a structure, or populate a missing
  entry at runtime.
- Do not add PHP, Apache, MariaDB, Docker, Web access control, queue behavior,
  or deployment state to Core.
- Do not reintroduce removed RC4 compatibility layers, `server_interface`, lazy
  population, profile namespaces, release receipts, attestation, or rollback
  machinery.
- No large reference download, full archive run, push, tag, release, or
  deployment without explicit user approval.
- For work involving the Web adapter or target server, read the sibling
  repository's `AGENTS.md`. If
  `../pdb2net-webserver/docs/server-state.local.md` exists, read it before any
  server or deployment planning; it is a local operational handoff and never
  overrides the product contracts here.

## Runtime

- Use `python3`, preferably Python 3.11 or 3.12.
- This repository normally has a local `.venv`. Prefer `.venv/bin/python` for tests and packaging checks when it exists. The system `python3` may point to a different interpreter without dev dependencies.
- Do not assume a `python` executable exists.
- The legacy script entry point is:

```bash
python3 -m pdb2net.main
```

- Prefer the stable CLI for new automation:

```bash
python3 -m pdb2net run --input-dir /path/to/inputs --output-dir /path/to/outputs --headless
```

- When validating the installed console script, use the package entry point:

```bash
pdb2net run --input-dir /path/to/inputs --output-dir /path/to/outputs --headless
```

- Keep Cytoscape optional for automated work. Prefer `open_in_cytoscape: false` for checks, fixtures, and headless runs.
- Keep PDB2Net standalone. Do not couple core code to PHP, MySQL, Apache, or a specific webserver repository.

## Configuration

Configuration is layered in this order:

1. `pdb2net/configs/config.base.json`
2. `pdb2net/configs/config.{linux|windows|darwin}.json`
3. `pdb2net/configs/config.local.json`
4. `PDB2NET_CONFIG_FILE`
5. environment variable overrides

Use `config.local.json` or environment variables for machine-specific paths. Do not commit local reference data paths, generated outputs, BLAST databases, or large datasets.

Packaging must include the distributable defaults in `pdb2net/configs/` (`config.base.json` and OS-specific configs), but must not include `config.local.json`. Verify installed-package config availability from outside the source tree so Python cannot import the working copy by accident.

```bash
python3 -m pip install --no-build-isolation --no-deps --target /tmp/pdb2net-package-check .
cd /tmp
PYTHONPATH=/tmp/pdb2net-package-check python3 -c 'import pdb2net, pathlib; root = pathlib.Path(pdb2net.__file__).parent; print((root / "configs/config.base.json").exists()); print((root / "configs/config.local.json").exists())'
```

Expected output is `True` for `config.base.json` and `False` for `config.local.json`.

## Output Contract

- Internal runs use timestamped folders with `combined/`, `protein/`, `chain/`, `distances/`, `runtime_analysis.txt`, `manifest.json`, and `run_summary.json`.
- Web `summary.json` should include copied outputs, counts, status, warnings/errors, package/contract versions, and the relevant runtime `config` snapshot.
- When adding output types, update `run_summary.json` and the web-output adapter.
- If `--web-output-dir` is used, preserve this stable structure:

```text
outputs/
├── summary.json
├── networks/
│   └── *.cx2
└── interactions/
    └── *.csv
```

## Safety

- Do not run full batch jobs on large input folders unless explicitly requested.
- Do not download large reference datasets unless explicitly requested.
- Do not require a running Cytoscape UI for verification.
- Prefer tiny fixtures, headless CX2 export, and focused import/unit checks when validating changes.
- BLAST work requires external `blastp`, `makeblastdb`, a Swiss-Prot FASTA, and a built BLAST database.
- Optional DIAMOND/UniRef90 fallback work requires external `diamond` and a local `.dmnd` database; do not download or build UniRef90 unless explicitly requested.

## Code Style

- Keep changes scoped to the requested behavior.
- Preserve the current config layering and headless behavior.
- Keep optional DIAMOND/UniRef90 fallback disabled by default and clearly marked as lower-confidence annotation unless explicitly configured otherwise.
- Preserve scientific behavior unless the user explicitly asks for scientific changes.
- Avoid introducing global side effects at import time beyond the existing config-loading pattern unless the surrounding code is being intentionally refactored.
- Prefer package-safe imports in new code. Existing script-style imports can be cleaned up in a dedicated import/packaging change.
- Keep biologist-friendly installation and usage documentation current.

## Good First Verification Targets

- `.venv/bin/python -m pytest` when `.venv` exists.
- `.venv/bin/python scripts/check_environment.py` or `python3 scripts/check_environment.py`.
- Config loading with a temporary config directory or `PDB2NET_CONFIG_FILE`.
- File extension and PDB ID parsing.
- Distance detection on small synthetic atom data.
- Headless CX2 export with tiny node/edge fixtures.
- Installed-package check for `pdb2net/configs/config.base.json`.
