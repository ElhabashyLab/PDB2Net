# AGENTS.md

## Project

PDB2Net extracts protein interaction networks from PDB/mmCIF structure files and exports Cytoscape-compatible CX2 networks. The pipeline combines Gemmi parsing, SIFTS/FASTA-based annotation, optional BLAST fallback, cKDTree distance searches, and headless or live Cytoscape export.

## Runtime

- Use `python3`, preferably Python 3.11 or 3.12.
- Do not assume a `python` executable exists.
- The legacy script entry point is:

```bash
python3 -m pdb2net.main
```

- Prefer the stable CLI for new automation:

```bash
python3 -m pdb2net run --input-dir /path/to/inputs --output-dir /path/to/outputs --headless
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

## Output Contract

- Internal runs use timestamped folders with `combined/`, `protein/`, `chain/`, `distances/`, `runtime_analysis.txt`, `manifest.json`, and `run_summary.json`.
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

## Code Style

- Keep changes scoped to the requested behavior.
- Preserve the current config layering and headless behavior.
- Preserve scientific behavior unless the user explicitly asks for scientific changes.
- Avoid introducing global side effects at import time beyond the existing config-loading pattern unless the surrounding code is being intentionally refactored.
- Prefer package-safe imports in new code. Existing script-style imports can be cleaned up in a dedicated import/packaging change.
- Keep biologist-friendly installation and usage documentation current.

## Good First Verification Targets

- `python3 scripts/check_environment.py`
- Config loading with a temporary config directory or `PDB2NET_CONFIG_FILE`.
- File extension and PDB ID parsing.
- Distance detection on small synthetic atom data.
- Headless CX2 export with tiny node/edge fixtures.
