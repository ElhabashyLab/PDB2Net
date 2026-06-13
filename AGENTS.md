# AGENTS.md

## Project

PDB2Net extracts protein interaction networks from PDB/mmCIF structure files and exports Cytoscape-compatible CX2 networks. The pipeline combines Gemmi parsing, SIFTS/FASTA-based annotation, optional BLAST fallback, cKDTree distance searches, and headless or live Cytoscape export.

## Runtime

- Use `python3`, preferably Python 3.11 or 3.12.
- Do not assume a `python` executable exists.
- The current script entry point is:

```bash
python3 -m pdb2net.main
```

- Keep Cytoscape optional for automated work. Prefer `open_in_cytoscape: false` for checks, fixtures, and headless runs.

## Configuration

Configuration is layered in this order:

1. `pdb2net/configs/config.base.json`
2. `pdb2net/configs/config.{linux|windows|darwin}.json`
3. `pdb2net/configs/config.local.json`
4. `PDB2NET_CONFIG_FILE`
5. environment variable overrides

Use `config.local.json` or environment variables for machine-specific paths. Do not commit local reference data paths, generated outputs, BLAST databases, or large datasets.

## Safety

- Do not run full batch jobs on large input folders unless explicitly requested.
- Do not download large reference datasets unless explicitly requested.
- Do not require a running Cytoscape UI for verification.
- Prefer tiny fixtures, headless CX2 export, and focused import/unit checks when validating changes.
- BLAST work requires external `blastp`, `makeblastdb`, a Swiss-Prot FASTA, and a built BLAST database.

## Code Style

- Keep changes scoped to the requested behavior.
- Preserve the current config layering and headless behavior.
- Avoid introducing global side effects at import time beyond the existing config-loading pattern unless the surrounding code is being intentionally refactored.
- Prefer package-safe imports in new code. Existing script-style imports can be cleaned up in a dedicated import/packaging change.

## Good First Verification Targets

- `python3 scripts/check_environment.py`
- Config loading with a temporary config directory or `PDB2NET_CONFIG_FILE`.
- File extension and PDB ID parsing.
- Distance detection on small synthetic atom data.
- Headless CX2 export with tiny node/edge fixtures.
