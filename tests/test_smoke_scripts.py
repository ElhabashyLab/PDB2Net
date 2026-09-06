"""Run standalone smoke entry points, which pytest collection alone skips."""
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize("script", ["test_batch_exitcode.py", "test_blast_cache.py"])
def test_smoke_script_is_self_contained(script, tmp_path):
    repository = Path(__file__).resolve().parents[1]
    outside_cache = tmp_path / "must-not-use.sqlite3"
    result = subprocess.run(
        [sys.executable, str(repository / "scripts" / script)],
        cwd=repository,
        env={**os.environ, "PDB2NET_BLAST_CACHE_PATH": str(outside_cache),
             "PDB2NET_OUTPUT": str(tmp_path), "PDB2NET_OPEN_IN_CYTOSCAPE": "false"},
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert not outside_cache.exists()
