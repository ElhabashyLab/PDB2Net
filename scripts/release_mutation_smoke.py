#!/usr/bin/env python3
"""Prove that the contract 1.2 regression suite kills critical mutations."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


@dataclass(frozen=True)
class Mutation:
    name: str
    scope: str
    relative_file: str
    original: str
    replacement: str
    test: str


MUTATIONS = (
    Mutation(
        "extended ID is no longer normalized to lowercase",
        "core",
        "pdb2net/structure_identity.py",
        "        return raw.lower()\n",
        "        return raw\n",
        "tests/test_contract_1_2_regressions.py::test_official_id_grammar_normalization_and_archive_shard_are_exact",
    ),
    Mutation(
        "ambiguous external SIFTS mapping is assigned arbitrarily",
        "core",
        "pdb2net/reference_data.py",
        "        if len(accessions) == 1:\n            pdb_to_uniprot[key] = next(iter(accessions))\n",
        "        if accessions:\n            pdb_to_uniprot[key] = next(iter(accessions))\n",
        "tests/test_reference_data.py::test_real_1ao7_chimera_segments_are_retained_and_compatibility_mapping_is_ambiguous",
    ),
    Mutation(
        "truncated gzip EOF is accepted",
        "core",
        "pdb2net/file_parser.py",
        (
            "                        if not chunk:\n"
            "                            raw_eof = True\n"
            "                            raise InputValidationError(\n"
            "                                \"INVALID_GZIP_INPUT\",\n"
            "                                f\"Gzip input is truncated or has no valid trailer: {Path(file_path).name}\",\n"
            "                            )\n"
        ),
        "                        if not chunk:\n                            raw_eof = True\n                            break\n",
        "tests/test_contract_1_2_regressions.py::test_every_truncated_gzip_trailer_is_rejected",
    ),
    Mutation(
        "annotation pipeline version is omitted from its cache profile",
        "core",
        "pdb2net/precomputed_store.py",
        '        "annotation_pipeline_version": ANNOTATION_PIPELINE_VERSION,\n',
        "",
        "tests/test_precomputed_store.py::test_annotation_pipeline_version_invalidates_only_annotations",
    ),
    Mutation(
        "tooltip selection invalidates scientific annotation cache entries",
        "core",
        "pdb2net/precomputed_store.py",
        '        "use_embedded_sifts": annotation_cfg["use_embedded_sifts"],\n',
        (
            '        "use_embedded_sifts": annotation_cfg["use_embedded_sifts"],\n'
            '        "tooltip_fields": annotation_cfg["tooltip_fields"],\n'
        ),
        "tests/test_extended_ids_and_nextgen.py::test_tooltip_selection_does_not_change_precompute_profiles",
    ),
    Mutation(
        "a CX1 nodeAttributes aspect is reintroduced",
        "core",
        "pdb2net/cx2_export.py",
        "        visual_props,\n        {\"status\": [{\"error\": \"\", \"success\": True}]},\n",
        (
            "        {\"nodeAttributes\": []},\n"
            "        visual_props,\n"
            "        {\"status\": [{\"error\": \"\", \"success\": True}]},\n"
        ),
        "tests/test_export_contracts.py::test_headless_cx2_uses_only_native_inline_attributes_and_layout",
    ),
    Mutation(
        "the required artifacts summary field is omitted from the Worker contract",
        "web",
        "worker/output_collector.py",
        '    "artifacts",\n',
        "",
        "tests/test_output_collector.py::test_real_contract_rejects_every_missing_summary_field",
    ),
)


def _copy_scope(destination: Path, source: Path, scope: str) -> None:
    if scope == "core":
        shutil.copytree(source / "pdb2net", destination / "pdb2net")
        shutil.copytree(source / "tests", destination / "tests")
        return
    shutil.copytree(source / "worker", destination / "worker")
    shutil.copytree(source / "tests", destination / "tests")


def _run_mutation(
    mutation: Mutation,
    *,
    core_root: Path,
    web_root: Path,
    python: Path,
) -> None:
    source = core_root if mutation.scope == "core" else web_root
    with tempfile.TemporaryDirectory(prefix="pdb2net-mutant-") as temporary:
        destination = Path(temporary)
        _copy_scope(destination, source, mutation.scope)
        target = destination / mutation.relative_file
        content = target.read_text(encoding="utf-8")
        occurrences = content.count(mutation.original)
        if occurrences != 1:
            raise RuntimeError(
                f"Mutation anchor for {mutation.name!r} occurred {occurrences} times."
            )
        target.write_text(
            content.replace(mutation.original, mutation.replacement, 1),
            encoding="utf-8",
        )

        environment = os.environ.copy()
        environment["PYTHONPATH"] = str(
            destination if mutation.scope == "core" else destination / "worker"
        )
        result = subprocess.run(
            [str(python), "-m", "pytest", "-q", mutation.test],
            cwd=destination,
            env=environment,
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        output = result.stdout + result.stderr
        if result.returncode != 1 or re.search(r"\b[1-9][0-9]* failed\b", output) is None:
            raise RuntimeError(
                f"Mutation survived or did not produce a clean assertion failure: {mutation.name}\n"
                f"pytest exit={result.returncode}\n{output[-4_000:]}"
            )
        print(f"killed: {mutation.name}")


def main() -> int:
    core_default = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core-root", type=Path, default=core_default)
    parser.add_argument(
        "--web-root",
        type=Path,
        default=core_default.parent / "pdb2net-webserver",
    )
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    args = parser.parse_args()

    core_root = args.core_root.expanduser().resolve()
    web_root = args.web_root.expanduser().resolve()
    # Preserve a virtualenv launcher symlink instead of resolving it to the
    # system interpreter, which may not carry the test dependencies.
    python = args.python.expanduser().absolute()
    if not (core_root / "pdb2net").is_dir():
        parser.error(f"Core source tree not found: {core_root}")
    if not (web_root / "worker").is_dir():
        parser.error(f"Webserver source tree not found: {web_root}")
    if not python.is_file():
        parser.error(f"Python interpreter not found: {python}")

    for mutation in MUTATIONS:
        _run_mutation(
            mutation,
            core_root=core_root,
            web_root=web_root,
            python=python,
        )
    print(f"Mutation smoke passed: {len(MUTATIONS)} critical mutants killed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
