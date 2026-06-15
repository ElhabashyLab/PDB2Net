"""Command-line interface for PDB2Net."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Any


def _add_bool_pair(parser: argparse.ArgumentParser, name: str, dest: str, help_text: str) -> None:
    group = parser.add_mutually_exclusive_group()
    group.add_argument(f"--{name}", dest=dest, action="store_true", help=help_text)
    group.add_argument(f"--no-{name}", dest=dest, action="store_false", help=f"Disable {help_text.lower()}")
    parser.set_defaults(**{dest: None})


def build_parser() -> argparse.ArgumentParser:
    """Create the top-level CLI parser."""
    parser = argparse.ArgumentParser(prog="pdb2net", description="Extract interaction networks from PDB/mmCIF files.")
    subparsers = parser.add_subparsers(dest="command")

    run = subparsers.add_parser("run", help="Run PDB2Net on a folder of structure files.")
    run.add_argument("--input-dir", required=True, help="Folder containing .pdb, .cif, or .mmcif files.")
    run.add_argument("--output-dir", required=True, help="Folder where timestamped internal run output is written.")
    run.add_argument("--web-output-dir", help="Optional stable output folder with summary.json, networks/, and interactions/.")
    run.add_argument("--config", help="Optional JSON config file loaded after local config.")
    run.add_argument("--headless", action="store_true", help="Disable live Cytoscape and write CX2 files only.")
    run.add_argument("--open-in-cytoscape", choices=["true", "false"], help="Explicitly enable or disable Cytoscape UI export.")
    run.add_argument("--pdb-fasta", help="Path to pdb_seqres.txt.")
    run.add_argument("--uniprot-fasta", help="Path to UniProt Swiss-Prot FASTA.")
    run.add_argument("--sifts-tsv", help="Path to SIFTS pdb_chain_uniprot.tsv.")
    run.add_argument("--blast-db", help="Folder containing the UniProt BLAST database prefix uniprot_db.")
    run.add_argument("--blastp", help="blastp executable name or path.")
    run.add_argument("--ca-radius", type=float, help="C-alpha distance radius in Angstrom.")
    run.add_argument("--all-atoms-radius", type=float, help="All-atoms distance radius in Angstrom.")
    _add_bool_pair(run, "chain-per-pdb", "chain_per_pdb", "Export per-PDB chain networks.")
    _add_bool_pair(run, "protein-per-pdb", "protein_per_pdb", "Export per-PDB protein networks.")
    _add_bool_pair(run, "combined-chain-network", "combined_chain_network", "Export combined chain networks.")
    _add_bool_pair(run, "combined-protein-network", "combined_protein_network", "Export combined protein networks.")
    _add_bool_pair(
        run,
        "export-detailed-interactions",
        "export_detailed_interactions",
        "Export detailed atom/residue interaction CSV files.",
    )

    return parser


def _apply_run_overrides(args: argparse.Namespace, config: dict[str, Any]) -> None:
    """Apply CLI values to the loaded config in place."""
    config["input_folder_path"] = str(Path(args.input_dir).expanduser())
    config["output_path"] = str(Path(args.output_dir).expanduser())

    path_options = {
        "pdb_fasta_path": args.pdb_fasta,
        "uniprot_fasta_path": args.uniprot_fasta,
        "sifts_tsv_path": args.sifts_tsv,
        "blast_db_path": args.blast_db,
        "blastp_executable": args.blastp,
    }
    for key, value in path_options.items():
        if value:
            config[key] = str(Path(value).expanduser()) if key != "blastp_executable" else value

    if args.headless:
        config["open_in_cytoscape"] = False
    elif args.open_in_cytoscape is not None:
        config["open_in_cytoscape"] = args.open_in_cytoscape == "true"

    config.setdefault("distance_thresholds", {})
    if args.ca_radius is not None:
        config["distance_thresholds"]["ca_radius"] = args.ca_radius
    if args.all_atoms_radius is not None:
        config["distance_thresholds"]["all_atoms_radius"] = args.all_atoms_radius

    config.setdefault("networks", {})
    for key in [
        "chain_per_pdb",
        "protein_per_pdb",
        "combined_chain_network",
        "combined_protein_network",
    ]:
        value = getattr(args, key)
        if value is not None:
            config["networks"][key] = value

    if args.export_detailed_interactions is not None:
        config["export_detailed_interactions"] = args.export_detailed_interactions


def run_command(args: argparse.Namespace) -> int:
    """Run the pipeline through the config-backed public API."""
    if args.config:
        os.environ["PDB2NET_CONFIG_FILE"] = args.config

    from .config_loader import config

    _apply_run_overrides(args, config)

    if config.get("open_in_cytoscape", True):
        from .cytoscape import ensure_cytoscape_running

        ensure_cytoscape_running()

    try:
        from .pipeline import run_pipeline

        output_paths = run_pipeline(config["input_folder_path"], web_output_dir=args.web_output_dir)
    except Exception as exc:
        print(f"PDB2Net run failed: {exc}", file=sys.stderr)
        if args.web_output_dir:
            print(f"Web output summary: {Path(args.web_output_dir) / 'summary.json'}", file=sys.stderr)
        return 1

    print(f"PDB2Net run complete: {output_paths.run_output_path}")
    print(f"Run summary: {output_paths.summary_file}")
    if args.web_output_dir:
        print(f"Web output summary: {Path(args.web_output_dir) / 'summary.json'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "run":
        return run_command(args)

    parser.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
