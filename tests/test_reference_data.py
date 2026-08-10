from pathlib import Path

from pdb2net.reference_data import (
    load_pdb_fasta_headers,
    load_sifts_mapping,
    load_sifts_segments,
)


def test_pdb_header_lookup_does_not_retain_reference_sequences(tmp_path: Path) -> None:
    fasta = tmp_path / "pdb_seqres.txt"
    fasta.write_text(
        ">1abc_A mol:protein length:6 Example protein\n"
        "AAAAAA\n"
        ">2xyz_B mol:na length:4 Example DNA\n"
        "ATGC\n"
        ">3def_C mol:protein length:3 Unrequested protein\n"
        "GGG\n",
        encoding="utf-8",
    )
    load_pdb_fasta_headers.cache_clear()

    headers = load_pdb_fasta_headers(str(fasta), ("1ABC", "2xyz"))

    assert headers == {
        "1abc_A": {
            "info": "mol:protein length:6 Example protein",
        },
        "2xyz_B": {
            "info": "mol:na length:4 Example DNA",
        },
    }


def test_real_1ao7_chimera_segments_are_retained_and_compatibility_mapping_is_ambiguous(
    tmp_path: Path,
) -> None:
    """Chain D reflects the two accessions in the official PDBe-SIFTS record."""
    mapping = tmp_path / "uniprot_segments_observed.tsv"
    mapping.write_text(
        "# PDB\tCHAIN\tSP_PRIMARY\tRES_BEG\tRES_END\tPDB_BEG\tPDB_END\tSP_BEG\tSP_END\n"
        "1ao7\tD\tA0A075B6T6\t1\t92\t1\t92\t23\t112\n"
        "1ao7\tD\tP01848\t118\t210\t118\t210\t1\t93\n",
        encoding="utf-8",
    )
    load_sifts_segments.cache_clear()
    load_sifts_mapping.cache_clear()

    segments = load_sifts_segments(str(mapping))
    compatibility = load_sifts_mapping(str(mapping))

    assert [segment["accession"] for segment in segments["1ao7_D"]] == [
        "A0A075B6T6",
        "P01848",
    ]
    assert "1ao7_D" not in compatibility


def test_repeated_external_sifts_segments_for_one_accession_remain_unique(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "uniprot_segments_observed.tsv"
    mapping.write_text(
        "# PDB\tCHAIN\tSP_PRIMARY\tRES_BEG\tRES_END\tPDB_BEG\tPDB_END\tSP_BEG\tSP_END\n"
        "8aly\tE\tQ9H2K2\t4\t71\t4\t71\t867\t934\n"
        "8aly\tE\tQ9H2K2\t99\t270\t99\t270\t962\t1133\n",
        encoding="utf-8",
    )
    load_sifts_segments.cache_clear()
    load_sifts_mapping.cache_clear()

    assert len(load_sifts_segments(str(mapping))["8aly_E"]) == 2
    assert load_sifts_mapping(str(mapping))["8aly_E"] == "Q9H2K2"
