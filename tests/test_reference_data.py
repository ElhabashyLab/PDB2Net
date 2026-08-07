from pathlib import Path

from pdb2net.reference_data import load_pdb_fasta_headers


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
