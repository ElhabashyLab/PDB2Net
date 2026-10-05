# 8ALY NextGen fixture provenance

`pdb_00008aly_xyz-enrich.cif` is a hand-reduced, static test fixture derived
from the official wwPDB NextGen example for entry 8ALY. It was not emitted by
PDB2Net production code. The reduction keeps the original data-block and entry
identity forms, quoted and semicolon-delimited mmCIF values, the genuine
`label_asym_id` A to author-chain E mapping, the UniProt segment, both Pfam
segments, and representative atom-site rows.

- Source: https://files-nextgen.wwpdb.org/pdb_nextgen/data/entries/divided/al/pdb_00008aly/pdb_00008aly_xyz-enrich.cif.gz
- Retrieved: 2026-08-10
- Original compressed size: 1,294,445 bytes
- Original compressed SHA-256: `ce1667a7fe5e7bb8fed6656edf7a49bd3cebf86f15c01187dca7d55223197453`
- Archive contract: https://www.wwpdb.org/ftp/pdb-nextgen-archive-site

The checked-in reduced file is 3,005 bytes as LF-normalized UTF-8 text, SHA-256
`c53dbbbc5629cc0899afcc16ae5770e0632e400bad625ebf63be930d51030661`.
The test normalizes text newlines when checking this digest so Windows Git
checkouts with CRLF line endings describe the same fixture.
Its two coordinate residues are ASP 875 and PHE 876 in author chain E
(label chain A). Expected embedded facts are UniProt Q9H2K2, PDB sequence
positions 4–299 mapped to UniProt 867–1162, and Pfam PF00644 at 99–270 plus
PF07647 at 8–71. This excerpt does not represent the full structure geometry.

`tests/test_nextgen_cli.py` checks those facts through the real headless CLI,
offline precompute and read-only assembly in the default test suite. Its tiny
FASTA/TSV lookups are synthetic test support, not biological reference snapshots.
