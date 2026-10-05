# Frozen Core reference dataset

This dataset contains 23 full structure files and three small, authentic
reference subsets. It is independent of production reference updates. Tests
do not download or refresh these files. Generated search databases, caches,
stores and outputs belong in pytest's temporary directory, never here.

Run the full suite with `.venv/bin/python -m pytest -q -ra`, or only these real
cases with `.venv/bin/python -m pytest tests/integration -q -ra`. Local BLAST+
is required for the real runs. `-m "not integration"` selects the smaller suite.
See `docs/development.md` for external copies, logs and candidate regeneration.

## Structure purposes

| Files | Purpose |
| --- | --- |
| `4hhb.pdb` | Repeated protein identities, interchain contacts, excluded haem. |
| `9rat.pdb` | Small single-chain protein; controlled search miss. |
| `1TUP.pdb`, `1tup.cif` | Protein/DNA contacts and paired-format equivalence. |
| `1a34.cif` | Protein/RNA contacts, two RNA chains and alternative positions. |
| `124d.cif` | DNA/RNA chains and hydrogen exclusion. |
| `170d.cif` | Hybrid DNA/RNA chains with CAR; F06 length regression. |
| `5rlu.pdb` | Larger protein chains, alternatives and excluded ligands. |
| `6w41.cif` | Insertion codes, antibodies and spike protein. |
| `6m17.cif` | Larger protein complex, excluded glycans and other nonpolymers. |
| `1a1t.pdb`, `1a1t.cif` | Protein/RNA, 25 models, first/all and model separation. |
| `1a8o.pdb`, `1a8o.cif` | MSE, water, zero occupancy and BLAST fallback. |
| `8aly.cif`, `pdb_00008aly_xyz-enrich.cif.gz` | Full regular/NextGen geometry, embedded annotation, gzip. |
| `AF-P69905-F1-model_v6.pdb`, `AF-P69905-F1-model_v6.cif` | AlphaFold identity and renamed-file annotation. |
| `123456.cif` | Deliberately edited 9JR2 identity conflict; expected rejection. |
| `1bna.cif` | Pure DNA duplex; exhaustive real contact/CSV check. |
| `1ehz.cif` | Pure RNA, modified nucleotides, no interchain edges; F06. |
| `2erk.cif` | Polymer phosphothreonine and phosphotyrosine. |
| `1ao7.cif` | Multiple proteins; two SIFTS accessions for chain D. |

Run format variants separately: they share structure identities. Do not feed
the whole directory into one raw-file job. The original small NextGen excerpt
and mini-UniRef90 fixture remain in their existing sibling directories.

## References and provenance

`provenance.json` records exact bytes, hashes, source URLs/dates where known,
reference selection and full-source hashes. Original download dates for the
older local Datensatz22 files are unverified; this is not concealed by assigning
a new download date. The eight September 2026 downloads and four new official
downloads have recorded provenance. The older, reduced NextGen excerpt has its
own August provenance; it is not the September full-file snapshot.

All three subsets come from the same available local reference corpus. The
SIFTS header identifies 27 July 2025 / UniProt 2025.03. This does not establish
the release dates of the two FASTA files: their exact source hashes identify
the snapshot. The files are frozen together, not claimed to share a release date.

All SIFTS rows for the selected PDB IDs are retained, including ambiguity and
segments. FASTA selection retains complete original records. Swiss-Prot contains
24 entries, including mouse alpha globin P01942 as a competitor for P69905.
Q75677 is present in SIFTS but absent from this Swiss-Prot source; that absence
is intentional, not filled with an invented sequence. The identity-conflict
fixture retains the 9JR2 reference records as part of the original selection.

Tests build a small BLAST database from this FASTA once per session. Tool
versions and reference hashes are recorded with the run. Scores from this small
database are not asserted to reproduce a full-production-database search.

## Scientific expectations

`config.json` freezes the scientific settings: 12 Å Cα radius, 5 Å heavy-atom
radius, ten Cα neighbours for protein/protein, other contact minima one, first
model by default, headless output and DIAMOND disabled. Explicit test variants
exercise all-model chain output and controlled fallback/reference omissions.
Machine-specific paths and environment overrides are isolated by test fixtures.

`expected/structures.json` records independently extracted structure facts and
reviewed annotation expectations. `expected/networks.json` records geometric
counts; `expected/search.json` records controlled search outcomes. Distances are
derived by the independent arithmetic in `tests/reference_oracle.py`, without
Core's atom selectors or KD-tree implementation. No Core export is treated as
an independent source of correct values. Expectations omit transient paths,
timestamps and layout coordinates.

The independent checks exposed F06 (1EHZ incorrectly reported 62 instead of
76 nt; 170D reported 11 instead of 12 nt per chain) and F07 (6M17 A/C counted
nonpolymer LEU 705 and reported 606 instead of 605 aa). Both are now corrected
by counting existing polymer roles/types. All length assertions run normally;
the independently derived fixture expectations were not changed. Raw and
read-only store exports are also checked against these correct values.

See the [validation summary](../../../docs/core-quality-audit.md) for the
current verification scope and [Development](../../../docs/development.md)
for reproducible checks.
