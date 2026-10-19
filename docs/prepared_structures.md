# Prepared structure analysis

The Python API in `pdb2net.structure_analysis` separates structure-dependent
measurements from reference-dependent annotation. An external adapter owns
storage, bounded file decoding and scheduling. Core owns scientific validation,
annotation, contact selection and the normal network/CSV exporters.

## Preparation

`prepare_structure(source, profile=None, contact_sink=None)` parses one official
PDB structure. It returns facts schema 1 and optionally sends each atom-contact
row to a caller-provided callback. No external annotation references or search
database are required. A callback failure aborts preparation; the caller must
discard incomplete outputs and enforce its streaming byte/row limits.

Embedded annotation extraction distinguishes structural chain-mapping tables
from SIFTS segment tables. Polymer-scheme and atom-site mappings read only their
two chain-ID columns, with a separate ceiling of 20 million rows per mapping
table. Genuine SIFTS segment categories retain their 200,000-row ceiling.
Large polymer tables therefore do not consume the annotation-segment budget.
Scheme precedence, ambiguous mappings and atom-site fallback are unchanged.

The fixed profile covers the asymmetric unit, first model, 5 Å polymer-heavy-atom
contacts and 12 Å C-alpha neighbours, using the current occupancy and alternative
site rules. Its identity includes parser/Gemmi, sequence/type, length and contact
policies, and excludes reference generations, sequence-search policy and rendering.

Facts contain canonical structure/chain identities, observed-polymer search
sequences, structural fallback types, polymer lengths, embedded source evidence,
input diagnostics, the SHA-256 of the expanded source bytes, and candidate-pair
counts. Final UniProt assignments, molecule names and interaction labels are
absent. Every supported chain with contact atoms remains, including isolated
and currently unassigned chains. Candidate pairs have at least one heavy-atom
contact and have not passed annotation-dependent filters.

Contact rows contain, in order:

```text
chain_a, chain_b,
residue_a_number, insertion_a, residue_a_label, atom_a,
residue_b_number, insertion_b, residue_b_label, atom_b,
minimum_distance
```

Residue numbers are strings. Distances are unrounded; one minimum witness is
retained per logical site pair. Rows follow candidate-pair order, then B-site
and A-site lexical order. `validate_contacts` checks strict ordering, endpoint
identity, finite in-range distances and agreement with every candidate count
without accumulating all contacts. Iterators must be consumed completely.

## Analysis

`analyze_prepared(structures, config, output_paths=None, contact_sources=None,
provenance=None, *, web_output_dir=None)` takes facts dictionaries and a complete
configuration snapshot. Run it in a dedicated process with one immutable
reference generation, just as for the ordinary CLI. It validates facts before
materializing mutable annotation records; supplied facts are not modified.

The function applies embedded/external SIFTS and FASTA annotation followed by
eligible BLAST/optional DIAMOND fallback. The shared contact selector then
applies current molecule types and minimum counts. It uses the ordinary CX2,
CSV, output-limit and diagnostic writers. Different model policies or radii
are rejected because the stored measurements cannot answer them exactly.

For detailed output, `contact_sources` maps each canonical PDB ID to a
zero-argument iterator factory. The full contact stream is validated, while
only accepted interactions are exported with their fresh annotations. No
source parsing or coordinate access occurs. Network-only requests need no
contact streams. Source warnings are retained in the generated report.

The storage adapter must bound input file sizes and decompression before JSON
decoding and must verify stored-file integrity. This API defines scientific
payloads, not a directory layout, archive downloader or persistence service.

## Public contracts and provenance

Package `0.2.0rc6` advertises facts schema 1 through capability schema 4. The
CLI contract 3 provides `run` and `capabilities`; persistence is an external
adapter responsibility. The former schema-3 Core store and its commands are
removed, with no migration layer.

Web-output contract 3.0 adds `prepared_inputs`. It is null for raw inputs and
failed reports. Prepared runs include facts schema, geometry profile, source
hashes and revision IDs; an adapter can additionally supply a catalog ID and
checked timestamp. Direct Python analysis without a catalog uses null for
those two snapshot fields. Provenance must match every analyzed input in order
and contains no filesystem paths. `facts_id` and `validate_prepared_provenance`
provide the shared revision-identity and provenance checks. Reference identity
continues to describe the references actually used during annotation.

The four CX2 views, CSV columns and stable output directories are unchanged.
An external adapter must update its capability/output validators together with
this package before admitting real jobs.

`tests/test_structure_analysis.py` checks frozen raw/prepared CX2 and CSV parity,
independent contact expectations, combined structures, fresh BLAST/SIFTS
reference generations, compression-independent source identity, corrupt payload
rejection and atomic detail-export failure. Full-archive capacity and server
publication are separate acceptance scopes.
