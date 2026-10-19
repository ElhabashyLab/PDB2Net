# Contact rules

PDB2Net reports direct geometric contacts between supported polymer chains in
the supplied coordinates. These rules define atom eligibility and counting;
[network projections](network-semantics-decision.md) define subsequent aggregation.

## Eligible atoms and chains

Protein, DNA and RNA polymer residues contribute their heavy atoms, including
modification atoms. Water, separate glycans, lipids, cofactors, ions and other
nonpolymer ligands are excluded. A blanket exclusion of HETATM records is not
appropriate: modified polymer residues can use those records.

Gemmi polymer roles and types supply the structural metadata. Hydrogen and
deuterium are removed. Explicit zero occupancy excludes a position; positive
occupancy and unknown occupancy retain it without weighting. Missing occupancy
does not produce a warning. Invalid occupancy of a relevant contact atom is a
file-content error.

A chain with no eligible contact atoms is omitted with a chain diagnostic. If
no usable chain remains, the raw-run file-exclusion rules apply. A retained chain
may legitimately have no contacts. Empty chain identifiers in selected models
are file errors; PDB2Net does not invent identifiers.

## Alternative positions and pair identity

All eligible alternative positions participate. A logical atom site is identified
within its structure, model and author chain by author residue number, insertion
code and atom name. Different altlocs or residue names at that site do not
multiply the contact count.

An atom-site pair contributes once if any eligible position pair is within the
cutoff. This is a union of possible contacts, not an assertion that all
alternative conformations occur simultaneously. Contacts are not restricted to
equal altloc letters across chains.

Distance cutoffs are inclusive and applied before rounding. Detailed CSV emits
one minimum-distance witness per logical pair. Ties use residue label, original
residue name, atom name, altloc and coordinates in deterministic order. Only the
CSV distance is rounded to two decimal places. Insertion codes remain part of
the residue label, for example `10A:ALA`.

At the numerical cutoff boundary, inclusion uses an exact squared-distance
comparison of the shortest round-trip decimal representations of the parsed
coordinates and radius. This prevents ordinary decimal boundary contacts from
disappearing through floating-point arithmetic. A coordinate-scale guard of eight
floating-point ULPs only selects candidates for that comparison; it does not
increase the scientific cutoff. Pairs just outside remain excluded. Original CIF
number tokens with precision beyond the parsed floats are not reconstructed.
The same decision applies to counts and detailed witnesses at both radii.

## Interaction filters and models

With the default profile, protein pairs must have at least ten C-alpha pairs
within 12 Å and at least one heavy-atom pair within 5 Å. Protein–nucleic-acid and
nucleic-acid pairs require at least one heavy-atom pair within 5 Å. The raw CLI
supports configured thresholds and filters; the fixed prepared-input profile is described
in [Prepared structures](prepared_structures.md).

Only interchain contacts within the same model are calculated. `first` selects
the first model; `all` keeps model-qualified chain identities without calculating
contacts between models. Neither policy generates biological assemblies or
applies AlphaFold confidence filtering.

Annotation sequences and displayed polymer lengths use observed residue records,
not the occupancy-filtered contact view. Missing residues are not reconstructed.
Protein sequences preserve canonical selenocysteine (`SEC` → `U`) and
pyrrolysine (`PYL` → `O`). Search preparation still masks `O` as `X` under its
existing policy; the stored structural sequence retains `O`.
See the [README](../README.md) and [architecture](architecture.md) for the
distinction between parsing, annotation, contact selection and compaction.

## Verification and compatibility

Small contact oracles independently check cutoffs, rigid transformations,
alternative-site unions, minimum-distance witnesses and insertion codes. The
frozen real-structure suite supplies independent expected contact counts; the
1BNA CSV is compared against every expected atom pair. See
[Development](development.md) and the [fixture provenance](../tests/fixtures/core_reference/README.md).

Parser and interaction fingerprints distinguish incompatible geometry rules.
Current version identifiers and rebuild requirements are maintained in
[Prepared structures](prepared_structures.md). A shared wrong result in raw and
prepared output is not, by itself, independent scientific validation.
