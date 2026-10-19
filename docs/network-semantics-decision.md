# Network projections

Network outputs are projections of already annotated chains and filtered
[contacts](contact-semantics-decision.md). Their weights count atom-site pairs;
they are not affinity, normalized interface size or independent evidence counts.

## Supplied coordinates and model policy

Raw runs analyze the supplied coordinates. Recognizing biological-assembly
operators does not create additional coordinate copies, select an assembly or
remove crystal contacts. A prepared assembly can be used when its intended
chains have usable, unique identifiers. The fixed offline store uses asymmetric
units; prepared-input analysis builds networks from those coordinates without
generating biological assemblies.

`first` is the default model policy. `all` calculates within-model contacts only,
qualifies chain identities by model and produces separate per-model chain
networks. Both protein-output options are rejected before analysis under `all`.
Reported model numbers are enumeration positions, not necessarily source MODEL
serial numbers.

## Four output views

| Output | Nodes and edges |
| --- | --- |
| Chain per structure | Supported retained chains, including unassigned and zero-contact chains; measured interchain contacts |
| Protein per structure | Assigned UniProt proteins plus their contacting nucleic-acid chains; summed contacts between distinct protein IDs or between a protein and a nucleic-acid chain |
| Combined chain | Components linked across input structures by shared UniProt identities, together with measured contacts and their connected chain neighbors |
| Combined protein | UniProt proteins from those linked chain components; summed protein–protein contacts between distinct IDs |

Protein aggregation merges chains sharing a UniProt accession. Contacts between
different chains with the same accession do not become self-loops. Unassigned
protein chains and their contacts remain available in chain outputs but do not
receive artificial UniProt nodes. A missing UniProt ID does not prevent contact
calculation when the molecule type is known; an unknown molecule type is not an
accepted interaction class.

Per-structure protein outputs include only nucleic-acid chains directly contacting
an assigned protein. Combined protein outputs contain no nucleic-acid nodes.
Sums may increase when additional similar structures are supplied and cannot
recover individual chain contributions without the chain networks.

## Identity links

Identity links connect chains sharing a UniProt accession across distinct input
structures. For each accession present in at least two structures, the exporter
uses a deterministic cross-structure spanning tree: all participating chains are
connected, no same-structure identity edge is added, and the group has `n - 1`
identity edges. Combined outputs are these linked components, not a union of
every input structure.

Identity edges have `interaction_type: "identity"` and contact count zero.
They do not represent measured distances between files. Ordinary contact edges
remain positive. The Python layout uses a separate historical attraction strength
for identity edges, so their displayed contact count does not control layout
strength. Live Cytoscape may calculate a different layout.

Store entries contain within-structure contacts. Identity links are generated
during export. All entries must satisfy the current complete store profile,
even though changing identity-edge presentation alone does not change stored
geometry.

## Verification

Tests check all four projections with hand-counted weights and fixed real
structure expectations. Identity checks verify permissible partners, connectivity,
zero weights and spanning-tree size independently of the production builder.
Examples are available in `docs/examples/q02_network_examples.py`; the maintained
test workflow and fixture sources are in [Development](development.md).
