"""One polymer/occupancy selection and logical-pair rule for counts and CSV."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Iterator

import numpy as np
from scipy.spatial import cKDTree

from .input_contract import InputValidationError


CONTACT_POLYMER_TYPES = {"PeptideL", "PeptideD", "CyclicPseudoPeptide", "Dna", "Rna", "DnaRnaHybrid"}
SiteKey = tuple[str, str, str]


@dataclass(frozen=True)
class Position:
    site: SiteKey
    residue: str
    atom_name: str
    altloc: str
    original_residue: str
    coordinates: tuple[float, float, float]

    @property
    def order(self):
        return self.residue, self.original_residue, self.atom_name, self.altloc, self.coordinates


def contact_positions(chain: dict[str, Any], extraction_type: str = "all_atoms") -> Iterator[Position]:
    """Keep raw residue data intact; omit only positions outside contact scope."""
    seen: dict[tuple[SiteKey, str, str], tuple[float, float, float]] = {}
    for index, residue in enumerate(chain.get("residues", [])):
        if not residue.get("is_polymer", True):
            continue
        kind = residue.get("polymer_type")
        if kind not in (None, "", "Unknown") and kind not in CONTACT_POLYMER_TYPES:
            continue
        number = residue.get("residue_number", f"internal-{index}")
        insertion = str(residue.get("insertion_code") or "")
        name = str(residue.get("residue_name") or "?")
        original = str(residue.get("original_residue_name") or name)
        for atom in residue.get("atoms", []):
            occupancy = atom.get("occupancy")
            if occupancy is not None:
                if not isinstance(occupancy, (int, float)) or not math.isfinite(occupancy) or occupancy < 0:
                    raise InputValidationError(
                        "INVALID_ATOM_OCCUPANCY", "A polymer contact atom has invalid occupancy."
                    )
                if occupancy == 0:
                    continue
            atom_name = str(atom.get("atom_name") or "")
            if extraction_type == "ca" and atom_name != "CA":
                continue
            if number is None or not atom_name:
                raise InputValidationError("AMBIGUOUS_ATOM_IDENTITY", "A contact atom has no usable residue/atom identity.")
            site = str(number), insertion, atom_name
            altloc = str(atom.get("altloc") or "")
            coordinates = tuple(float(value) for value in atom["coordinates"])
            identity = site, original, altloc
            previous = seen.get(identity)
            if previous is not None:
                if previous != coordinates:
                    raise InputValidationError(
                        "AMBIGUOUS_ATOM_IDENTITY", "Different contact positions share the same atom and alternative-location identity."
                    )
                continue
            seen[identity] = coordinates
            yield Position(site, f"{residue.get('residue_number', '?')}{insertion}:{name}",
                           atom_name, altloc, original, coordinates)


class ContactView:
    def __init__(self, chain: dict[str, Any], extraction_type: str = "all_atoms"):
        self.positions = list(contact_positions(chain, extraction_type))
        self.coordinates = np.array([p.coordinates for p in self.positions], dtype=float).reshape(-1, 3)
        self.tree = cKDTree(self.coordinates) if self.positions else None
        groups: dict[SiteKey, list[int]] = {}
        for index, position in enumerate(self.positions):
            groups.setdefault(position.site, []).append(index)
        self.groups = groups
        self.has_alternatives = len(groups) != len(self.positions)


def count_contacts(left: ContactView, right: ContactView, radius: float) -> int:
    if left.tree is None or right.tree is None:
        return 0
    if not left.has_alternatives and not right.has_alternatives:
        return int(left.tree.count_neighbors(right.tree, radius))
    total = 0
    for indices in right.groups.values():
        neighbours: set[SiteKey] = set()
        for index in indices:
            neighbours.update(left.positions[i].site for i in left.tree.query_ball_point(right.coordinates[index], radius))
        total += len(neighbours)
    return total


def unique_contacts(left: ContactView, right: ContactView, radius: float, *, budget=None):
    """Yield one minimum-distance witness per site pair, bounded by one right site."""
    if left.tree is None or right.tree is None:
        return
    maximum_left_positions = max(map(len, left.groups.values()))
    for site_b in sorted(right.groups):
        minima: dict[SiteKey, tuple[float, int, int]] = {}
        indices_b = right.groups[site_b]
        for index_b in indices_b:
            if budget is not None:
                size = int(left.tree.query_ball_point(right.coordinates[index_b], radius, return_length=True))
                # A site has at most this many positions. This lower bound is
                # exact without left alternatives and never charges duplicate
                # positions as CSV rows. Check it before allocating the list.
                budget.ensure_rows_available((size + maximum_left_positions - 1) // maximum_left_positions)
            for index_a in left.tree.query_ball_point(right.coordinates[index_b], radius):
                atom_a, atom_b = left.positions[index_a], right.positions[index_b]
                distance = float(np.linalg.norm(left.coordinates[index_a] - right.coordinates[index_b]))
                previous = minima.get(atom_a.site)
                if previous is None:
                    if budget is not None:
                        budget.ensure_rows_available(len(minima) + 1)
                    minima[atom_a.site] = distance, index_a, index_b
                else:
                    old_distance, old_a, old_b = previous
                    if (distance, atom_a.order, atom_b.order) < (
                        old_distance, left.positions[old_a].order, right.positions[old_b].order
                    ):
                        minima[atom_a.site] = distance, index_a, index_b
        for site_a in sorted(minima):
            distance, index_a, index_b = minima[site_a]
            yield left.positions[index_a], right.positions[index_b], distance
