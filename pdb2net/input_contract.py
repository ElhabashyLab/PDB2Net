"""Shared input-contract errors for folder-based processing."""

from __future__ import annotations


class InputValidationError(ValueError):
    """Raised when the configured folder-based input contract is violated."""

    def __init__(self, code: str, message: str) -> None:
        self.code = str(code)
        self.message = str(message)
        super().__init__(f"{self.code}: {self.message}")

    def __reduce__(self):
        """Keep the structured code when an error crosses a process boundary."""
        return self.__class__, (self.code, self.message)


# Only content problems attributable to one input may be omitted from a raw
# run. Resource, integrity, tool and unexpected implementation failures remain
# fatal. Keep this separate from the generic parser's optional error collector.
INPUT_EXCLUSION_REASONS = {
    "EMPTY_INPUT_FILE": "The file is empty. Replace it with a structure file or remove it.",
    "INVALID_STRUCTURE_TEXT_ENCODING": "The file is not valid UTF-8 text. Export it again as a PDB/mmCIF text file.",
    "INVALID_PDB_INPUT": "The PDB structure could not be read. Check the file or export it again.",
    "INVALID_MMCIF_INPUT": "The mmCIF structure could not be read. Check the file or export it again.",
    "INVALID_MMCIF_DATA_BLOCK_COUNT": "Exactly one mmCIF data block is required. Submit each structure as a separate file.",
    "INVALID_STRUCTURE_IDENTITY": "The structure contains an invalid identity. Correct the identity or export the file again.",
    "CONFLICTING_STRUCTURE_IDENTITY": "The structure contains conflicting identities. Correct the file's identity fields.",
    "MISSING_CHAIN_ID": "A chain in a selected model has no identifier. The file was omitted. Correct the chain identifiers or export the structure again.",
    "NO_VALID_CHAINS": "No supported chains were found. Check the structure and model selection or remove this file.",
    "INVALID_ATOM_OCCUPANCY": "A polymer contact atom has invalid occupancy. Correct nonnumeric, nonfinite or negative occupancy values and export the file again.",
    "AMBIGUOUS_ATOM_IDENTITY": "Atom positions or metadata cannot be assigned unambiguously. Check residue numbers, insertion codes, atom names and alternative-location labels.",
}


def is_excludable_input_error(error: Exception) -> bool:
    return isinstance(error, InputValidationError) and error.code in INPUT_EXCLUSION_REASONS
