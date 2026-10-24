"""Amino-acid sequences and point mutations (UniProt numbering throughout)."""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

_SHORT = re.compile(r"^([A-Z])(\d+)([A-Z])$")
_HGVS = re.compile(r"^(?:p\.)?([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


class SequenceError(ValueError):
    pass


def is_residue(code: object) -> bool:
    """True for exactly one canonical amino-acid letter."""
    return isinstance(code, str) and len(code) == 1 and code in AMINO_ACIDS


@dataclass(frozen=True, order=True)
class Mutation:
    position: int
    wt: str
    mt: str

    def __post_init__(self) -> None:
        if not (is_residue(self.wt) and is_residue(self.mt)):
            raise SequenceError(f"non-canonical residue in {self.wt}{self.position}{self.mt}")
        if self.wt == self.mt:
            raise SequenceError(f"synonymous substitution {self.wt}{self.position}{self.mt}")
        if self.position < 1:
            raise SequenceError(f"position must be >= 1, got {self.position}")

    @property
    def id(self) -> str:
        return f"{self.wt}{self.position}{self.mt}"

    @classmethod
    def parse(cls, text: str) -> "Mutation":
        """Accepts ``R856H`` or HGVS protein notation ``p.Arg856His``."""
        text = text.strip()
        if m := _SHORT.match(text.upper()):
            return cls(int(m.group(2)), m.group(1), m.group(3))
        if m := _HGVS.match(text):
            try:
                return cls(int(m.group(2)), THREE_TO_ONE[m.group(1)], THREE_TO_ONE[m.group(3)])
            except KeyError as exc:
                raise SequenceError(f"unknown residue code in {text!r}") from exc
        raise SequenceError(f"cannot parse mutation {text!r}")


def clean_sequence(seq: str, max_length: int | None = None) -> str:
    """Normalise a user/database sequence and reject anything ESMFold can't fold."""
    seq = re.sub(r"\s+", "", seq).upper()
    if seq.startswith(">"):
        raise SequenceError("pass the bare sequence, not a FASTA header")
    if not seq:
        raise SequenceError("empty sequence")
    bad = sorted(set(seq) - set(AMINO_ACIDS))
    if bad:
        raise SequenceError(f"non-canonical residues: {''.join(bad)}")
    if max_length is not None and len(seq) > max_length:
        raise SequenceError(f"sequence has {len(seq)} residues; limit is {max_length}")
    return seq


def apply_mutations(seq: str, mutations: list[Mutation], offset: int = 1) -> str:
    """Apply substitutions to ``seq`` whose first residue has UniProt number ``offset``."""
    chars = list(seq)
    for mut in mutations:
        idx = mut.position - offset
        if not 0 <= idx < len(chars):
            raise SequenceError(f"{mut.id} lies outside residues {offset}-{offset + len(seq) - 1}")
        if chars[idx] != mut.wt:
            raise SequenceError(f"{mut.id}: reference residue at {mut.position} is {chars[idx]}")
        chars[idx] = mut.mt
    return "".join(chars)


def sequence_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()[:16]
