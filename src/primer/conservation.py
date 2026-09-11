"""Conservation analysis and degenerate primer construction.

The premise of the whole project.

A diagnostic PCR primer is a short sequence that must match the pathogen it is looking
for. Viruses mutate. When a mutation lands under a primer binding site, the primer
stops binding and the test returns **negative** — not "inconclusive", not an error. The
assay goes blind while continuing to report confidently.

Cotton Leaf Curl Virus is a live example: resistance-breaking strains circulating in
Punjab have defeated both resistant cultivars and the PCR primers used to detect them.

Two defences, and this module implements both:

  **Target conserved regions.** Positions that have not varied across every observed
  strain are the least likely to vary next season.

  **Degenerate bases.** Where variation exists, an IUPAC code such as `R` (A or G)
  matches both. Degeneracy costs specificity — every degenerate position multiplies the
  number of distinct molecules in the tube — so it is spent deliberately, not sprinkled.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass

# IUPAC ambiguity codes: symbol -> the set of bases it matches.
IUPAC: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),  # puRine
    "Y": frozenset("CT"),  # pYrimidine
    "S": frozenset("GC"),  # Strong (3 H-bonds)
    "W": frozenset("AT"),  # Weak (2 H-bonds)
    "K": frozenset("GT"),  # Keto
    "M": frozenset("AC"),  # aMino
    "B": frozenset("CGT"),  # not A
    "D": frozenset("AGT"),  # not C
    "H": frozenset("ACT"),  # not G
    "V": frozenset("ACG"),  # not T
    "N": frozenset("ACGT"),
}

# Reverse lookup: the smallest code covering a set of bases.
_BY_SET = {bases: code for code, bases in IUPAC.items()}


class AlignmentError(ValueError):
    pass


def code_for(bases: Sequence[str]) -> str:
    """The IUPAC symbol covering exactly these bases."""
    observed = frozenset(b.upper() for b in bases if b.upper() in "ACGT")
    if not observed:
        return "N"
    return _BY_SET.get(observed, "N")


def matches(code: str, base: str) -> bool:
    return base.upper() in IUPAC.get(code.upper(), frozenset())


def degeneracy(sequence: str) -> int:
    """Number of distinct molecules a degenerate primer represents.

    The cost of degeneracy, made explicit. A primer with six `N` positions is 4096
    different oligos sharing one tube, so each is at 1/4096 of the nominal
    concentration — which is why an over-degenerate primer simply fails.
    """
    total = 1
    for symbol in sequence.upper():
        total *= len(IUPAC.get(symbol, frozenset("ACGT")))
    return total


@dataclass(frozen=True)
class Position:
    index: int
    counts: dict[str, int]
    consensus: str
    conservation: float  # 1.0 = invariant
    entropy: float  # 0.0 = invariant, 2.0 = uniform over 4 bases
    gaps: int


def check_alignment(sequences: Sequence[str]) -> None:
    if not sequences:
        raise AlignmentError("no sequences")
    lengths = {len(s) for s in sequences}
    if len(lengths) != 1:
        # Unequal lengths mean the sequences were never aligned, and column-wise
        # analysis of unaligned sequences is meaningless rather than merely inaccurate.
        raise AlignmentError(f"sequences are not aligned: lengths {sorted(lengths)}")


def analyse(alignment: Sequence[str]) -> list[Position]:
    """Per-column conservation across an aligned set of strains."""
    check_alignment(alignment)
    positions: list[Position] = []

    for i in range(len(alignment[0])):
        column = [s[i].upper() for s in alignment]
        gaps = sum(1 for b in column if b in "-.")
        bases = [b for b in column if b in "ACGT"]
        counts = Counter(bases)

        if not bases:
            positions.append(Position(i, {}, "N", 0.0, 2.0, gaps))
            continue

        consensus, top = counts.most_common(1)[0]
        # Conservation is the majority share among *called* bases. Gaps are reported
        # separately rather than counted as disagreement, because a gap is missing
        # data, not an observed difference.
        conservation = top / len(bases)

        entropy = -sum((n / len(bases)) * math.log2(n / len(bases)) for n in counts.values())

        positions.append(
            Position(
                index=i,
                counts=dict(counts),
                consensus=consensus,
                conservation=round(conservation, 6),
                entropy=round(entropy, 6),
                gaps=gaps,
            )
        )

    return positions


def conserved_windows(
    positions: Sequence[Position],
    *,
    length: int = 20,
    min_conservation: float = 0.99,
    max_gap_fraction: float = 0.0,
) -> list[tuple[int, float]]:
    """Windows where every position clears the conservation threshold.

    Returns (start, mean conservation). Requiring *every* position to clear the bar,
    rather than the window mean, is deliberate: one highly variable position under a
    primer is enough to break it, and a mean happily hides it behind nineteen
    invariant neighbours.
    """
    if length <= 0:
        raise ValueError("length must be positive")

    n_seqs_gap_limit = max_gap_fraction
    out: list[tuple[int, float]] = []

    for start in range(len(positions) - length + 1):
        window = positions[start : start + length]
        if any(p.conservation < min_conservation for p in window):
            continue
        total_bases = sum(sum(p.counts.values()) + p.gaps for p in window) or 1
        if sum(p.gaps for p in window) / total_bases > n_seqs_gap_limit:
            continue
        out.append((start, round(sum(p.conservation for p in window) / length, 6)))

    return out


def degenerate_consensus(
    alignment: Sequence[str], start: int, length: int, *, min_frequency: float = 0.05
) -> str:
    """Build a degenerate primer covering the variation in a window.

    `min_frequency` excludes singletons: a base seen once in four hundred strains is
    usually a sequencing error, and covering it costs real specificity for no gain.
    That threshold is the entire difference between a usable degenerate primer and a
    useless one.
    """
    check_alignment(alignment)
    out: list[str] = []

    for i in range(start, start + length):
        column = [s[i].upper() for s in alignment]
        bases = [b for b in column if b in "ACGT"]
        if not bases:
            out.append("N")
            continue
        counts = Counter(bases)
        kept = [b for b, n in counts.items() if n / len(bases) >= min_frequency]
        out.append(code_for(kept or [counts.most_common(1)[0][0]]))

    return "".join(out)
