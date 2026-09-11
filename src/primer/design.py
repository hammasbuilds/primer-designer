"""Candidate generation, filtering and pair selection.

Design is a constrained search, not an optimisation: most of the work is *rejecting*
candidates, and the constraints come from the assay rather than from the algorithm.

A pair is scored on four things, and the ordering between them is a judgement worth
stating. Coverage comes first — a primer that no longer matches the pathogen is worth
nothing however elegant its thermodynamics. Tm match comes second, because a mismatched
pair amplifies inconsistently in a way that looks like operator error. Then structural
cleanliness, then degeneracy, which is a cost to be minimised rather than a feature.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

from .conservation import Position, conserved_windows, degenerate_consensus, degeneracy
from .thermo import (
    SequenceError,
    gc_content,
    hairpin_stem,
    has_gc_clamp,
    has_gc_run,
    melting_temperature,
    reverse_complement,
    self_complementarity,
    three_prime_dimer,
)
from .validate import CoverageReport, coverage


@dataclass
class Constraints:
    min_length: int = 18
    max_length: int = 24
    min_tm: float = 55.0
    max_tm: float = 65.0
    min_gc: float = 0.40
    max_gc: float = 0.60
    max_degeneracy: int = 16
    max_hairpin_stem: int = 4
    max_self_complementarity: int = 5
    max_three_prime_dimer: int = 3
    require_gc_clamp: bool = True
    allow_gc_run: bool = False


@dataclass
class Candidate:
    sequence: str
    start: int
    length: int
    tm: float
    gc: float
    degeneracy: int
    conservation: float
    is_reverse: bool = False
    rejections: list[str] = field(default_factory=list)

    @property
    def viable(self) -> bool:
        return not self.rejections


def evaluate(
    sequence: str, *, start: int, conservation: float, constraints: Constraints,
    is_reverse: bool = False,
) -> Candidate:
    """Check one candidate against every constraint.

    Every failure is collected rather than returning on the first. A designer needs to
    know that a candidate failed on Tm *and* GC content — fixing one at a time is how
    a design session takes a day.
    """
    rejections: list[str] = []
    deg = degeneracy(sequence)

    if not constraints.min_length <= len(sequence) <= constraints.max_length:
        rejections.append(f"length {len(sequence)} outside "
                          f"{constraints.min_length}-{constraints.max_length}")

    if deg > constraints.max_degeneracy:
        rejections.append(f"degeneracy {deg} above {constraints.max_degeneracy}")

    gc = gc_content(sequence.replace("N", ""))
    if not constraints.min_gc <= gc <= constraints.max_gc:
        rejections.append(f"GC {gc:.0%} outside "
                          f"{constraints.min_gc:.0%}-{constraints.max_gc:.0%}")

    if constraints.require_gc_clamp and not has_gc_clamp(sequence):
        rejections.append("no GC clamp at the 3' end")

    if not constraints.allow_gc_run and has_gc_run(sequence):
        rejections.append("run of 5 or more G/C")

    if (stem := hairpin_stem(sequence)) > constraints.max_hairpin_stem:
        rejections.append(f"hairpin stem {stem} above {constraints.max_hairpin_stem}")

    if (sc := self_complementarity(sequence)) > constraints.max_self_complementarity:
        rejections.append(f"self-complementarity {sc} above "
                          f"{constraints.max_self_complementarity}")

    if (tp := three_prime_dimer(sequence)) > constraints.max_three_prime_dimer:
        rejections.append(f"3' dimer {tp} - extendable primer-dimer")

    tm = 0.0
    try:
        # Tm is computed on the most-common resolution of a degenerate primer; the
        # spread across variants is small relative to the constraint window.
        resolved = "".join(
            b if b in "ACGT" else {"R": "A", "Y": "C", "S": "G", "W": "A", "K": "G",
                                   "M": "A", "B": "C", "D": "A", "H": "A",
                                   "V": "A", "N": "A"}[b]
            for b in sequence.upper()
        )
        tm = melting_temperature(resolved).tm
        if not constraints.min_tm <= tm <= constraints.max_tm:
            rejections.append(f"Tm {tm:.1f} outside "
                              f"{constraints.min_tm}-{constraints.max_tm}")
    except (SequenceError, KeyError) as exc:
        rejections.append(f"thermodynamics unavailable: {exc}")

    return Candidate(
        sequence=sequence, start=start, length=len(sequence), tm=tm, gc=gc,
        degeneracy=deg, conservation=conservation, is_reverse=is_reverse,
        rejections=rejections,
    )


def generate(
    alignment: Sequence[str], positions: Sequence[Position], *,
    constraints: Constraints | None = None, min_conservation: float = 0.98,
) -> list[Candidate]:
    """All viable candidates across every conserved window and allowed length."""
    constraints = constraints or Constraints()
    seen: set[tuple[int, int]] = set()
    out: list[Candidate] = []

    for length in range(constraints.min_length, constraints.max_length + 1):
        for start, mean_conservation in conserved_windows(
            positions, length=length, min_conservation=min_conservation
        ):
            if (start, length) in seen:
                continue
            seen.add((start, length))
            sequence = degenerate_consensus(alignment, start, length)
            candidate = evaluate(
                sequence, start=start, conservation=mean_conservation,
                constraints=constraints,
            )
            if candidate.viable:
                out.append(candidate)

    return sorted(out, key=lambda c: (-c.conservation, c.degeneracy))


@dataclass
class PrimerPair:
    forward: Candidate
    reverse: Candidate
    amplicon_length: int
    tm_difference: float
    score: float
    forward_coverage: CoverageReport | None = None
    reverse_coverage: CoverageReport | None = None

    def summary(self) -> dict:
        return {
            "forward": self.forward.sequence,
            "reverse": self.reverse.sequence,
            "amplicon": self.amplicon_length,
            "tm_forward": round(self.forward.tm, 2),
            "tm_reverse": round(self.reverse.tm, 2),
            "tm_difference": round(self.tm_difference, 2),
            "degeneracy": self.forward.degeneracy * self.reverse.degeneracy,
            "score": round(self.score, 4),
            "forward_coverage": (
                self.forward_coverage.coverage if self.forward_coverage else None
            ),
            "reverse_coverage": (
                self.reverse_coverage.coverage if self.reverse_coverage else None
            ),
        }


def pair_candidates(
    candidates: Sequence[Candidate], *, min_amplicon: int = 100,
    max_amplicon: int = 400, max_tm_difference: float = 2.0,
    strains: Sequence[str] | None = None, limit: int = 5,
) -> list[PrimerPair]:
    """Pair forward and reverse candidates into usable assays.

    `max_tm_difference` defaults to 2 °C. Two primers that melt 6 °C apart cannot share
    an annealing temperature: one is too loose and amplifies non-specifically, or the
    other is too tight and does not amplify at all.
    """
    pairs: list[PrimerPair] = []

    for forward in candidates:
        for reverse_source in candidates:
            amplicon = reverse_source.start + reverse_source.length - forward.start
            if not min_amplicon <= amplicon <= max_amplicon:
                continue

            reverse = Candidate(
                sequence=reverse_complement(reverse_source.sequence),
                start=reverse_source.start, length=reverse_source.length,
                tm=reverse_source.tm, gc=reverse_source.gc,
                degeneracy=reverse_source.degeneracy,
                conservation=reverse_source.conservation, is_reverse=True,
            )

            tm_difference = abs(forward.tm - reverse.tm)
            if tm_difference > max_tm_difference:
                continue

            forward_cov = coverage(forward.sequence, strains) if strains else None
            reverse_cov = (
                coverage(reverse_source.sequence, strains) if strains else None
            )

            # Coverage first: a primer that no longer matches the pathogen is worth
            # nothing however good its thermodynamics.
            coverage_term = (
                min(forward_cov.coverage, reverse_cov.coverage)
                if forward_cov and reverse_cov else 1.0
            )
            conservation_term = (forward.conservation + reverse.conservation) / 2
            tm_term = 1 - tm_difference / max_tm_difference
            degeneracy_term = 1 / (1 + (forward.degeneracy * reverse.degeneracy) / 16)

            score = (
                0.45 * coverage_term
                + 0.30 * conservation_term
                + 0.15 * tm_term
                + 0.10 * degeneracy_term
            )

            pairs.append(PrimerPair(
                forward=forward, reverse=reverse, amplicon_length=amplicon,
                tm_difference=tm_difference, score=score,
                forward_coverage=forward_cov, reverse_coverage=reverse_cov,
            ))

    pairs.sort(key=lambda p: -p.score)
    return pairs[:limit]
