"""In-silico validation: will this primer still find the virus?

The mechanism the whole project turns on, and the reason a diagnostic quietly stops
working:

**A mismatch at the 3' end is not like a mismatch anywhere else.**

Taq polymerase extends from the 3' terminus. A mismatch in the middle of a primer
lowers the melting temperature a little and the reaction usually still works. A
mismatch at the terminal base blocks extension almost completely — amplification fails
while the primer still binds, so the assay reports a clean negative.

That asymmetry is why a virus can drift a single nucleotide and render a validated PCR
test blind, with no error, no warning, and no change in the control.

So mismatches are weighted by distance from the 3' end rather than counted. A primer
with three mismatches at the 5' end is fine. A primer with one at position −1 is dead.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

from .conservation import matches
from .thermo import reverse_complement

# Penalty by distance from the 3' end. Index 0 is the terminal base.
# Steep on purpose: this curve *is* the biology, and flattening it would produce a
# tool that cheerfully approves primers which cannot extend.
POSITION_PENALTY = (10.0, 6.0, 3.0, 1.5, 1.0)
DISTAL_PENALTY = 0.3

# Above this, extension is considered blocked.
BLIND_THRESHOLD = 8.0


def mismatch_penalty(distance_from_3prime: int) -> float:
    if distance_from_3prime < len(POSITION_PENALTY):
        return POSITION_PENALTY[distance_from_3prime]
    return DISTAL_PENALTY


@dataclass
class Alignment:
    target: str
    mismatches: list[int] = field(default_factory=list)  # distances from the 3' end
    penalty: float = 0.0
    found: bool = True

    @property
    def blind(self) -> bool:
        """Would this primer fail to amplify this target?"""
        return not self.found or self.penalty >= BLIND_THRESHOLD

    @property
    def terminal_mismatch(self) -> bool:
        return 0 in self.mismatches


def score_against(primer: str, target_region: str) -> Alignment:
    """Compare a primer against the sequence it is meant to bind.

    `target_region` must already be the correct strand and length. Degenerate primer
    codes match any of their bases.
    """
    if len(primer) != len(target_region):
        raise ValueError("primer and target region must be the same length")

    mismatches = []
    penalty = 0.0
    last = len(primer) - 1

    for i, (code, base) in enumerate(zip(primer.upper(), target_region.upper(), strict=False)):
        if base not in "ACGT" or not matches(code, base):
            distance = last - i
            mismatches.append(distance)
            penalty += mismatch_penalty(distance)

    return Alignment(target=target_region, mismatches=mismatches, penalty=round(penalty, 3))


# Above this fraction of mismatched bases, the primer does not anneal at all and the
# position-weighted penalty stops being meaningful.
MAX_MISMATCH_FRACTION = 0.3


def scan_strain(
    primer: str,
    strain: str,
    *,
    reverse: bool = False,
    max_mismatch_fraction: float = MAX_MISMATCH_FRACTION,
) -> Alignment:
    """Find the best binding site for a primer in a whole strain sequence.

    Best means lowest penalty, not fewest mismatches: a site with one terminal
    mismatch is worse than a site with three internal ones, and ranking by count picks
    the wrong site.

    But the penalty model only applies to sequences that *anneal*. Sites with more
    than `max_mismatch_fraction` mismatched bases are discarded before ranking —
    otherwise a stretch of unrelated sequence with sixteen mismatches can accumulate a
    lower weighted penalty than the true site carrying one terminal mismatch, and the
    report names a binding site that does not exist. That is a real failure this
    project hit: coverage was right by accident while the diagnosis was nonsense.

    A primer with no plausible site does not bind, which is `found=False` rather than
    a low-quality match.
    """
    probe = reverse_complement(primer) if reverse else primer
    if len(probe) > len(strain):
        return Alignment(target="", found=False, penalty=float("inf"))

    limit = max(1, int(len(probe) * max_mismatch_fraction))
    best: Alignment | None = None

    for start in range(len(strain) - len(probe) + 1):
        candidate = score_against(probe, strain[start : start + len(probe)])
        if len(candidate.mismatches) > limit:
            continue
        if best is None or candidate.penalty < best.penalty:
            best = candidate
            if best.penalty == 0.0:
                break

    return best or Alignment(target="", found=False, penalty=float("inf"))


@dataclass
class CoverageReport:
    primer: str
    strains_tested: int
    perfect: int
    tolerated: int
    blind: int
    terminal_mismatches: int
    failing_strains: list[str] = field(default_factory=list)

    @property
    def coverage(self) -> float:
        """Fraction of strains this primer can still amplify."""
        if not self.strains_tested:
            return 0.0
        return round((self.strains_tested - self.blind) / self.strains_tested, 6)

    def verdict(self) -> str:
        if self.coverage >= 0.99:
            return "healthy"
        if self.coverage >= 0.95:
            return "degrading - monitor"
        if self.coverage >= 0.80:
            return "failing - redesign"
        return "blind - withdraw"


def coverage(
    primer: str,
    strains: Sequence[str],
    *,
    names: Sequence[str] | None = None,
    reverse: bool = False,
) -> CoverageReport:
    """How much of the observed population this primer still detects.

    This is the primer-health check. Run monthly against new sequences, it says when a
    deployed assay has started going blind — *before* the field notices that the tests
    have stopped finding a virus that is plainly still there.
    """
    names = list(names) if names else [f"strain-{i}" for i in range(len(strains))]

    perfect = tolerated = blind = terminal = 0
    failing: list[str] = []

    for name, strain in zip(names, strains, strict=False):
        alignment = scan_strain(primer, strain, reverse=reverse)
        if alignment.terminal_mismatch:
            terminal += 1
        if alignment.blind:
            blind += 1
            failing.append(name)
        elif not alignment.mismatches:
            perfect += 1
        else:
            tolerated += 1

    return CoverageReport(
        primer=primer,
        strains_tested=len(strains),
        perfect=perfect,
        tolerated=tolerated,
        blind=blind,
        terminal_mismatches=terminal,
        failing_strains=failing,
    )


def primer_health_alert(
    report: CoverageReport,
    *,
    previous_coverage: float | None = None,
    drop_threshold: float = 0.02,
) -> dict | None:
    """Raise an alert when a deployed primer's coverage falls.

    The *trend* matters more than the level. Coverage sliding from 99% to 96% over a
    season is a virus actively escaping the assay, and it deserves attention well
    before the absolute number looks alarming.
    """
    alerts = []

    if report.coverage < 0.95:
        alerts.append(f"coverage {report.coverage:.1%} - {report.verdict()}")

    if report.terminal_mismatches:
        share = report.terminal_mismatches / max(report.strains_tested, 1)
        alerts.append(
            f"{share:.1%} of strains carry a 3' terminal mismatch - "
            "extension is blocked, not merely weakened"
        )

    if previous_coverage is not None:
        drop = previous_coverage - report.coverage
        if drop >= drop_threshold:
            alerts.append(
                f"coverage fell {drop:.1%} since the last check "
                f"({previous_coverage:.1%} -> {report.coverage:.1%})"
            )

    if not alerts:
        return None

    return {
        "primer": report.primer,
        "coverage": report.coverage,
        "verdict": report.verdict(),
        "alerts": alerts,
        "failing_strains": report.failing_strains[:10],
    }
