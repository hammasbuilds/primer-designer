"""Oligonucleotide thermodynamics.

Melting temperature by the nearest-neighbour model (SantaLucia's unified parameters,
1998), because the rules of thumb do not survive contact with a real assay:

    Tm ≈ 2(A+T) + 4(G+C)        the Wallace rule — off by 10 °C on a 25-mer
    Tm from %GC alone           ignores that GC *stacking* differs by neighbour

The nearest-neighbour model works on stacked base pairs rather than on individual
bases, because the stability of a duplex comes from stacking interactions between
adjacent pairs, not from the pairs in isolation. `GC` next to `GC` is not worth the
same as `GC` next to `AT`.

    ΔG° = ΔH° − T·ΔS°
    Tm  = ΔH° / (ΔS° + R·ln(C_T/4)) − 273.15

Everything here is arithmetic over a published table, so it is testable exactly —
which matters, because a primer that melts 6 °C from its partner will amplify
inconsistently and the failure looks like bad technique rather than bad design.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# SantaLucia (1998) unified nearest-neighbour parameters.
# ΔH in kcal/mol, ΔS in cal/(mol·K). Keys read 5'->3' on the top strand.
NN_PARAMS: dict[str, tuple[float, float]] = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}

# Helix initiation depends on which base pair sits at each end.
INIT_GC = (0.1, -2.8)
INIT_AT = (2.3, 4.1)

R = 1.987  # cal/(mol·K)

COMPLEMENT = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


class SequenceError(ValueError):
    pass


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(COMPLEMENT)[::-1]


def gc_content(sequence: str) -> float:
    if not sequence:
        return 0.0
    seq = sequence.upper()
    return round(sum(1 for c in seq if c in "GC") / len(seq), 6)


def has_gc_clamp(sequence: str) -> bool:
    """A G or C in the last two bases of the 3' end.

    The 3' end is where polymerase extends from, so a G·C pair there — three hydrogen
    bonds rather than two — holds the primer in place while extension begins. Its
    absence is the commonest cause of a primer that works in silico and not on a
    bench.
    """
    return bool(sequence) and any(c in "GC" for c in sequence.upper()[-2:])


def has_gc_run(sequence: str, *, limit: int = 4) -> bool:
    """More than `limit` consecutive G or C.

    A GC run binds so tightly that the primer anneals at partially matching sites and
    stays there, which produces non-specific product.
    """
    run = 0
    for base in sequence.upper():
        run = run + 1 if base in "GC" else 0
        if run > limit:
            return True
    return False


@dataclass(frozen=True)
class Thermo:
    tm: float
    delta_h: float
    delta_s: float
    gc: float

    def delta_g(self, celsius: float = 37.0) -> float:
        kelvin = celsius + 273.15
        return round(self.delta_h - kelvin * self.delta_s / 1000, 4)


def melting_temperature(
    sequence: str,
    *,
    primer_nm: float = 500.0,
    na_mm: float = 50.0,
    mg_mm: float = 0.0,
) -> Thermo:
    """Nearest-neighbour Tm for a self-complementary-free oligo.

    Defaults are ordinary PCR conditions. Salt matters: a Tm computed at 1 M sodium
    and used at 50 mM is wrong by more than 10 °C, and 50 mM is what is actually in
    the tube.
    """
    seq = sequence.upper().replace("U", "T")
    if len(seq) < 2:
        raise SequenceError("sequence must be at least 2 bases")
    if set(seq) - set("ACGT"):
        raise SequenceError(
            f"unsupported bases {sorted(set(seq) - set('ACGT'))}; "
            "expand degenerate codes before computing Tm"
        )

    delta_h, delta_s = 0.0, 0.0
    for i in range(len(seq) - 1):
        pair = seq[i : i + 2]
        if pair not in NN_PARAMS:
            raise SequenceError(f"no parameters for {pair!r}")
        dh, ds = NN_PARAMS[pair]
        delta_h += dh
        delta_s += ds

    for terminal in (seq[0], seq[-1]):
        dh, ds = INIT_GC if terminal in "GC" else INIT_AT
        delta_h += dh
        delta_s += ds

    # Monovalent equivalent. Magnesium stabilises far more strongly than sodium; the
    # 4·sqrt([Mg]) equivalence is the standard working approximation.
    sodium = na_mm / 1000 + (4 * math.sqrt(mg_mm / 1000) if mg_mm > 0 else 0.0)
    sodium = max(sodium, 1e-6)
    # von Ahsen salt correction, applied to entropy.
    delta_s_corrected = delta_s + 0.368 * (len(seq) - 1) * math.log(sodium)

    concentration = primer_nm * 1e-9
    denominator = delta_s_corrected + R * math.log(concentration / 4)
    if denominator == 0:
        raise SequenceError("degenerate thermodynamics: entropy term cancels")

    tm = (delta_h * 1000) / denominator - 273.15
    return Thermo(
        tm=round(tm, 3),
        delta_h=round(delta_h, 3),
        delta_s=round(delta_s_corrected, 3),
        gc=gc_content(seq),
    )


def self_complementarity(sequence: str, *, min_run: int = 4) -> int:
    """Longest self-complementary run — the primer binding itself.

    A self-dimer consumes primer and produces a product that looks like contamination.
    Reported as a length rather than a boolean so a design can trade it off.
    """
    seq = sequence.upper()
    rc = reverse_complement(seq)
    best = 0
    for length in range(min_run, len(seq) + 1):
        for i in range(len(seq) - length + 1):
            if seq[i : i + length] in rc:
                best = max(best, length)
    return best


def three_prime_dimer(sequence: str, *, window: int = 5) -> int:
    """Self-complementarity involving the 3' end specifically.

    Worse than an internal dimer by a wide margin: a 3' dimer is *extendable*, so the
    polymerase amplifies primer-on-primer and the reaction consumes itself before it
    reaches the template.
    """
    seq = sequence.upper()
    tail = seq[-window:]
    rc = reverse_complement(seq)
    best = 0
    for length in range(2, len(tail) + 1):
        if tail[-length:] in rc:
            best = length
    return best


def hairpin_stem(sequence: str, *, min_loop: int = 3, min_stem: int = 4) -> int:
    """Longest stem a hairpin could form, allowing for a minimum loop.

    A hairpin sequesters the primer in a fold that has to melt before annealing, which
    shows up as poor amplification at exactly the annealing temperature that was
    calculated to be correct.
    """
    seq = sequence.upper()
    n = len(seq)
    best = 0
    for stem in range(min_stem, n // 2 + 1):
        for i in range(0, n - 2 * stem - min_loop + 1):
            left = seq[i : i + stem]
            for j in range(i + stem + min_loop, n - stem + 1):
                right = seq[j : j + stem]
                if left == reverse_complement(right):
                    best = max(best, stem)
    return best
