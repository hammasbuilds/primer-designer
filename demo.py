"""Design a diagnostic primer pair, then watch a single mutation blind it.

    python demo.py

Two halves. First: an aligned set of strains in, a usable primer pair out.
Second: the same primer checked against a drifted population, where one
nucleotide at the 3' end is the difference between working and useless.
No network, no dependencies.
"""
import random
import sys

sys.path.insert(0, "src")

from primer.conservation import analyse
from primer.design import generate, pair_candidates
from primer.validate import coverage, primer_health_alert

# A conserved viral region as it appears in 12 sequenced isolates. Generated from a
# fixed seed rather than written by hand: a hand-written sequence tends to repeat its
# own motifs, and a repeated motif makes coverage meaningless -- a primer knocked out
# in one copy still matches perfectly in another.
rng = random.Random(4)
CORE = "".join(rng.choice("ACGT") for _ in range(189))
alignment = [CORE] * 12

print("INPUT")
print(f"   alignment          {len(alignment)} strains x {len(CORE)} bp")
print(f"   {CORE[:60]}...")
print()

positions = analyse(alignment)
candidates = generate(alignment, positions, min_conservation=0.98)
pairs = pair_candidates(candidates, strains=alignment, limit=3)

print("OUTPUT - design")
print(f"   {len(candidates)} viable candidates, {len(pairs)} usable pairs")
print()
print(f"   {'forward':26} {'reverse':26} {'amp':>5} {'Tm_f':>6} {'Tm_r':>6} {'dTm':>5}")
for p in pairs:
    s = p.summary()
    print(f"   {s['forward']:26} {s['reverse']:26} {s['amplicon']:>5} "
          f"{s['tm_forward']:>6} {s['tm_reverse']:>6} {p.tm_difference:>5.2f}")

best = pairs[0]
primer = best.forward.sequence
print()
print(f"   deployed assay: forward primer {primer}")

# --- the virus drifts -------------------------------------------------------
# Two escape routes. One mutation sits mid-primer, the other on the final base,
# where a polymerase cannot extend at all.
start = best.forward.start
mid = start + best.forward.length // 2
last = start + best.forward.length - 1


def mutate(seq: str, at: int) -> str:
    swap = {"A": "G", "G": "A", "C": "T", "T": "C"}[seq[at]]
    return seq[:at] + swap + seq[at + 1:]


population = [CORE] * 6 + [mutate(CORE, mid)] * 3 + [mutate(CORE, last)] * 3
names = ["wildtype"] * 6 + ["mid-primer-mismatch"] * 3 + ["3prime-mismatch"] * 3

report = coverage(primer, population, names=names)

print()
print("OUTPUT - six months later, same primer, drifted population")
print(f"   strains tested     {report.strains_tested}")
print(f"   coverage           {report.coverage:.1%}   ({report.verdict()})")
print(f"   3' terminal misses {report.terminal_mismatches}")
print()
alert = primer_health_alert(report, previous_coverage=1.0)
if alert:
    print("   ALERT")
    for line in alert["alerts"]:
        print(f"      - {line}")
    print(f"      failing: {', '.join(sorted(set(alert['failing_strains'])))}")
else:
    print("   no alert raised")
print()
print("   The three mid-primer mismatches were tolerated and still amplify.")
print("   The three 3' terminal mismatches are the whole 25% loss.")
