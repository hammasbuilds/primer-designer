"""Design primers against real sequenced strains, not a generated one.

`demo.py` builds its alignment from a fixed random seed. That is right for a demonstration -
it is reproducible, it needs no network, and it makes the *mechanism* obvious. It is wrong
as evidence, because a synthetic alignment has whatever conservation the generator gave it,
and the number this project turns on is how conserved a real virus actually is.

The input here is 210 Cotton leaf curl Multan virus haplotypes from NCBI GenBank, aligned
and committed to `data/`. Two things about how that set was built matter more than its size:

**One sequence per haplotype, not one per record.** 250 records collapse to 210 haplotypes;
a single 2021 Punjab submission contributes eight identical genomes. Conservation counted
per record would let that submission vote eight times on how safe a site is, and a primer
is then designed against a number inflated by its own sampling.

**One species.** Mixing Multan with Kokhran and Burewala virus measures the distance between
species and reports it as within-species variation - which would make every site look
variable and every primer look impossible.

    python scripts/real_data.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from primer.conservation import analyse  # noqa: E402
from primer.design import generate, pair_candidates  # noqa: E402
from primer.validate import coverage  # noqa: E402

FASTA = ROOT / "data" / "clcuv_aligned.fasta"


def read_fasta(path: Path) -> tuple[list[str], list[str]]:
    """Returns (headers, sequences). Lines starting with ';' are provenance, not data."""
    headers: list[str] = []
    seqs: list[str] = []
    current: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(";") or not line.strip():
            continue
        if line.startswith(">"):
            if current:
                seqs.append("".join(current))
                current = []
            headers.append(line[1:].strip())
        else:
            current.append(line.strip())
    if current:
        seqs.append("".join(current))
    return headers, seqs


def main() -> int:
    if not FASTA.exists():
        print(f"missing {FASTA}", file=sys.stderr)
        print("Regenerate it from clcuv-surveillance:", file=sys.stderr)
        print(
            "  uv run python scripts/export_alignment.py "
            "../primer-designer/data/clcuv_aligned.fasta",
            file=sys.stderr,
        )
        return 2

    headers, alignment = read_fasta(FASTA)
    width = len(alignment[0])
    print("INPUT")
    print(f"   {len(alignment)} CLCuV haplotypes x {width} aligned columns (NCBI GenBank)")
    print(f"   e.g. {headers[0]}")
    print()

    positions = analyse(alignment)
    fully = sum(1 for p in positions if p.conservation >= 0.999)
    over98 = sum(1 for p in positions if p.conservation >= 0.98)
    gapped = sum(1 for p in positions if p.gaps > len(alignment) * 0.1)
    print("CONSERVATION")
    print(f"   {fully:>5} of {width} columns are invariant across all {len(alignment)}")
    print(f"   {over98:>5} are >= 98% conserved")
    print(f"   {gapped:>5} are gapped in more than 10% of strains")
    print()

    for threshold in (0.99, 0.98, 0.95):
        candidates = generate(alignment, positions, min_conservation=threshold)
        pairs = pair_candidates(candidates, strains=alignment, limit=3)
        print(f"DESIGN at min_conservation={threshold:.2f}")
        print(f"   {len(candidates)} viable candidates, {len(pairs)} usable pairs")
        for pair in pairs[:2]:
            fwd = coverage(pair.forward.sequence, alignment)
            rev = coverage(pair.reverse.sequence, alignment, reverse=True)
            # An assay needs BOTH primers to bind. Reporting them separately, and the
            # pair as the weaker of the two, avoids the flattering version where a 99%
            # forward hides a 90% reverse.
            print(
                f"     {pair.forward.sequence:24} {pair.reverse.sequence:24} "
                f"amp={pair.amplicon_length:>4}  fwd={fwd.coverage:.1%} rev={rev.coverage:.1%} "
                f"pair={min(fwd.coverage, rev.coverage):.1%}"
            )
        print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
