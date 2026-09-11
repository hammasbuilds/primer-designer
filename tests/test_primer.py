"""Primer designer tests.

Thermodynamics is arithmetic over a published table and mismatch penalties are a fixed
curve, so most of this is exact. Where a property is easier to state than a value —
"GC-rich melts hotter than AT-rich" — the property is tested rather than a memorised
constant, which is the honest way round.
"""

from __future__ import annotations

import pytest

from primer.conservation import (
    AlignmentError,
    analyse,
    code_for,
    conserved_windows,
    degeneracy,
    degenerate_consensus,
    matches,
)
from primer.design import Constraints, evaluate, generate, pair_candidates
from primer.thermo import (
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
from primer.validate import (
    coverage,
    mismatch_penalty,
    primer_health_alert,
    scan_strain,
    score_against,
)

CORE = "ATGGCTAAGCGTCCAGGATCCAAGTTCGATCCGTTAACGGCTAAGGCATCGTAGCTAGCATCG"


class TestSequenceBasics:
    def test_reverse_complement(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_reverse_complement_is_an_involution(self):
        assert reverse_complement(reverse_complement(CORE)) == CORE

    def test_gc_content(self):
        assert gc_content("GGCC") == 1.0
        assert gc_content("ATAT") == 0.0
        assert gc_content("ATGC") == 0.5

    def test_gc_clamp_looks_at_the_three_prime_end(self):
        """Where the polymerase extends from."""
        assert has_gc_clamp("AAAAG")
        assert not has_gc_clamp("GGGAA")

    def test_gc_runs_are_detected(self):
        assert has_gc_run("AAGGCCGAA")
        assert not has_gc_run("AAGGCAGGA")


class TestThermodynamics:
    def test_gc_rich_melts_hotter_than_at_rich(self):
        hot = melting_temperature("GCGCGCGCGCGCGCGCGCGC").tm
        cold = melting_temperature("ATATATATATATATATATAT").tm
        assert hot > cold + 30

    def test_a_typical_primer_lands_in_the_usable_range(self):
        """If a 20-mer at 50% GC did not come out near 55-60 °C, the parameters are
        wrong somewhere."""
        tm = melting_temperature("TTGACCTGCAAGTGCATCGA").tm
        assert 50 < tm < 65

    def test_longer_is_hotter(self):
        short = melting_temperature("TTGACCTGCAAGTGC").tm
        long = melting_temperature("TTGACCTGCAAGTGCATCGATCGAT").tm
        assert long > short

    def test_salt_raises_melting_temperature(self):
        """A Tm computed at 1 M sodium and used at 50 mM is wrong by over 10 °C, and
        50 mM is what is in the tube."""
        low = melting_temperature("TTGACCTGCAAGTGCATCGA", na_mm=50).tm
        high = melting_temperature("TTGACCTGCAAGTGCATCGA", na_mm=500).tm
        assert high > low + 5

    def test_magnesium_stabilises_more_strongly_than_sodium(self):
        without = melting_temperature("TTGACCTGCAAGTGCATCGA", na_mm=50).tm
        with_mg = melting_temperature("TTGACCTGCAAGTGCATCGA", na_mm=50, mg_mm=3).tm
        assert with_mg > without

    def test_nearest_neighbour_beats_a_base_count(self):
        """Same composition, different order, different Tm — which is the entire
        reason for using a stacking model."""
        a = melting_temperature("GGGGCCCCAAAATTTT").tm
        b = melting_temperature("GCGCGCGCATATATAT").tm
        assert a != b

    def test_delta_g_is_negative_for_a_stable_duplex(self):
        assert melting_temperature("TTGACCTGCAAGTGCATCGA").delta_g(37.0) < 0

    def test_degenerate_bases_are_refused(self):
        """Rather than silently guessing a resolution."""
        with pytest.raises(SequenceError):
            melting_temperature("ATGCRYATGC")

    def test_a_single_base_is_refused(self):
        with pytest.raises(SequenceError):
            melting_temperature("A")


class TestSecondaryStructure:
    def test_a_hairpin_is_detected(self):
        assert hairpin_stem("GGGGAAAACCCC") >= 4

    def test_a_clean_sequence_has_no_hairpin(self):
        assert hairpin_stem("AAAGAAAGAAAGAAAG") == 0

    def test_self_complementarity_is_measured(self):
        assert self_complementarity("ACGTACGT") >= 4

    def test_a_three_prime_dimer_is_detected(self):
        """Worse than an internal dimer: it is extendable, so the reaction amplifies
        primer-on-primer and consumes itself."""
        assert three_prime_dimer("AAAAAAAAAAGCGCGC") >= 2


class TestIUPAC:
    @pytest.mark.parametrize(
        ("bases", "code"),
        [
            ("AG", "R"),
            ("CT", "Y"),
            ("GC", "S"),
            ("AT", "W"),
            ("GT", "K"),
            ("AC", "M"),
            ("ACGT", "N"),
            ("A", "A"),
        ],
    )
    def test_codes_cover_exactly_their_bases(self, bases, code):
        assert code_for(list(bases)) == code

    def test_matching_is_by_code(self):
        assert matches("R", "A") and matches("R", "G")
        assert not matches("R", "C")

    def test_degeneracy_is_multiplicative(self):
        """The cost made explicit: six N positions is 4096 oligos sharing one tube,
        each at 1/4096 of the nominal concentration."""
        assert degeneracy("NNN") == 64
        assert degeneracy("ACGT") == 1
        assert degeneracy("RYACGT") == 4

    def test_an_empty_base_set_becomes_n(self):
        assert code_for([]) == "N"


class TestConservation:
    def test_an_invariant_column_scores_one(self):
        assert analyse(["ACGT", "ACGT", "ACGT"])[0].conservation == 1.0

    def test_a_variable_column_is_found(self):
        positions = analyse(["ACGT", "ACGT", "TCGT"])
        assert positions[0].conservation == pytest.approx(2 / 3)
        assert positions[0].counts == {"A": 2, "T": 1}

    def test_entropy_is_zero_when_invariant(self):
        assert analyse(["ACGT"] * 5)[0].entropy == 0.0

    def test_entropy_is_maximal_when_uniform(self):
        assert analyse(["AC", "CC", "GC", "TC"])[0].entropy == pytest.approx(2.0)

    def test_gaps_are_missing_data_not_disagreement(self):
        """A gap is an absence of observation. Counting it as a difference makes every
        poorly sequenced region look variable."""
        positions = analyse(["ACGT", "-CGT", "ACGT"])
        assert positions[0].conservation == 1.0
        assert positions[0].gaps == 1

    def test_unaligned_sequences_are_refused(self):
        """Column-wise analysis of unaligned sequences is meaningless, not merely
        inaccurate."""
        with pytest.raises(AlignmentError):
            analyse(["ACGT", "ACG"])

    def test_empty_input_is_refused(self):
        with pytest.raises(AlignmentError):
            analyse([])

    def test_a_window_needs_every_position_conserved(self):
        """One variable position under a primer is enough to break it, and a window
        mean hides it behind nineteen invariant neighbours."""
        alignment = ["A" * 10 + "C" + "A" * 10, "A" * 10 + "T" + "A" * 10]
        positions = analyse(alignment)
        windows = conserved_windows(positions, length=5, min_conservation=0.99)
        assert all(not (start <= 10 < start + 5) for start, _ in windows)

    def test_degenerate_consensus_covers_real_variation(self):
        alignment = ["ACGT", "ACGT", "AGGT", "AGGT"]
        assert degenerate_consensus(alignment, 0, 4) == "ASGT"

    def test_singletons_are_excluded_from_the_consensus(self):
        """A base seen once in four hundred strains is usually a sequencing error, and
        covering it costs specificity for no gain."""
        alignment = ["ACGT"] * 99 + ["ATGT"]
        assert degenerate_consensus(alignment, 0, 4, min_frequency=0.05) == "ACGT"


class TestMismatchAsymmetry:
    def test_a_terminal_mismatch_costs_far_more_than_an_internal_one(self):
        """The mechanism the project turns on: Taq extends from the 3' terminus, so a
        mismatch there blocks extension while the primer still binds — and the assay
        reports a clean negative."""
        assert mismatch_penalty(0) > mismatch_penalty(10) * 20

    def test_penalty_decreases_with_distance_from_the_three_prime_end(self):
        penalties = [mismatch_penalty(d) for d in range(6)]
        assert penalties == sorted(penalties, reverse=True)

    def test_one_terminal_mismatch_blinds_a_primer(self):
        perfect = "ACGTACGTACGTACGTACGT"
        drifted = perfect[:-1] + "A"
        assert score_against(perfect, drifted).blind

    def test_several_distal_mismatches_do_not(self):
        perfect = "ACGTACGTACGTACGTACGT"
        drifted = "TTT" + perfect[3:]
        assert not score_against(perfect, drifted).blind

    def test_degenerate_primers_match_their_variants(self):
        assert score_against("ASGT", "ACGT").mismatches == []
        assert score_against("ASGT", "AGGT").mismatches == []

    def test_length_mismatch_is_refused(self):
        with pytest.raises(ValueError):
            score_against("ACGT", "ACG")


class TestScanning:
    def test_the_true_site_is_found(self):
        assert scan_strain(CORE[5:25], CORE).mismatches == []

    def test_an_implausible_site_is_not_reported_as_a_binding_site(self):
        """A real bug this caught: an unrelated stretch with sixteen mismatches
        accumulated a *lower* weighted penalty than the true site's single terminal
        mismatch, so coverage was right by accident while the diagnosis was nonsense.
        """
        primer = CORE[5:25]
        drifted = CORE[:24] + "A" + CORE[25:]
        best = scan_strain(primer, drifted)
        assert best.terminal_mismatch
        assert len(best.mismatches) == 1

    def test_a_primer_with_no_plausible_site_does_not_bind(self):
        assert not scan_strain(CORE[5:25], "T" * 60).found

    def test_a_primer_longer_than_the_strain_does_not_bind(self):
        assert not scan_strain("ACGTACGTACGT", "ACGT").found


class TestCoverage:
    def test_a_conserved_primer_covers_everything(self):
        assert coverage(CORE[5:25], [CORE] * 20).coverage == 1.0

    def test_drift_at_the_three_prime_end_halves_coverage(self):
        primer = CORE[5:25]
        drifted = CORE[:24] + "A" + CORE[25:]
        report = coverage(primer, [CORE if i % 2 else drifted for i in range(50)])
        assert report.coverage == 0.5
        assert report.terminal_mismatches == 25
        assert report.verdict() == "blind - withdraw"

    def test_failing_strains_are_named(self):
        primer = CORE[5:25]
        drifted = CORE[:24] + "A" + CORE[25:]
        report = coverage(primer, [CORE, drifted], names=["ok", "escaped"])
        assert report.failing_strains == ["escaped"]

    def test_verdicts_are_graded(self):
        primer = CORE[5:25]
        drifted = CORE[:24] + "A" + CORE[25:]
        strains = [drifted] * 3 + [CORE] * 97
        assert coverage(primer, strains).verdict() == "degrading - monitor"

    def test_no_strains_means_no_coverage(self):
        assert coverage(CORE[5:25], []).coverage == 0.0


class TestHealthAlerts:
    def test_a_healthy_primer_raises_nothing(self):
        assert primer_health_alert(coverage(CORE[5:25], [CORE] * 50)) is None

    def test_a_falling_trend_alerts_before_the_level_looks_alarming(self):
        """Coverage sliding from 99% to 96% over a season is a virus actively
        escaping the assay."""
        primer = CORE[5:25]
        drifted = CORE[:24] + "A" + CORE[25:]
        report = coverage(primer, [drifted] * 3 + [CORE] * 97)
        alert = primer_health_alert(report, previous_coverage=0.99)
        assert alert is not None
        assert any("fell" in a for a in alert["alerts"])

    def test_the_terminal_mismatch_mechanism_is_named_in_the_alert(self):
        """An operator needs to know *why*, not just that a number moved."""
        primer = CORE[5:25]
        drifted = CORE[:24] + "A" + CORE[25:]
        alert = primer_health_alert(coverage(primer, [drifted] * 50))
        assert any("extension is blocked" in a for a in alert["alerts"])


class TestDesign:
    def test_every_constraint_failure_is_collected(self):
        """Fixing one at a time is how a design session takes a day."""
        candidate = evaluate(
            "AAAAAAAAAAAAAAAAAAAA",
            start=0,
            conservation=1.0,
            constraints=Constraints(),
        )
        assert len(candidate.rejections) >= 2
        assert not candidate.viable

    def test_a_good_candidate_is_viable(self):
        candidate = evaluate(
            "TGTCGAGCGACGGAATTAGA",
            start=0,
            conservation=1.0,
            constraints=Constraints(),
        )
        assert candidate.viable, candidate.rejections

    def test_over_degenerate_candidates_are_rejected(self):
        candidate = evaluate(
            "NNNNACCTGCAAGTGCATCG",
            start=0,
            conservation=1.0,
            constraints=Constraints(max_degeneracy=16),
        )
        assert any("degeneracy" in r for r in candidate.rejections)

    def test_candidates_are_generated_from_conserved_regions(self):
        alignment = [CORE] * 10
        candidates = generate(alignment, analyse(alignment))
        assert candidates
        assert all(c.viable for c in candidates)

    def test_pairs_respect_the_amplicon_window(self):
        alignment = [CORE] * 10
        candidates = generate(alignment, analyse(alignment))
        pairs = pair_candidates(candidates, min_amplicon=20, max_amplicon=60)
        assert all(20 <= p.amplicon_length <= 60 for p in pairs)

    def test_pairs_respect_the_melting_temperature_difference(self):
        """Two primers melting 6 °C apart cannot share an annealing temperature."""
        alignment = [CORE] * 10
        candidates = generate(alignment, analyse(alignment))
        pairs = pair_candidates(candidates, min_amplicon=20, max_amplicon=60, max_tm_difference=2.0)
        assert all(p.tm_difference <= 2.0 for p in pairs)

    def test_a_pair_summary_is_reportable(self):
        alignment = [CORE] * 10
        candidates = generate(alignment, analyse(alignment))
        pairs = pair_candidates(candidates, min_amplicon=20, max_amplicon=60)
        if pairs:
            summary = pairs[0].summary()
            assert {"forward", "reverse", "amplicon", "tm_difference"} <= set(summary)
