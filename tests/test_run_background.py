"""Pure logic of the run-background measure: contrast against the run's noise."""

import pytest
from run_background import classify_variant, alt_fraction, summarise_noise


# --- which counter does this variant read? ---------------------------------

@pytest.mark.parametrize("ref,alt,kind", [
    ("A", "T", "snv"),
    ("C", "G", "snv"),
    ("GAAAAA", "G", "del"),      # c.2626-8_2626-4del
    ("A", "AT", "ins"),          # dup / insertion
    ("AT", "GC", "mnv"),         # not counted: no single counter fits
])
def test_classify_variant(ref, alt, kind):
    assert classify_variant(ref, alt) == kind


# --- the alt fraction reads the matching counter ----------------------------

def test_alt_fraction_uses_the_alt_base_for_snv():
    counts = {"depth": 500, "A": 10, "C": 0, "G": 490, "T": 0, "ins": 0, "del": 3}
    assert alt_fraction(counts, "G", "A") == pytest.approx(0.02)


def test_alt_fraction_uses_the_del_counter_for_deletions():
    counts = {"depth": 552, "A": 1, "C": 0, "G": 551, "T": 0, "ins": 0, "del": 13}
    assert alt_fraction(counts, "GAAA", "G") == pytest.approx(13 / 552)


def test_alt_fraction_is_none_without_depth():
    """No reads is not zero noise — it is no measurement."""
    assert alt_fraction({"depth": 0, "A": 0, "C": 0, "G": 0, "T": 0,
                         "ins": 0, "del": 0}, "G", "A") is None


# --- summarising the run ----------------------------------------------------

def _counts(depth, **kw):
    base = {"depth": depth, "A": 0, "C": 0, "G": 0, "T": 0, "ins": 0, "del": 0}
    return {**base, **kw}


def test_clean_site_yields_a_high_contrast():
    """Variant_htz_propre: patient 50%, no other sample carries the allele."""
    s = summarise_noise(
        patient=_counts(1516, T=751, C=764),
        others={"S2": _counts(396, C=395, A=1)},
        ref="C", alt="T", min_other_samples=1,
    )
    assert s["n_with_alt"] == 0
    assert s["contrast"] > 100


def test_shared_noise_collapses_the_contrast():
    """Artefact_called_chez_S3_pas_chez_S2: 2.4% against 2.0% — a ratio near 1."""
    s = summarise_noise(
        patient=_counts(552, G=551, A=1, **{"del": 13}),
        others={"S2": _counts(556, G=556, **{"del": 11})},
        ref="GA", alt="G", min_other_samples=1,
    )
    assert s["contrast"] == pytest.approx(1.19, abs=0.05)


def test_genotype_carriers_are_excluded_from_the_noise():
    """Another sample at 45% carries the variant; it is not background noise."""
    s = summarise_noise(
        patient=_counts(500, A=250, G=250),
        others={"S2": _counts(500, A=225, G=275),      # carrier, excluded
                "S3": _counts(500, G=500)},            # truly clean
        ref="G", alt="A", genotype_floor=0.25, min_other_samples=1,
    )
    assert s["carriers"] == ["S2"]
    assert s["n_with_alt"] == 0        # S3 has none, S2 does not count as noise


def test_measurement_needs_enough_samples():
    """Below min_other_samples the median is meaningless: report nothing."""
    s = summarise_noise(patient=_counts(500, A=50, G=450),
                        others={"S2": _counts(500, G=500)},
                        ref="G", alt="A", min_other_samples=5)
    assert s["contrast"] is None


# --- clustering des positions : un parcours de BAM au lieu de N -------------

def test_positions_group_into_regions():
    from run_background import cluster_positions
    pos = [("chr9", 100), ("chr9", 150), ("chr9", 9000), ("chr16", 20)]
    assert cluster_positions(pos, gap=1000) == [
        ("chr16", 20, 20), ("chr9", 100, 150), ("chr9", 9000, 9000)]


def test_a_wide_gap_splits_the_region():
    """Bridging a 40 kb gap would pile up 40 000 useless columns."""
    from run_background import cluster_positions
    assert cluster_positions([("chr9", 1), ("chr9", 50_000)], gap=2000) == [
        ("chr9", 1, 1), ("chr9", 50_000, 50_000)]


def test_sample_comes_from_the_column_when_present(tmp_path):
    """An aggregate holds many patients in one sheet: the filename cannot say
    which row belongs to whom, only the SAMPLE column can."""
    import pandas as pd
    from openpyxl import load_workbook
    from run_background import annotate_workbooks

    p = tmp_path / "run_2_all_variants_20260903.xlsx"
    pd.DataFrame([
        {"SAMPLE": "A_S1", "CHR": "chr9", "POS": 100, "REF": "A", "ALT": "T",
         "SYMBOL": "TSC1", "DP": 10},
    ]).to_excel(p, index=False)

    seen = {}
    import run_background as rb
    real = rb.count_matrix
    rb.count_matrix = lambda v, b, **k: {s: {(c, pos): dict(rb._EMPTY, depth=10, T=5)
                                            for c, pos, _ in v} for s in b}
    rb.raw_bams = lambda d: {"A_S1": "x", "B_S2": "y"}
    try:
        annotate_workbooks([p], tmp_path, log=lambda *a: None)
    finally:
        rb.count_matrix = real

    ws = load_workbook(tmp_path / "run_2_all_variants_20260903_bg.xlsx").active
    dp = [c for c in ws[1] if c.value == "DP"][0].column
    assert "A_S1 (patient)" in ws.cell(2, dp).comment.text


# --- ratio AF primaire / AF brute, et code couleur --------------------------

def test_primary_ratio_collapses_on_structural_artefacts():
    """A real variant keeps its fraction when supplementary reads are dropped;
    a structural artefact loses it. Validated on the run's positive control."""
    from run_background import primary_ratio
    real = {"depth": 100, "T": 0, "G": 10, "A": 0, "C": 0, "ins": 0, "del": 0,
            "pri": {"depth": 98, "T": 0, "G": 10, "A": 0, "C": 0, "ins": 0, "del": 0}}
    fake = {"depth": 100, "T": 0, "G": 10, "A": 0, "C": 0, "ins": 0, "del": 0,
            "pri": {"depth": 70, "T": 0, "G": 1, "A": 0, "C": 0, "ins": 0, "del": 0}}
    assert primary_ratio(real, "T", "G") > 0.9
    assert primary_ratio(fake, "T", "G") < 0.3


def test_primary_ratio_is_none_without_primary_counts():
    from run_background import primary_ratio
    assert primary_ratio({"depth": 10, "G": 1}, "T", "G") is None


@pytest.mark.parametrize("contrast,ratio,expected", [
    (500.0, 1.00, "green"),    # causal: stands out, reads are primary
    (498.0, 0.08, "red"),      # G4 artefact: high contrast but reads vanish
    (1.2,   1.00, "red"),      # indistinguishable from the run's noise
    (3.1,   1.00, "yellow"),   # above the noise but the site is chronically busy
    (25.0,  0.95, "green"),    # TPOS mosaic
    # The run's positive control measures x7.3: the only mosaic known to be real
    # here must come out green, which sets contrast_good at 5 rather than 10.
    (7.3,   1.02, "green"),
    (5.0,   0.70, "orange"),   # unclear on both axes
    (None,  None, "grey"),     # not measured
])
def test_dp_colour_combines_both_detectors(contrast, ratio, expected):
    """Either detector alone can condemn: shared noise and structural artefact
    are independent failure modes."""
    from run_background import dp_colour
    assert dp_colour(contrast, ratio) == expected
