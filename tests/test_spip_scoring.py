"""SPiP contributes to the biological axis, cumulatively with SpliceAI.

The score ranks candidates to open in IGV; it is not an ACMG tally. Two methods
agreeing on a splice mechanism raise the odds a variant is worth looking at,
even though ClinGen would grant a single PP3 for both.
"""

import pytest
from vcf_filter import _score_variant, load_config
from pathlib import Path


@pytest.fixture
def scoring():
    return load_config(Path(__file__).resolve().parent.parent / "config.yaml")["scoring"]


def _row(**kw):
    return {"CHR": "chr9", "POS": 100, "REF": "A", "ALT": "T",
            "IMPACT": "MODIFIER", "DP": 500, "AF_PCT": 45.0, "SB": 0.5,
            "CALLNB": 4, "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0, **kw}


def _bio(row, scoring):
    """Biological axis only — index 0 of (bio, tech, rarity, combined, comment)."""
    return _score_variant(row, scoring)[0]


def test_spip_alter_alone_scores_above_ntr(scoring):
    """A variant recruited by SPiP alone must not score zero on the splice axis."""
    alter = _bio(_row(SPiP_Interpretation="Alter_by_create_New_splice_site"), scoring)
    ntr = _bio(_row(SPiP_Interpretation="NTR"), scoring)
    assert alter > ntr


def test_spip_and_spliceai_accumulate(scoring):
    """Both positive beats either alone — convergence is worth more for triage."""
    both = _bio(_row(SpliceAI_num=0.95, SPiP_Interpretation="Alter_by_SPiCE"), scoring)
    spliceai_only = _bio(_row(SpliceAI_num=0.95, SPiP_Interpretation="NTR"), scoring)
    spip_only = _bio(_row(SpliceAI_num=0.0, SPiP_Interpretation="Alter_by_SPiCE"), scoring)
    assert both > spliceai_only
    assert both > spip_only


def test_missing_spip_does_not_penalise(scoring):
    """No SPiP column: the axis is scored without it, not scored as a failure."""
    absent = _bio(_row(SpliceAI_num=0.95), scoring)
    ntr = _bio(_row(SpliceAI_num=0.95, SPiP_Interpretation="NTR"), scoring)
    assert absent > ntr


# --- CADD is uninformative on deep intronic variants ------------------------

def test_cadd_is_skipped_for_modifier_impact(scoring):
    """A MODIFIER variant must not be dragged down by a CADD it cannot score.

    Deep intronic variants sit at CADD 3-6 whatever their splice impact, so
    counting CADD against them dilutes the very axis FQ_MODIFIER exists to
    surface. Skipping it leaves bio_max untouched, exactly like an absent value.
    """
    low_cadd = _bio(_row(IMPACT="MODIFIER", CADD_PHRED=3.65), scoring)
    no_cadd = _bio(_row(IMPACT="MODIFIER"), scoring)
    assert low_cadd == no_cadd


def test_cadd_still_counts_for_coding_impact(scoring):
    """The skip is scoped: CADD keeps discriminating where it is calibrated."""
    strong = _bio(_row(IMPACT="MODERATE", CADD_PHRED=35.0), scoring)
    weak = _bio(_row(IMPACT="MODERATE", CADD_PHRED=3.65), scoring)
    assert strong > weak


# ---------------------------------------------------------------------------
# skip_when_consequence — un predicteur d'epissage ne juge pas un frameshift
# ---------------------------------------------------------------------------


def test_splice_predictors_skipped_on_frameshift(scoring):
    """SpliceAI bas + SPiP NTR ne doivent pas penaliser un decalage du cadre.

    Le frameshift porte deja sa perte de fonction ; le predicteur d'epissage
    mesure autre chose. Meme principe que cadd.skip_when_impact — le critere
    est non evalue, pas rate, donc il ne rentre pas au denominateur.
    """
    common = dict(IMPACT="HIGH", CLINVAR="not_provided",
                  SpliceAI_num=0.01, SPiP_Interpretation="NTR")
    frameshift = _bio(_row(Consequence="frameshift_variant", **common), scoring)
    missense = _bio(_row(Consequence="missense_variant", **common), scoring)
    assert frameshift > missense


def test_splice_predictors_skipped_on_stop_gained(scoring):
    row = _row(Consequence="stop_gained", IMPACT="HIGH", CLINVAR="not_provided",
               SpliceAI_num=0.01, SPiP_Interpretation="NTR")
    comment = _score_variant(row, scoring)[4]
    assert "SpliceAI" not in comment
    assert "SPiP" not in comment


def test_compound_consequence_matches_each_term(scoring):
    """VEP concatene les consequences : 'frameshift_variant&splice_region_variant'."""
    row = _row(Consequence="frameshift_variant&splice_region_variant",
               IMPACT="HIGH", SpliceAI_num=0.01, SPiP_Interpretation="NTR")
    assert "SpliceAI" not in _score_variant(row, scoring)[4]


def test_splice_predictors_still_count_on_splice_variants(scoring):
    """La regle ne doit pas neutraliser l'axe la ou il sert."""
    row = _row(Consequence="intron_variant", SpliceAI_num=0.95,
               SPiP_Interpretation="Alter_ESR")
    assert "SpliceAI=0.95" in _score_variant(row, scoring)[4]


def test_nan_annotation_counts_as_absent_not_as_failure(scoring):
    """pandas rend une cellule vide en NaN et float(nan) reussit.

    Sans garde, un CADD non annote entrait au denominateur de l'axe bio comme
    critere rate. Le chemin openpyxl (cellule vide = None) n'avait pas ce bug,
    d'ou une divergence entre la couleur d'une feuille et le SCORE de l'agregat.
    """
    base = _row(Consequence="missense_variant", IMPACT="MODERATE")
    without = _bio(base, scoring)
    with_nan = _bio({**base, "CADD_PHRED": float("nan"), "REVEL": float("nan")}, scoring)
    assert with_nan == without
