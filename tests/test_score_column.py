"""SCORE column: computed on the DataFrame so the sheets can be ranked."""

import pandas as pd

from vcf_filter import add_score_column


def _row(**over):
    base = {
        "IMPACT": "MODIFIER", "Consequence": "intron_variant", "CLINVAR": None,
        "SB": 0.95, "CALLNB": 1, "RNCNT_TOTAL": 5, "PJCNT_TOTAL": 8,
        "CADD_PHRED": 2.0, "REVEL": None, "SpliceAI_num": None,
        "MAX_AF_GNOMADEX2.1": 0.004, "HGVSc": "c.1A>G",
    }
    base.update(over)
    return base


_GOOD = _row(**{"IMPACT": "HIGH", "Consequence": "stop_gained",
                "CLINVAR": "Pathogenic", "SB": 0.50, "CALLNB": 4,
                "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0, "CADD_PHRED": 35.0,
                "MAX_AF_GNOMADEX2.1": None, "HGVSc": "c.2T>A"})


def test_adds_score_between_zero_and_one():
    out = add_score_column(pd.DataFrame([_GOOD]))
    assert "SCORE" in out.columns
    assert 0.0 <= out["SCORE"].iloc[0] <= 1.0


def test_sorts_best_first():
    out = add_score_column(pd.DataFrame([_row(), _GOOD]))
    assert out["HGVSc"].tolist() == ["c.2T>A", "c.1A>G"]
    assert out["SCORE"].iloc[0] > out["SCORE"].iloc[1]


def test_ties_keep_original_order():
    rows = [_row(HGVSc="c.1A>G"), _row(HGVSc="c.2A>G"), _row(HGVSc="c.3A>G")]
    out = add_score_column(pd.DataFrame(rows))
    assert out["SCORE"].nunique() == 1
    assert out["HGVSc"].tolist() == ["c.1A>G", "c.2A>G", "c.3A>G"]


def test_noop_when_scoring_disabled():
    df = pd.DataFrame([_GOOD])
    out = add_score_column(df, {"enabled": False})
    assert "SCORE" not in out.columns


def test_noop_on_empty_frame():
    df = pd.DataFrame()
    assert add_score_column(df).empty
