"""The aggregation must score on the profile's gnomAD column, not the default.

config.yaml (069) declares af_filter_column: GNOMAD4_AF, config_legacy.yaml
(063) declares MAX_AF_GNOMADEX2.1. Both columns exist in the sheets, so using
the wrong one is silent: the aggregate scores variants differently from the
per-patient workbooks they came from.
"""

import pandas as pd
from aggregate_annotations import aggregate_all


def _variant(**kw):
    return {"CHR": "chr9", "POS": 100, "REF": "A", "ALT": "T", "SYMBOL": "TSC1",
            "IMPACT": "HIGH", "DP": 500, "AF_PCT": 45.0, "SB": 0.5, "CALLNB": 4,
            "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0, **kw}


def _one_file(tmp_path, rows):
    path = tmp_path / "run_2_Singularity_100007_S6_filtered.xlsx"
    with pd.ExcelWriter(path, engine="openpyxl") as w:
        pd.DataFrame({"Key": ["date"], "Value": ["2026-09-03"]}).to_excel(
            w, sheet_name="_summary", index=False)
        pd.DataFrame(rows).to_excel(w, sheet_name="HIGH_IMPACT", index=False)
    return path


def test_rarity_follows_the_configured_af_column(tmp_path):
    """A variant absent from gnomAD v4 but common in ex2.1 must rank as rare
    under the 069 profile — the two columns disagree on purpose here."""
    path = _one_file(tmp_path, [
        _variant(POS=100, GNOMAD4_AF=0.0, **{"MAX_AF_GNOMADEX2.1": 0.40}),
    ])
    v4 = aggregate_all(tmp_path, files=[path], every_variant=True,
                       af_column="GNOMAD4_AF")
    ex21 = aggregate_all(tmp_path, files=[path], every_variant=True,
                         af_column="MAX_AF_GNOMADEX2.1")
    assert v4.loc[0, "SCORE"] > ex21.loc[0, "SCORE"]
