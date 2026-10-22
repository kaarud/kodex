import pytest

from vcf_filter import SchemaMismatchError, load_variants

ALIASES_069 = {"GENE": "SYMBOL", "CSQ": "Consequence", "REGION": "REGIONS"}


def test_load_069_applies_aliases(make_xlsx_069):
    path = make_xlsx_069([{"SYMBOL": "TSC2"}])  # SYMBOL→stored as GENE by factory
    df = load_variants(
        path,
        genes=["TSC1", "TSC2"],
        column_aliases=ALIASES_069,
        synthesize_user_columns=True,
        required_columns=["SYMBOL", "GNOMAD4_AF"],
    )
    assert "SYMBOL" in df.columns
    assert "Consequence" in df.columns
    assert "REGIONS" in df.columns
    assert "GENE" not in df.columns
    assert df["SYMBOL"].iloc[0] == "TSC2"


def test_load_069_synthesizes_empty_user_columns(make_xlsx_069):
    path = make_xlsx_069()
    df = load_variants(
        path,
        genes=["TSC1", "TSC2"],
        column_aliases=ALIASES_069,
        synthesize_user_columns=True,
        required_columns=["SYMBOL", "GNOMAD4_AF"],
    )
    for col in ("USER_CLASS", "USER_ANNOT", "USER_COM"):
        assert col in df.columns
        assert df[col].isna().all()


def test_load_069_missing_required_column_raises(make_xlsx_069):
    # No aliases passed → GENE never renamed to SYMBOL → mismatch
    path = make_xlsx_069()
    with pytest.raises(SchemaMismatchError) as exc:
        load_variants(
            path,
            genes=["TSC1", "TSC2"],
            column_aliases=None,
            synthesize_user_columns=False,
            required_columns=["SYMBOL"],
        )
    assert "SYMBOL" in str(exc.value)
    assert "--legacy" in str(exc.value)


def test_load_063_unchanged_under_legacy_defaults(make_xlsx):
    # Legacy path: no aliases, no synthesis — behaves exactly as before.
    path = make_xlsx([{"SYMBOL": "TSC1"}, {"SYMBOL": "TSC2"}])
    df = load_variants(
        path,
        genes=["TSC1", "TSC2"],
        column_aliases={},
        synthesize_user_columns=False,
        required_columns=["SYMBOL", "MAX_AF_GNOMADEX2.1"],
    )
    assert set(df["SYMBOL"].unique()) == {"TSC1", "TSC2"}
    assert "USER_CLASS" in df.columns  # 063 already has it


def test_load_069_blank_pjcnt_normalised_to_zero(make_xlsx_069):
    # Singularity writes a project count of 0 as an empty cell. After load it
    # must be an integer 0 (not NaN), so display/highlight/scoring match 063.
    path = make_xlsx_069([
        {"SYMBOL": "TSC2", "RNCNT_TOTAL": 6, "RNCNT_LOW": 6,
         "PJCNT_TOTAL": None, "PJCNT_LOW": None, "PJCNT_HTZ": None},
    ])
    df = load_variants(
        path,
        genes=["TSC1", "TSC2"],
        column_aliases=ALIASES_069,
        synthesize_user_columns=True,
        required_columns=["SYMBOL", "GNOMAD4_AF"],
    )
    for col in ("PJCNT_TOTAL", "PJCNT_LOW", "PJCNT_HTZ"):
        assert df[col].notna().all(), f"{col} still has NaN"
        assert df[col].iloc[0] == 0
        assert df[col].dtype == "int64"
    # Run-level counts pass through untouched
    assert df["RNCNT_TOTAL"].iloc[0] == 6


from vcf_filter import apply_exclusions, classify_sheets


def _base_cfg():
    return {
        "exclusions": {
            "gnomad_max_af": 0.005,
            "rncnt_max_total": 6,
            "pjcnt_max_total": 10,
            "exclude_hom": True,
            "exclude_clinvar_benign": True,
            "mosaic_artifact": {"enabled": False},
        },
        "af_filter_column": "GNOMAD4_AF",
    }


def test_exclusions_use_af_filter_column(make_xlsx_069):
    import pandas as pd
    # Variant common in GNOMAD4_AF but rare in popmax → must be excluded
    df = pd.DataFrame([
        {"GNOMAD4_AF": 0.02, "MAX_AF_GNOMADEX2.1": 0.0001,
         "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0, "RNCNT_HOM": 0,
         "RNCNT_LOW": 0, "AF_PCT": 48.0, "CLINVAR": None},
        {"GNOMAD4_AF": 0.0001, "MAX_AF_GNOMADEX2.1": 0.0001,
         "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0, "RNCNT_HOM": 0,
         "RNCNT_LOW": 0, "AF_PCT": 48.0, "CLINVAR": None},
    ])
    out, report = apply_exclusions(df, _base_cfg())
    assert len(out) == 1
    assert report["gnomad"] == 1


def test_classify_sheets_use_af_filter_column():
    import pandas as pd
    cfg = {
        "af_filter_column": "GNOMAD4_AF",
        "sheets": {"fq_moderate": {"enabled": True, "gnomad_max_af": 0.001}},
    }
    df = pd.DataFrame([
        # rare in GNOMAD4_AF + MODERATE → in FQ_MODERATE
        {"GNOMAD4_AF": 0.0001, "MAX_AF_GNOMADEX2.1": 0.5, "IMPACT": "MODERATE",
         "CLINVAR": None},
        # common in GNOMAD4_AF → excluded from FQ_MODERATE
        {"GNOMAD4_AF": 0.5, "MAX_AF_GNOMADEX2.1": 0.0001, "IMPACT": "MODERATE",
         "CLINVAR": None},
    ])
    sheets = classify_sheets(df, cfg)
    assert len(sheets["FQ_MODERATE"]) == 1


from vcf_filter import _cell_fill, _score_variant


def test_cell_fill_honours_af_column():
    # With af_column=GNOMAD4_AF, that column gets the AF coloring...
    assert _cell_fill("GNOMAD4_AF", None, af_column="GNOMAD4_AF") is not None
    # ...and the default popmax column no longer matches.
    assert _cell_fill("MAX_AF_GNOMADEX2.1", 0.0, af_column="GNOMAD4_AF") is None


def test_score_variant_reads_af_column():
    row = {
        "SB": 0.50, "CALLNB": 4, "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0,
        "DP": 200, "AF_PCT": 48.0, "CLINVAR": "Pathogenic", "IMPACT": "HIGH",
        "Consequence": "stop_gained", "CADD_PHRED": 35.0, "REVEL": 0.9,
        "GNOMAD4_AF": None, "MAX_AF_GNOMADEX2.1": 0.5,  # popmax common, v4 novel
    }
    *_, comment = _score_variant(row, af_column="GNOMAD4_AF")
    assert "Absent gnomAD (novel)" in comment


import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent


def _run(args):
    return subprocess.run(
        [sys.executable, "vcf_filter.py", *args],
        capture_output=True, text=True, cwd=str(REPO),
    )


def test_cli_069_default_config_produces_output(make_xlsx_069):
    rows = [
        {"IMPACT": "HIGH", "Consequence": "stop_gained",
         "GNOMAD4_AF": 0.0, "MAX_AF_GNOMADEX2.1": 0.0,
         "CHR": "chr9", "POS": 200, "REF": "C", "ALT": "T"},
    ]
    path = make_xlsx_069(rows)
    result = _run([str(path)])  # no --config, no --legacy → 069 default
    assert result.returncode == 0, f"STDERR: {result.stderr}"
    out = path.with_name(path.stem + "_filtered.xlsx")
    assert out.exists()
    sheets = pd.ExcelFile(out, engine="openpyxl").sheet_names
    assert "HIGH_IMPACT" in sheets


def test_cli_063_legacy_flag_still_works(make_xlsx):
    rows = [
        {"IMPACT": "HIGH", "Consequence": "stop_gained",
         "CHR": "chr9", "POS": 200, "REF": "C", "ALT": "T"},
    ]
    path = make_xlsx(rows)
    result = _run(["--legacy", str(path)])
    assert result.returncode == 0, f"STDERR: {result.stderr}"
    out = path.with_name(path.stem + "_filtered.xlsx")
    assert out.exists()


def test_cli_063_under_default_config_exits_clean(make_xlsx):
    # 063 file lacks GNOMAD4_AF → default (069) config must fail cleanly.
    path = make_xlsx([{"SYMBOL": "TSC1"}])
    result = _run([str(path)])
    assert result.returncode != 0
    assert "GNOMAD4_AF" in (result.stdout + result.stderr)
    assert "--legacy" in (result.stdout + result.stderr)
    assert "Traceback" not in result.stderr


def test_cli_069_under_legacy_exits_clean(make_xlsx_069):
    # 069 file under legacy config: GENE never aliased → SYMBOL missing.
    path = make_xlsx_069()
    result = _run(["--legacy", str(path)])
    assert result.returncode != 0
    assert "SYMBOL" in (result.stdout + result.stderr)
    assert "Traceback" not in result.stderr
