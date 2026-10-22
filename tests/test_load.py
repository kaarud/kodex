from vcf_filter import load_variants


def test_load_strips_sample_prefix(make_xlsx):
    path = make_xlsx()
    df = load_variants(path)
    assert "AF_PCT" in df.columns
    assert "DP" in df.columns
    assert "CALLFILTER_mutect2" in df.columns
    assert "SYMBOL" in df.columns
    # No column should still contain the sample prefix
    for col in df.columns:
        assert "\n" not in col, f"Column still has newline prefix: {col!r}"


def test_load_dedup_ref_trans(make_xlsx):
    """Positions with a VRAI transcript → keep only VRAI; others untouched."""
    path = make_xlsx([
        # Same position: 1 VRAI + 3 FAUX → keep only VRAI
        {"CHR": "chr16", "POS": 2083692, "REF": "CA", "ALT": "AT", "REF_TRANS": "VRAI"},
        {"CHR": "chr16", "POS": 2083692, "REF": "CA", "ALT": "AT", "REF_TRANS": "FAUX"},
        {"CHR": "chr16", "POS": 2083692, "REF": "CA", "ALT": "AT", "REF_TRANS": "FAUX"},
        # Different position, no VRAI → keep all
        {"CHR": "chr9", "POS": 135786850, "REF": "C", "ALT": "T", "REF_TRANS": "FAUX"},
        {"CHR": "chr9", "POS": 135786850, "REF": "C", "ALT": "T", "REF_TRANS": "FAUX"},
    ])
    df = load_variants(path)
    assert len(df) == 3  # 1 VRAI + 2 FAUX (different position)
    vrai_rows = df[df["REF_TRANS"].astype(str).str.upper() == "VRAI"]
    assert len(vrai_rows) == 1
    assert vrai_rows.iloc[0]["POS"] == 2083692


def test_load_dedup_ref_trans_yes_no(make_xlsx):
    """069 schema uses Yes/No: keep only Yes when present; else keep all."""
    path = make_xlsx([
        # Same position: 1 Yes + 3 No → keep only Yes (the RefSeq canonical)
        {"CHR": "chr9", "POS": 132896301, "REF": "CG", "ALT": "C", "REF_TRANS": "Yes"},
        {"CHR": "chr9", "POS": 132896301, "REF": "CG", "ALT": "C", "REF_TRANS": "No"},
        {"CHR": "chr9", "POS": 132896301, "REF": "CG", "ALT": "C", "REF_TRANS": "No"},
        {"CHR": "chr9", "POS": 132896301, "REF": "CG", "ALT": "C", "REF_TRANS": "No"},
        # Different position, no Yes → keep all
        {"CHR": "chr16", "POS": 2083692, "REF": "C", "ALT": "T", "REF_TRANS": "No"},
        {"CHR": "chr16", "POS": 2083692, "REF": "C", "ALT": "T", "REF_TRANS": "No"},
    ])
    df = load_variants(path)
    assert len(df) == 3  # 1 Yes + 2 No (different position)
    yes_rows = df[df["REF_TRANS"].astype(str).str.upper() == "YES"]
    assert len(yes_rows) == 1
    assert yes_rows.iloc[0]["POS"] == 132896301


def test_load_filters_genes(make_xlsx):
    path = make_xlsx([
        {"SYMBOL": "TSC1"},
        {"SYMBOL": "TSC2"},
        {"SYMBOL": "BRCA1"},
    ])
    df = load_variants(path, genes=["TSC1", "TSC2"])
    assert set(df["SYMBOL"].unique()) == {"TSC1", "TSC2"}
    assert len(df) == 2
