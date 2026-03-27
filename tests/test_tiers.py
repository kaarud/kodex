from vcf_filter import classify_sheets


def _cfg():
    return {
        "sheets": {
            "high_impact": {"enabled": True},
            "clinvar_p_lp": {"enabled": True},
            "clinvar_conflicting": {"enabled": True},
            "clinvar_vus": {"enabled": True},
            "clinvar_not_provided": {"enabled": True},
            "fq_moderate": {"enabled": True, "gnomad_max_af": 0.001},
            "fq_low": {"enabled": True, "gnomad_max_af": 0.001},
            "fq_modifier": {"enabled": True, "gnomad_max_af": 0.001, "spliceai_threshold": 0.50},
            "four_callers": {"enabled": True, "min_callers": 4},
            "lovd": {"enabled": True},
            "mosaic": {"enabled": True, "af_min_pct": 5.0, "af_max_pct": 35.0, "rncnt_max_total": 2},
        }
    }


def test_sheet_high_impact(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"IMPACT": "HIGH", "Consequence": "stop_gained"},
        {"IMPACT": "MODERATE", "Consequence": "missense_variant"},
        {"IMPACT": "HIGH", "Consequence": "frameshift_variant"},
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["HIGH_IMPACT"]) == 2


def test_sheet_clinvar_p_lp(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"CLINVAR": "Pathogenic", "IMPACT": "MODIFIER"},
        {"CLINVAR": "Likely_pathogenic", "IMPACT": "MODIFIER"},
        {"CLINVAR": "Conflicting_interpretations_of_pathogenicity", "IMPACT": "MODIFIER"},
        {"CLINVAR": None, "IMPACT": "MODIFIER"},
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["CLINVAR_P_LP"]) == 2  # excludes Conflicting


def test_sheet_clinvar_conflicting(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"CLINVAR": "Conflicting_interpretations_of_pathogenicity", "IMPACT": "MODIFIER"},
        {"CLINVAR": "Pathogenic", "IMPACT": "MODIFIER"},
        {"CLINVAR": None, "IMPACT": "MODIFIER"},
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["CLINVAR_CONFLICTING"]) == 1


def test_sheet_clinvar_vus(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"CLINVAR": "Uncertain_significance", "IMPACT": "MODIFIER"},
        {"CLINVAR": "Pathogenic", "IMPACT": "MODIFIER"},
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["CLINVAR_VUS"]) == 1


def test_sheet_clinvar_not_provided(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"CLINVAR": "not_provided", "IMPACT": "MODIFIER"},
        {"CLINVAR": "Pathogenic", "IMPACT": "MODIFIER"},
        {"CLINVAR": None, "IMPACT": "MODIFIER"},
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["CLINVAR_NOT_PROVIDED"]) == 1


def test_sheet_fq_moderate(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"IMPACT": "MODERATE", "MAX_AF_GNOMADEX2.1": 0.0005},  # pass
        {"IMPACT": "MODERATE", "MAX_AF_GNOMADEX2.1": 0.002},   # gnomAD too high
        {"IMPACT": "LOW", "MAX_AF_GNOMADEX2.1": 0.0005},       # wrong impact
        {"IMPACT": "MODERATE", "MAX_AF_GNOMADEX2.1": None},    # novel -> pass (0 < 0.001)
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["FQ_MODERATE"]) == 2  # rows 0, 3


def test_sheet_fq_low(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"IMPACT": "LOW", "MAX_AF_GNOMADEX2.1": 0.0005},
        {"IMPACT": "LOW", "MAX_AF_GNOMADEX2.1": 0.002},    # gnomAD too high
        {"IMPACT": "MODERATE", "MAX_AF_GNOMADEX2.1": 0.0005},  # wrong impact
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["FQ_LOW"]) == 1


def test_sheet_fq_modifier_spliceai(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"IMPACT": "MODIFIER", "MAX_AF_GNOMADEX2.1": 0.0005, "SpliceAI_num": 0.60},  # pass via spliceAI
        {"IMPACT": "MODIFIER", "MAX_AF_GNOMADEX2.1": 0.0005, "SPiP_Interpretation": "Altered"},  # pass via SPiP
        {"IMPACT": "MODIFIER", "MAX_AF_GNOMADEX2.1": 0.0005, "SpliceAI_num": 0.10},  # splice too low
        {"IMPACT": "MODIFIER", "MAX_AF_GNOMADEX2.1": 0.002, "SpliceAI_num": 0.60},   # gnomAD too high
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["FQ_MODIFIER"]) == 2  # rows 0, 1


def test_sheet_four_callers(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"CALLNB": 4},
        {"CALLNB": 3},
        {"CALLNB": 4},
        {"CALLNB": None},
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["4_CALLERS"]) == 2


def test_sheet_lovd(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"LOVD": "some_entry"},
        {"LOVD": "."},
        {"LOVD": None},
        {"LOVD": "another_entry"},
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["LOVD"]) == 2


def test_sheet_mosaic(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"AF_PCT": 15.0, "RNCNT_TOTAL": 1},   # pass
        {"AF_PCT": 48.0, "RNCNT_TOTAL": 0},    # AF > 35% -> excluded (germinal HTZ)
        {"AF_PCT": 3.0, "RNCNT_TOTAL": 0},     # AF < 5% -> excluded
        {"AF_PCT": 15.0, "RNCNT_TOTAL": 5},    # RNCNT > 2 -> excluded
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    assert len(sheets["MOSAIC"]) == 1  # row 0 only


def test_overlap_between_sheets(make_xlsx):
    """A variant can appear in multiple sheets."""
    from vcf_filter import load_variants
    path = make_xlsx([
        {
            "IMPACT": "HIGH",
            "CLINVAR": "Pathogenic",
            "CALLNB": 4,
            "LOVD": "entry",
            "AF_PCT": 15.0,
            "RNCNT_TOTAL": 1,
        },
    ])
    df = load_variants(path)
    sheets = classify_sheets(df, _cfg())
    # Should appear in HIGH_IMPACT, CLINVAR_P_LP, 4_CALLERS, LOVD (AF=15% < 35% → MOSAIC too)
    assert len(sheets["HIGH_IMPACT"]) == 1
    assert len(sheets["CLINVAR_P_LP"]) == 1
    assert len(sheets["4_CALLERS"]) == 1
    assert len(sheets["LOVD"]) == 1
    assert len(sheets["MOSAIC"]) == 1
