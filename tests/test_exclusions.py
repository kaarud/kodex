from vcf_filter import apply_exclusions


def _cfg():
    """Minimal config dict matching config.yaml structure."""
    return {
        "exclusions": {
            "gnomad_max_af": 0.005,
            "rncnt_max_total": 4,
            "pjcnt_max_total": 10,
            "exclude_hom": True,
            "exclude_clinvar_benign": True,
            "mosaic_artifact": {
                "enabled": True,
                "rncnt_low_max": 5,
                "af_max_pct": 35.0,
            },
        }
    }


def test_exclude_gnomad_too_frequent(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"MAX_AF_GNOMADEX2.1": 0.006},  # > 0.005 -> excluded
        {"MAX_AF_GNOMADEX2.1": 0.004},  # kept
        {"MAX_AF_GNOMADEX2.1": None},   # novel -> kept
    ])
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2
    assert report["gnomad"] == 1


def test_exclude_rncnt_total(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"RNCNT_TOTAL": 5},  # > 4 -> excluded (run actuel)
        {"RNCNT_TOTAL": 4},  # kept (= seuil, pas exclu)
        {"RNCNT_TOTAL": 0},  # kept
    ])
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2
    assert report["rncnt_total"] == 1


def test_exclude_pjcnt_total(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"PJCNT_TOTAL": 11},  # > 10 -> excluded (méta-analyse 30 runs)
        {"PJCNT_TOTAL": 10},  # kept (= seuil, pas exclu)
    ])
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 1
    assert report["pjcnt_total"] == 1


def test_exclude_hom(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"RNCNT_HOM": 1},   # > 0 -> excluded
        {"RNCNT_HOM": 0},   # kept
        {"RNCNT_HOM": None}, # treated as 0 -> kept
    ])
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2
    assert report["hom"] == 1


def test_exclude_clinvar_benign(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"CLINVAR": "Benign"},
        {"CLINVAR": "Likely_benign"},
        {"CLINVAR": "Benign/Likely_benign"},
        {"CLINVAR": "Pathogenic"},
        {"CLINVAR": None},
    ])
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2  # Pathogenic + None
    assert report["clinvar_benign"] == 3


def test_exclude_mosaic_artifact(make_xlsx):
    """Recurrent low-AF variant = sequencing hotspot."""
    from vcf_filter import load_variants
    path = make_xlsx([
        {"RNCNT_LOW": 6, "AF_PCT": 12.0},  # > 5 & < 35% -> excluded
        {"RNCNT_LOW": 6, "AF_PCT": 48.0},  # > 5 but AF > 35% -> kept
        {"RNCNT_LOW": 3, "AF_PCT": 12.0},  # <= 5 -> kept
    ])
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2
    assert report["mosaic_artifact"] == 1
