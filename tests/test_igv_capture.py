def test_load_annotated_variants_single_sheet(make_filtered_xlsx):
    from igv_capture import load_annotated_variants
    path = make_filtered_xlsx({
        "HIGH_IMPACT": [
            {"SYMBOL": "TSC1", "HGVSc": "c.100C>T", "CHR": "chr9",
             "POS": 135786850, "REF": "C", "ALT": "T", "DP": 120,
             "USER_CLASS": "5", "USER_ANNOT": None, "USER_COM": None},
            {"SYMBOL": "TSC1", "HGVSc": "c.200G>A", "CHR": "chr9",
             "POS": 135786900, "REF": "G", "ALT": "A", "DP": 80,
             "USER_CLASS": None, "USER_ANNOT": None, "USER_COM": None},
        ]
    })
    variants = load_annotated_variants(path)
    assert len(variants) == 1
    assert variants[0]["HGVSc"] == "c.100C>T"


def test_load_annotated_variants_dedup_across_sheets(make_filtered_xlsx):
    from igv_capture import load_annotated_variants
    same_row = {"SYMBOL": "TSC2", "HGVSc": "c.500A>G", "CHR": "chr16",
                "POS": 2000000, "REF": "A", "ALT": "G", "DP": 200,
                "USER_CLASS": None, "USER_ANNOT": "check", "USER_COM": None}
    path = make_filtered_xlsx({
        "HIGH_IMPACT": [same_row],
        "LOVD":        [same_row],
    })
    variants = load_annotated_variants(path)
    assert len(variants) == 1  # deduplicated


def test_load_annotated_variants_user_com_triggers(make_filtered_xlsx):
    from igv_capture import load_annotated_variants
    path = make_filtered_xlsx({
        "MOSAIC": [
            {"SYMBOL": "TSC1", "HGVSc": "c.300T>C", "CHR": "chr9",
             "POS": 135787000, "REF": "T", "ALT": "C", "DP": 50,
             "USER_CLASS": None, "USER_ANNOT": None, "USER_COM": "à revoir"},
        ]
    })
    variants = load_annotated_variants(path)
    assert len(variants) == 1
