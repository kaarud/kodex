import pytest
from unittest.mock import patch


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


from pathlib import Path


def test_extract_sample_id():
    from igv_capture import extract_sample_id
    assert extract_sample_id(Path("TSC_063_345822_S2_filtered.xlsx")) == "345822_S2"
    assert extract_sample_id(Path("TSC_063_346769_S1_filtered.xlsx")) == "346769_S1"


def test_find_bams_all_present(tmp_path):
    from igv_capture import find_bams
    sid = "345822_S2"
    for name in [f"{sid}.bam", f"{sid}.bam.bai",
                 f"{sid}.mutect2.bam", f"{sid}.mutect2.bam.bai",
                 f"{sid}.chim.bam", f"{sid}.chim.bam.bai"]:
        (tmp_path / name).touch()
    bams = find_bams(sid, tmp_path)
    assert set(bams.keys()) == {"raw", "mutect2", "chim"}
    assert bams["raw"]["bam"].name == f"{sid}.bam"
    assert bams["raw"]["bai"].name == f"{sid}.bam.bai"


def test_find_bams_bai_fallback(tmp_path):
    """Accept {stem}.bai (primary) if {file}.bam.bai absent."""
    from igv_capture import find_bams
    sid = "345822_S2"
    (tmp_path / f"{sid}.bam").touch()
    (tmp_path / f"{sid}.bai").touch()          # .bai not .bam.bai — primary lookup
    (tmp_path / f"{sid}.mutect2.bam").touch()
    (tmp_path / f"{sid}.mutect2.bam.bai").touch()
    (tmp_path / f"{sid}.chim.bam").touch()
    (tmp_path / f"{sid}.chim.bam.bai").touch()
    bams = find_bams(sid, tmp_path)
    assert bams["raw"]["bai"].name == f"{sid}.bai"


def test_find_bams_missing_returns_none(tmp_path):
    from igv_capture import find_bams
    sid = "345822_S2"
    # Only raw BAM present, others missing
    (tmp_path / f"{sid}.bam").touch()
    (tmp_path / f"{sid}.bam.bai").touch()
    assert find_bams(sid, tmp_path) is None


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


from unittest.mock import patch, MagicMock


def _mock_get(status=200, text="OK"):
    m = MagicMock()
    m.status_code = status
    m.text = text
    return m


def test_check_igv_reachable():
    from igv_capture import check_igv
    with patch("requests.get", return_value=_mock_get(200)):
        assert check_igv(60151) is True


def test_check_igv_unreachable():
    from igv_capture import check_igv
    with patch("requests.get", side_effect=Exception("refused")):
        assert check_igv(60151) is False


def test_igv_cmd_builds_correct_url():
    from igv_capture import igv_cmd
    with patch("requests.get", return_value=_mock_get()) as mock_get:
        igv_cmd(60151, "goto", locus="chr9:100-200")
        mock_get.assert_called_once_with(
            "http://localhost:60151/goto",
            params={"locus": "chr9:100-200"},
            timeout=10,
        )


def test_igv_cmd_setPreference():
    from igv_capture import igv_cmd
    with patch("requests.get", return_value=_mock_get()) as mock_get:
        igv_cmd(60151, "setPreference", name="SAM.COLOR_BY", value="READ_STRAND")
        mock_get.assert_called_once_with(
            "http://localhost:60151/setPreference",
            params={"name": "SAM.COLOR_BY", "value": "READ_STRAND"},
            timeout=10,
        )


def test_sanitize_filename_special_chars():
    from igv_capture import sanitize_filename
    assert sanitize_filename("c.100C>T") == "c.100CgtT"
    assert sanitize_filename("c.1234+1G>A") == "c.1234+1GgtA"
    assert sanitize_filename("c.*145del") == "c.star145del"
    assert sanitize_filename("NM_000368.5:c.100C>T") == "NM_000368.5_c.100CgtT"
    assert sanitize_filename("c.100 del") == "c.100_del"


def test_sanitize_filename_empty_fallback():
    from igv_capture import build_variant_name
    v = {"SYMBOL": "TSC1", "HGVSc": None, "CHR": "chr9",
         "POS": 135786850, "REF": "C", "ALT": "T"}
    assert build_variant_name(v) == "TSC1_chr9_135786850_C_T"


def test_build_variant_name_normal():
    from igv_capture import build_variant_name
    v = {"SYMBOL": "TSC2", "HGVSc": "c.500A>G", "CHR": "chr16",
         "POS": 2000000, "REF": "A", "ALT": "G"}
    assert build_variant_name(v) == "TSC2_c.500AgtG"


def test_normalize_chr():
    from igv_capture import normalize_chr
    assert normalize_chr("9") == "chr9"
    assert normalize_chr("chr9") == "chr9"
    assert normalize_chr("chrX") == "chrX"
    assert normalize_chr("X") == "chrX"


def test_init_igv_sends_correct_commands(tmp_path):
    from igv_capture import init_igv
    sid = "345822_S2"
    bams = {
        "raw":    {"bam": tmp_path / f"{sid}.bam",    "bai": tmp_path / f"{sid}.bam.bai"},
        "mutect2":{"bam": tmp_path / f"{sid}.mutect2.bam", "bai": tmp_path / f"{sid}.mutect2.bam.bai"},
        "chim":   {"bam": tmp_path / f"{sid}.chim.bam",   "bai": tmp_path / f"{sid}.chim.bam.bai"},
    }
    calls = []
    with patch("igv_capture.igv_cmd", side_effect=lambda p, c, **kw: calls.append((c, kw))):
        init_igv(60151, bams)
    commands = [c for c, _ in calls]
    assert commands[0] == "new"
    assert commands[1] == "genome"
    assert commands.count("load") == 3


def test_main_no_annotated_variants(make_filtered_xlsx, capsys):
    from igv_capture import main
    path = make_filtered_xlsx({
        "HIGH_IMPACT": [
            {"SYMBOL": "TSC1", "HGVSc": "c.100C>T", "CHR": "chr9",
             "POS": 135786850, "REF": "C", "ALT": "T", "DP": 120,
             "USER_CLASS": None, "USER_ANNOT": None, "USER_COM": None},
        ]
    })
    with pytest.raises(SystemExit) as exc:
        main([str(path)])
    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "aucun variant annoté" in captured.out.lower()


def test_main_igv_unreachable(make_filtered_xlsx, capsys):
    from igv_capture import main
    path = make_filtered_xlsx({
        "HIGH_IMPACT": [
            {"SYMBOL": "TSC1", "HGVSc": "c.100C>T", "CHR": "chr9",
             "POS": 135786850, "REF": "C", "ALT": "T", "DP": 120,
             "USER_CLASS": "5", "USER_ANNOT": None, "USER_COM": None},
        ]
    })
    with patch("igv_capture.check_igv", return_value=False), \
         pytest.raises(SystemExit) as exc:
        main([str(path)])
    assert exc.value.code == 1


def test_main_bams_missing(make_filtered_xlsx, tmp_path, capsys):
    from igv_capture import main
    path = make_filtered_xlsx({
        "HIGH_IMPACT": [
            {"SYMBOL": "TSC1", "HGVSc": "c.100C>T", "CHR": "chr9",
             "POS": 135786850, "REF": "C", "ALT": "T", "DP": 120,
             "USER_CLASS": "5", "USER_ANNOT": None, "USER_COM": None},
        ]
    })
    # bam_dir exists but is empty — find_bams returns None
    with patch("igv_capture.check_igv", return_value=True), \
         pytest.raises(SystemExit) as exc:
        main([str(path), "--bam-dir", str(tmp_path)])
    assert exc.value.code == 1


def test_capture_variant_sequence(tmp_path):
    from igv_capture import capture_variant
    variant = {
        "SYMBOL": "TSC1", "HGVSc": "c.100C>T",
        "CHR": "chr9", "POS": 135786850, "REF": "C", "ALT": "T",
        "DP": 100,
    }
    calls = []
    with patch("igv_capture.igv_cmd", side_effect=lambda p, c, **kw: calls.append((c, kw))), \
         patch("time.sleep"):
        result = capture_variant(60151, variant, tmp_path / "out", delay=0)

    commands = [c for c, _ in calls]
    # Panel height set before goto
    assert "maxPanelHeight" in commands
    assert "goto" in commands
    # Strand capture before soft clips
    strand_idx = next(i for i,(c,kw) in enumerate(calls)
                      if c == "setPreference" and kw.get("value") == "READ_STRAND")
    softclip_idx = next(i for i,(c,kw) in enumerate(calls)
                        if c == "setPreference" and kw.get("name") == "SAM.SHOW_SOFT_CLIPPED"
                        and kw.get("value") == "true")
    assert strand_idx < softclip_idx
    # Exactly 2 snapshots taken
    assert commands.count("snapshot") == 2
    # Returns two PNG paths
    assert result["strand"].suffix == ".png"
    assert result["softclip"].suffix == ".png"
