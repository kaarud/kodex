import pytest
from unittest.mock import patch


def _touch_run(tmp_path):
    # 100008_S2 is the patient of the make_filtered_xlsx fixture
    # (run_1_100008_S2_filtered.xlsx → extract_sample_id → 100008_S2),
    # so its BAMs must exist for the xlsx-row test to reach a successful goto.
    for sid in ("100008_S1", "100008_S2", "999111_S2"):
        for name in (f"{sid}.bam", f"{sid}.bam.bai",
                     f"{sid}.chim.bam", f"{sid}.chim.bam.bai",
                     f"{sid}.mutect2.bam", f"{sid}.mutect2.bam.bai"):
            (tmp_path / name).touch()
    return tmp_path


def test_main_sends_load_and_goto(tmp_path):
    from igv_goto import main
    bam_dir = _touch_run(tmp_path)
    calls = []
    with patch("igv_goto.check_igv", return_value=True), \
         patch("igv_goto.igv_cmd", side_effect=lambda p, c, **kw: calls.append((c, kw))):
        main(["100008_S1", "chr9:1000", "--bam-dir", str(bam_dir),
              "--reference-lane", "S2", "--half-window", "145"])
    commands = [c for c, _ in calls]
    assert commands[0] == "new"
    assert commands[1] == "genome"
    assert commands.count("load") == 4           # RAW, CHIM, MUTECT2, REF_S2
    assert commands[-1] == "goto"
    goto_kw = calls[-1][1]
    assert goto_kw["locus"] == "chr9:855-1145"


def test_main_igv_unreachable_exits_1(tmp_path):
    from igv_goto import main
    bam_dir = _touch_run(tmp_path)
    with patch("igv_goto.check_igv", return_value=False), \
         pytest.raises(SystemExit) as exc:
        main(["100008_S1", "chr9:1000", "--bam-dir", str(bam_dir)])
    assert exc.value.code == 1


def test_main_missing_bams_exits_1(tmp_path):
    from igv_goto import main
    with patch("igv_goto.check_igv", return_value=True), \
         pytest.raises(SystemExit) as exc:
        main(["100008_S1", "chr9:1000", "--bam-dir", str(tmp_path)])  # empty dir
    assert exc.value.code == 1


def test_main_xlsx_row_mode(tmp_path, make_filtered_xlsx):
    from igv_goto import main
    bam_dir = _touch_run(tmp_path)
    xlsx = make_filtered_xlsx({
        "HIGH_IMPACT": [
            {"SYMBOL": "TSC1", "CHR": "chr9", "POS": 1000, "REF": "C", "ALT": "T"},
            {"SYMBOL": "TSC2", "CHR": "chr16", "POS": 2098100, "REF": "A", "ALT": "G"},
        ]
    })
    calls = []
    with patch("igv_goto.check_igv", return_value=True), \
         patch("igv_goto.igv_cmd", side_effect=lambda p, c, **kw: calls.append((c, kw))):
        main([str(xlsx), "--sheet", "HIGH_IMPACT", "--row", "1",
              "--bam-dir", str(bam_dir), "--reference-lane", "S2"])
    goto_kw = calls[-1][1]
    # row index 1 (0-based, second data row) → chr16:2098100 ± 145
    assert goto_kw["locus"] == "chr16:2097955-2098245"
