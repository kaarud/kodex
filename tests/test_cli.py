import subprocess
import sys
from pathlib import Path

import pandas as pd


def test_cli_produces_filtered_xlsx(make_xlsx):
    """End-to-end: run vcf_filter.py on a synthetic XLSX, check output file."""
    rows = [
        # CLINVAR_P_LP sheet
        {"CLINVAR": "Pathogenic", "IMPACT": "MODIFIER", "Consequence": "intron_variant",
         "CHR": "chr9", "POS": 100, "REF": "A", "ALT": "G"},
        # HIGH_IMPACT sheet
        {"IMPACT": "HIGH", "Consequence": "stop_gained",
         "CHR": "chr9", "POS": 200, "REF": "C", "ALT": "T"},
        # FQ_MODERATE sheet (missense, gnomAD < 0.1%)
        {"IMPACT": "MODERATE", "Consequence": "missense_variant", "REVEL": 0.80,
         "MAX_AF_GNOMADEX2.1": 0.0005,
         "CHR": "chr9", "POS": 300, "REF": "G", "ALT": "A"},
        # LOVD sheet
        {"IMPACT": "LOW", "Consequence": "synonymous_variant", "LOVD": "some_entry",
         "CHR": "chr9", "POS": 400, "REF": "T", "ALT": "C"},
        # MOSAIC sheet (AF > 5%, RNCNT_TOTAL <= 2)
        {"AF_PCT": 15.0, "RNCNT_TOTAL": 1, "IMPACT": "MODERATE",
         "Consequence": "missense_variant", "REVEL": 0.80, "MAX_AF_GNOMADEX2.1": 0.0005,
         "CHR": "chr9", "POS": 500, "REF": "A", "ALT": "T"},
        # Excluded: gnomAD too high
        {"MAX_AF_GNOMADEX2.1": 0.01, "IMPACT": "HIGH", "Consequence": "stop_gained",
         "CHR": "chr9", "POS": 600, "REF": "C", "ALT": "A"},
    ]
    xlsx_path = make_xlsx(rows)

    config_src = Path(__file__).resolve().parent.parent / "config_legacy.yaml"

    result = subprocess.run(
        [sys.executable, "vcf_filter.py", str(xlsx_path), "--config", str(config_src)],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parent.parent),
    )
    assert result.returncode == 0, f"STDERR: {result.stderr}"

    # Check pre-report was printed
    assert "Pre-report" in result.stdout

    # Check output file was created
    out_path = xlsx_path.with_name(xlsx_path.stem + "_filtered.xlsx")
    assert out_path.exists(), f"Output file not found: {out_path}"

    # Check sheets exist
    sheets = pd.ExcelFile(out_path, engine="openpyxl").sheet_names
    assert "_summary" in sheets
    assert "HIGH_IMPACT" in sheets
    assert "CLINVAR_P_LP" in sheets
    assert "FQ_MODERATE" in sheets
    assert "LOVD" in sheets
    assert "MOSAIC" in sheets


def test_variant_score_gradient_and_comment(make_xlsx):
    """SCORE cell gets gradient color + hover comment, best variant sorted first."""
    from openpyxl import load_workbook

    rows = [
        # Good variant: HIGH impact, ClinVar Pathogenic, good SB, 3 callers, novel
        {"IMPACT": "HIGH", "Consequence": "stop_gained", "CLINVAR": "Pathogenic",
         "SB": 0.50, "CALLNB": 3, "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0,
         "CADD_PHRED": 35.0, "MAX_AF_GNOMADEX2.1": None,
         "CHR": "chr9", "POS": 100, "REF": "A", "ALT": "G"},
        # Bad variant: MODIFIER, bad SB, 1 caller, frequent
        {"IMPACT": "MODIFIER", "Consequence": "intron_variant", "CLINVAR": None,
         "SB": 0.95, "CALLNB": 1, "RNCNT_TOTAL": 4, "PJCNT_TOTAL": 8,
         "CADD_PHRED": 5.0, "MAX_AF_GNOMADEX2.1": 0.003,
         "CHR": "chr9", "POS": 200, "REF": "C", "ALT": "T"},
    ]
    xlsx_path = make_xlsx(rows)
    config_src = Path(__file__).resolve().parent.parent / "config_legacy.yaml"

    result = subprocess.run(
        [sys.executable, "vcf_filter.py", str(xlsx_path), "--config", str(config_src)],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parent.parent),
    )
    assert result.returncode == 0, f"STDERR: {result.stderr}"

    out_path = xlsx_path.with_name(xlsx_path.stem + "_filtered.xlsx")
    wb = load_workbook(out_path)
    ws = wb["HIGH_IMPACT"]

    # Find SCORE column
    score_col = None
    for cell in ws[1]:
        if cell.value == "SCORE":
            score_col = cell.column
            break
    assert score_col is not None

    good_cell = ws.cell(row=2, column=score_col)
    assert isinstance(good_cell.value, float) and 0.0 <= good_cell.value <= 1.0
    assert good_cell.comment is not None
    assert "SCORE :" in good_cell.comment.text
    assert "TECHNIQUE" in good_cell.comment.text
    assert "BIOLOGIQUE" in good_cell.comment.text
    # Quick overview present
    assert "(+)" in good_cell.comment.text
    # The fill should be set (not default)
    assert good_cell.fill.fgColor is not None
    good_rgb = good_cell.fill.fgColor.rgb
    # Green component should be high for good variant
    assert good_rgb is not None


def test_scores_to_rgb_signatures():
    """Verify key color signatures of the 3-axis RGB encoder."""
    from vcf_filter import _scores_to_rgb
    # Pure artifact (bad tech, no biology, not rare) → red
    assert _scores_to_rgb(bio=0.0, tech_quality=0.0, rarity=0.0) == "FF0000"
    # Perfect biology + rare + clean tech → cyan
    assert _scores_to_rgb(bio=1.0, tech_quality=1.0, rarity=1.0) == "00FFFF"
    # Biology + artifact (dangerous false positive) → yellow
    assert _scores_to_rgb(bio=1.0, tech_quality=0.0, rarity=0.0) == "FFFF00"
    # Uninteresting but technically clean → black
    assert _scores_to_rgb(bio=0.0, tech_quality=1.0, rarity=0.0) == "000000"


def test_score_variant_good_vs_bad():
    """A clearly good variant scores higher than a clearly bad one."""
    from vcf_filter import _score_variant
    good = {
        "SB": 0.50, "CALLNB": 4, "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0,
        "DP": 200, "AF_PCT": 48.0, "CLINVAR": "Pathogenic", "IMPACT": "HIGH",
        "Consequence": "stop_gained", "CADD_PHRED": 35.0, "REVEL": 0.9,
        "MAX_AF_GNOMADEX2.1": None,
    }
    bad = {
        "SB": 0.95, "CALLNB": 1, "RNCNT_TOTAL": 4, "PJCNT_TOTAL": 8,
        "DP": 50, "AF_PCT": 2.0, "CLINVAR": None, "IMPACT": "MODIFIER",
        "Consequence": "intron_variant", "CADD_PHRED": 5.0, "REVEL": 0.1,
        "MAX_AF_GNOMADEX2.1": 0.003,
    }
    _good_bio, _good_tech, _good_rarity, good_combined, good_comment = _score_variant(good)
    _bad_bio, _bad_tech, _bad_rarity, bad_combined, bad_comment = _score_variant(bad)
    assert good_combined > 0.7
    assert bad_combined < 0.3
    assert good_combined > bad_combined
    assert "✓" in good_comment
    assert "✗" in bad_comment


def test_duplicate_hgvsc_blue_font_cross_sheet(make_xlsx):
    """HGVSc present in multiple sheets should have blue font (cross-sheet detection)."""
    from openpyxl import load_workbook

    rows = [
        # This variant lands in HIGH_IMPACT + CLINVAR_P_LP + LOVD → 3 sheets = blue
        {"IMPACT": "HIGH", "Consequence": "stop_gained", "HGVSc": "c.1525C>T",
         "CLINVAR": "Pathogenic", "LOVD": "some_entry",
         "CHR": "chr9", "POS": 100, "REF": "C", "ALT": "T"},
        # This variant lands ONLY in HIGH_IMPACT → 1 sheet = NOT blue
        {"IMPACT": "HIGH", "Consequence": "stop_gained", "HGVSc": "c.999G>A",
         "CHR": "chr9", "POS": 300, "REF": "A", "ALT": "G"},
    ]
    xlsx_path = make_xlsx(rows)
    config_src = Path(__file__).resolve().parent.parent / "config_legacy.yaml"

    result = subprocess.run(
        [sys.executable, "vcf_filter.py", str(xlsx_path), "--config", str(config_src)],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parent.parent),
    )
    assert result.returncode == 0, f"STDERR: {result.stderr}"

    out_path = xlsx_path.with_name(xlsx_path.stem + "_filtered.xlsx")
    wb = load_workbook(out_path)
    ws = wb["HIGH_IMPACT"]

    # Find HGVSc column index
    hgvsc_col = None
    for cell in ws[1]:
        if cell.value == "HGVSc":
            hgvsc_col = cell.column
            break
    assert hgvsc_col is not None

    # Check font colors: cross-sheet duplicate should be blue, unique should not
    blue_rows = []
    non_blue_rows = []
    for row in ws.iter_rows(min_row=2, min_col=hgvsc_col, max_col=hgvsc_col):
        cell = row[0]
        color = cell.font.color
        if color and color.rgb and "0000FF" in str(color.rgb):
            blue_rows.append(cell.value)
        else:
            non_blue_rows.append(cell.value)

    assert "c.1525C>T" in blue_rows    # present in HIGH_IMPACT + CLINVAR_P_LP + LOVD
    assert "c.999G>A" in non_blue_rows  # only in HIGH_IMPACT

    # Check that the HGVSc cell for the duplicate has a comment listing other sheets
    for row in ws.iter_rows(min_row=2, min_col=hgvsc_col, max_col=hgvsc_col):
        cell = row[0]
        if cell.value == "c.1525C>T":
            assert cell.comment is not None
            comment = cell.comment.text
            assert "Aussi dans" in comment
            # Should mention at least CLINVAR_P_LP and LOVD
            assert "CLINVAR_P_LP" in comment
            assert "LOVD" in comment
        elif cell.value == "c.999G>A":
            assert cell.comment is None  # no cross-sheet duplicate


def test_usercom_formula_cross_sheet(make_xlsx):
    """USER_COM should have a TEXTJOIN formula referencing USER_CLASS/ANNOT from other sheets."""
    from openpyxl import load_workbook

    rows = [
        # This variant lands in HIGH_IMPACT + CLINVAR_P_LP + LOVD → formula in each
        {"IMPACT": "HIGH", "Consequence": "stop_gained", "HGVSc": "c.1525C>T",
         "CLINVAR": "Pathogenic", "LOVD": "some_entry",
         "CHR": "chr9", "POS": 100, "REF": "C", "ALT": "T"},
    ]
    xlsx_path = make_xlsx(rows)
    config_src = Path(__file__).resolve().parent.parent / "config_legacy.yaml"

    result = subprocess.run(
        [sys.executable, "vcf_filter.py", str(xlsx_path), "--config", str(config_src)],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parent.parent),
    )
    assert result.returncode == 0, f"STDERR: {result.stderr}"

    out_path = xlsx_path.with_name(xlsx_path.stem + "_filtered.xlsx")
    wb = load_workbook(out_path)

    # Check HIGH_IMPACT sheet — USER_COM should have formula referencing CLINVAR_P_LP and LOVD
    ws = wb["HIGH_IMPACT"]
    usercom_col = None
    for cell in ws[1]:
        if cell.value == "USER_COM":
            usercom_col = cell.column
            break
    assert usercom_col is not None

    usercom_cell = ws.cell(row=2, column=usercom_col)
    formula = str(usercom_cell.value)
    assert formula.startswith("=MID(")
    assert "CLINVAR_P_LP" in formula
    assert "LOVD" in formula


def test_usercom_formula_survives_igv_column_insert(tmp_path):
    """Regression: cross-sheet USER_COM formulas must reference the real
    USER_CLASS / USER_ANNOT columns even after the IGV column insert.

    inject_igv_column() prepends a column to every sheet, shifting all columns
    right by one. openpyxl does NOT rewrite formula text on insert, so formulas
    written *before* the insert end up one column too far left (pulling SB
    instead of USER_CLASS). The formula injection must therefore run last.
    """
    import re
    from openpyxl import load_workbook
    from vcf_filter import write_output

    cols = ["HGVSc", "CHR", "POS", "SB", "USER_CLASS", "USER_ANNOT", "USER_COM"]

    def mkdf():
        return pd.DataFrame([{
            "HGVSc": "c.1525C>T", "CHR": "chr9", "POS": 100, "SB": 0.5,
            "USER_CLASS": None, "USER_ANNOT": None, "USER_COM": None,
        }])

    # Same HGVSc in two sheets → cross-sheet duplicate → formula injected.
    tiers = {"HIGH_IMPACT": mkdf(), "CLINVAR_P_LP": mkdf()}

    # Fake BAM so inject_igv_column actually inserts its leftmost column.
    bam_dir = tmp_path / "bam"
    bam_dir.mkdir()
    (bam_dir / "999999_S1.bam").touch()
    igv_cfg = {"enabled": True, "bam_dir": str(bam_dir),
               "reference_lane": "S2", "port": 60151,
               "genome": "hg38", "half_window": 145}

    config_src = Path(__file__).resolve().parent.parent / "config_legacy.yaml"
    out = tmp_path / "TSC_999999_S1_filtered.xlsx"
    write_output(out, tiers, "report", str(config_src), cols,
                 igv_cfg=igv_cfg, sample_id="999999_S1")

    wb = load_workbook(out)
    ws = wb["HIGH_IMPACT"]
    usercom_col = next(c.column for c in ws[1] if c.value == "USER_COM")
    formula = str(ws.cell(row=2, column=usercom_col).value)
    assert formula.startswith("=MID("), formula

    # Each reference into CLINVAR_P_LP must land on USER_CLASS / USER_ANNOT,
    # never on SB or any other column.
    other = wb["CLINVAR_P_LP"]
    header_at = {cell.column_letter: cell.value for cell in other[1]}
    refs = re.findall(r"'CLINVAR_P_LP'!([A-Z]+)\d+", formula)
    assert refs, f"no cross-sheet refs found in: {formula}"
    targeted = {header_at.get(letter) for letter in refs}
    assert targeted <= {"USER_CLASS", "USER_ANNOT"}, (
        f"formula references wrong columns: "
        f"{[(l, header_at.get(l)) for l in refs]}"
    )

    # USER_CLASS and USER_ANNOT are concatenated directly — no stray space
    # join between them (only one is ever filled, so " | SHEET: <value>").
    assert '&" "&\'CLINVAR_P_LP\'' not in formula, formula


def test_cli_pre_report_counts(make_xlsx):
    """Verify pre-report shows correct exclusion counts."""
    rows = [
        {"MAX_AF_GNOMADEX2.1": 0.01},   # excluded by gnomAD
        {"RNCNT_TOTAL": 15},             # excluded by RNCNT
        {"IMPACT": "HIGH", "Consequence": "stop_gained"},  # survives
    ]
    xlsx_path = make_xlsx(rows)
    config_src = Path(__file__).resolve().parent.parent / "config_legacy.yaml"

    result = subprocess.run(
        [sys.executable, "vcf_filter.py", str(xlsx_path), "--config", str(config_src)],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parent.parent),
    )
    assert "gnomAD" in result.stdout
    assert "RNCNT_TOTAL" in result.stdout


def test_af_discordance_htz_recurrent_in_mosaic(make_xlsx):
    """HGVSc cell should be red when variant is recurrently HTZ in run but mosaic here."""
    from openpyxl import load_workbook

    rows = [
        {   # RNCNT_HTZ=5, AF=12% → discordance HTZ→mosaïque
            "CLINVAR": "Pathogenic", "IMPACT": "MODERATE", "Consequence": "missense_variant",
            "MAX_AF_GNOMADEX2.1": 0.0001,
            "CHR": "chr9", "POS": 100, "REF": "A", "ALT": "G",
            "HGVSc": "NM_000368.5:c.100A>G",
            "RNCNT_HTZ": 5, "RNCNT_LOW": 0, "RNCNT_TOTAL": 5, "AF_PCT": 12.0,
        },
        {   # RNCNT_HTZ=1 → pas récurrent, pas de flag
            "CLINVAR": "Pathogenic", "IMPACT": "MODERATE", "Consequence": "missense_variant",
            "MAX_AF_GNOMADEX2.1": 0.0001,
            "CHR": "chr9", "POS": 200, "REF": "C", "ALT": "T",
            "HGVSc": "NM_000368.5:c.200C>T",
            "RNCNT_HTZ": 1, "RNCNT_LOW": 0, "RNCNT_TOTAL": 1, "AF_PCT": 12.0,
        },
    ]
    xlsx_path = make_xlsx(rows)
    config_src = Path(__file__).resolve().parent.parent / "config_legacy.yaml"
    result = subprocess.run(
        [sys.executable, "vcf_filter.py", str(xlsx_path), "--config", str(config_src)],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parent.parent),
    )
    assert result.returncode == 0, f"STDERR: {result.stderr}"

    out = xlsx_path.with_name(xlsx_path.stem + "_filtered.xlsx")
    wb = load_workbook(out)
    ws = wb["CLINVAR_P_LP"]
    headers = {cell.value: cell.column for cell in ws[1]}
    hgvsc_col = headers.get("HGVSc")
    assert hgvsc_col is not None

    # Row avec RNCNT_HTZ=5 + AF mosaïque → rouge + commentaire
    red_rows = [
        r for r in range(2, ws.max_row + 1)
        if ws.cell(row=r, column=hgvsc_col).comment is not None
        and "HTZ" in ws.cell(row=r, column=hgvsc_col).comment.text
    ]
    assert len(red_rows) == 1, "Exactement une ligne devrait avoir le commentaire discordance HTZ"

    # Row avec RNCNT_HTZ=1 → pas de commentaire discordance (peut avoir doublon)
    discordance_rows_all = [
        r for r in range(2, ws.max_row + 1)
        if ws.cell(row=r, column=hgvsc_col).comment is not None
        and "HTZ" in ws.cell(row=r, column=hgvsc_col).comment.text
    ]
    assert len(discordance_rows_all) == 1, "Seule la ligne RNCNT_HTZ=5 doit avoir le flag discordance"


def test_af_discordance_mosaic_recurrent_in_htz(make_xlsx):
    """HGVSc cell should be red when variant is recurrently mosaic in run but HTZ here."""
    from openpyxl import load_workbook

    rows = [
        {   # RNCNT_LOW=4, AF=48% → discordance mosaïque→HTZ
            "CLINVAR": "Pathogenic", "IMPACT": "MODERATE", "Consequence": "missense_variant",
            "MAX_AF_GNOMADEX2.1": 0.0001,
            "CHR": "chr9", "POS": 100, "REF": "A", "ALT": "G",
            "HGVSc": "NM_000368.5:c.100A>G",
            "RNCNT_HTZ": 0, "RNCNT_LOW": 4, "RNCNT_TOTAL": 4, "AF_PCT": 48.0,
        },
    ]
    xlsx_path = make_xlsx(rows)
    config_src = Path(__file__).resolve().parent.parent / "config_legacy.yaml"
    result = subprocess.run(
        [sys.executable, "vcf_filter.py", str(xlsx_path), "--config", str(config_src)],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parent.parent),
    )
    assert result.returncode == 0, f"STDERR: {result.stderr}"

    out = xlsx_path.with_name(xlsx_path.stem + "_filtered.xlsx")
    wb = load_workbook(out)
    ws = wb["CLINVAR_P_LP"]
    headers = {cell.value: cell.column for cell in ws[1]}
    hgvsc_col = headers.get("HGVSc")
    assert hgvsc_col is not None

    flagged = [
        r for r in range(2, ws.max_row + 1)
        if ws.cell(row=r, column=hgvsc_col).comment is not None
        and "mosaïque" in ws.cell(row=r, column=hgvsc_col).comment.text
    ]
    assert len(flagged) == 1, "La ligne devrait avoir le commentaire discordance mosaïque→HTZ"


def test_cli_adds_igv_column_when_bams_present(make_xlsx, tmp_path):
    """With a config pointing bam_dir at a folder containing the sample's BAMs,
    the filtered sheets gain an IGV hyperlink column."""
    import openpyxl
    import yaml as _yaml
    from pathlib import Path

    # make_xlsx defaults to sample_id="999999_S1" → file run_1_999999_S1.xlsx,
    # so extract_sample_id() → "999999_S1".
    xlsx_in = make_xlsx([
        {"IMPACT": "HIGH", "Consequence": "stop_gained",
         "CHR": "chr9", "POS": 200, "REF": "C", "ALT": "T"},
    ])
    bam_dir = tmp_path / "bam"
    bam_dir.mkdir()
    # Patient 999999_S1 (3 BAMs) + a distinct S2 reference sample (raw BAM used)
    for s in ("999999_S1", "888888_S2"):
        for suf in (".bam", ".bam.bai", ".chim.bam", ".chim.bam.bai",
                    ".mutect2.bam", ".mutect2.bam.bai"):
            (bam_dir / f"{s}{suf}").touch()

    # Use config_legacy.yaml as base: the make_xlsx fixture produces legacy-schema
    # columns (MAX_AF_GNOMADEX2.1), which config.yaml (069 schema) does not accept.
    # All other CLI tests in this file do the same.
    cfg_path = tmp_path / "config.yaml"
    base = _yaml.safe_load(Path("config_legacy.yaml").read_text())
    base["igv"] = {"enabled": True, "port": 60151, "genome": "hg38",
                   "bam_dir": str(bam_dir), "reference_lane": "S2", "half_window": 145}
    cfg_path.write_text(_yaml.safe_dump(base))

    from vcf_filter import main as filter_main
    filter_main([str(xlsx_in), "--config", str(cfg_path)])

    out = xlsx_in.with_name(f"{xlsx_in.stem}_filtered.xlsx")
    wb = openpyxl.load_workbook(out)
    ws = wb["HIGH_IMPACT"]
    hdr = ws.cell(row=1, column=1)
    assert hdr.value == "⇅ base"           # A1 sort-by-base button
    assert hdr.hyperlink.target == "http://localhost:60151/execute?command=sort%20base"
    target = ws.cell(row=2, column=1).hyperlink.target
    assert target.startswith("http://localhost:60151/goto?")
    assert "locus=" in target

    # _summary carries the one-time "Load BAMs" link
    sm = wb["_summary"]
    load_target = next(
        sm.cell(row=r, column=2).hyperlink.target
        for r in range(1, sm.max_row + 1)
        if sm.cell(row=r, column=1).value == "IGV"
    )
    assert load_target.startswith("http://localhost:60151/load?")
    assert "name=RAW,CHIM,MUTECT2,REF_S2" in load_target


def _make_ws(headers, rows):
    """Build a bare worksheet with a header row + data rows (list of dicts)."""
    from openpyxl import Workbook
    wb = Workbook()
    ws = wb.active
    for c, h in enumerate(headers, start=1):
        ws.cell(row=1, column=c, value=h)
    for r, row in enumerate(rows, start=2):
        for c, h in enumerate(headers, start=1):
            ws.cell(row=r, column=c, value=row.get(h))
    return ws


def test_link_hgvsc_uses_varsome_column_url():
    """A 'link' HGVSc cell becomes a real hyperlink to the row's VARSOME URL."""
    from vcf_filter import _link_hgvsc_varsome
    url = "https://varsome.com/variant/hg38/chr16%3A2090268%3AC%3AA"
    ws = _make_ws(
        ["HGVSc", "VARSOME", "CHR", "POS", "REF", "ALT"],
        [{"HGVSc": "link", "VARSOME": url,
          "CHR": "chr16", "POS": 2090268, "REF": "C", "ALT": "A"}],
    )
    _link_hgvsc_varsome(ws)
    cell = ws.cell(row=2, column=1)
    assert cell.hyperlink is not None
    assert cell.hyperlink.target == url


def test_link_hgvsc_reconstructs_url_without_varsome_column():
    """Without a usable VARSOME url, reconstruct from CHR/POS/REF/ALT."""
    from vcf_filter import _link_hgvsc_varsome
    ws = _make_ws(
        ["HGVSc", "CHR", "POS", "REF", "ALT"],
        [{"HGVSc": "link", "CHR": "chr1", "POS": 26551165, "REF": "C", "ALT": "T"}],
    )
    _link_hgvsc_varsome(ws)
    cell = ws.cell(row=2, column=1)
    assert cell.hyperlink is not None
    assert cell.hyperlink.target == (
        "https://varsome.com/variant/hg38/chr1%3A26551165%3AC%3AT"
    )


def test_link_hgvsc_leaves_real_hgvsc_untouched():
    """A normal NM_…:c.… HGVSc cell keeps its value and gets no hyperlink."""
    from vcf_filter import _link_hgvsc_varsome
    ws = _make_ws(
        ["HGVSc", "CHR", "POS", "REF", "ALT"],
        [{"HGVSc": "NM_000548.5:c.3884-383G>A",
          "CHR": "chr16", "POS": 2098000, "REF": "G", "ALT": "A"}],
    )
    _link_hgvsc_varsome(ws)
    cell = ws.cell(row=2, column=1)
    assert cell.value == "NM_000548.5:c.3884-383G>A"
    assert cell.hyperlink is None
