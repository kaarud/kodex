import sys
from pathlib import Path

import pandas as pd
import pytest

# Allow imports from parent directory (07_ngs_reanalysis/)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# Columns that are sample-prefixed in the NiourK XLSX (with \n separator)
SAMPLE_COLUMNS = [
    "AF_PCT", "DP", "SB", "CALLNB",
    "CALLAF_mutect2|strelka|gatkHC|deepvariant",
    "CALLAD_mutect2", "CALLAD_strelka", "CALLAD_gatkHC", "CALLAD_deepvariant",
    "CALLFILTER_mutect2", "CALLFILTER_strelka", "CALLFILTER_gatkHC", "CALLFILTER_deepvariant",
]

# Shared annotation columns (no prefix)
SHARED_COLUMNS = [
    "Import DefGen", "SYMBOL", "REF_TRANS", "REGIONS", "HGVSc", "HGVSp",
    "USER_CLASS", "USER_ANNOT", "USER_COM",
    "IMPACT", "Consequence", "CLINVAR",
    "RNCNT_TOTAL", "PJCNT_TOTAL", "RNCNT_LOW", "PJCNT_LOW",
    "RNCNT_HTZ", "PJCNT_HTZ", "RNCNT_HOM", "PJCNT_HOM",
    "LOVD", "GNOMAD3_AF", "GNOMADEX2.1_AF", "MAX_AF_GNOMADEX2.1",
    "SPiP_Interpretation", "SPiP_InterConfident",
    "SpliceAI", "SpliceAI_num", "CADD_PHRED", "REVEL", "SIFT", "PolyPhen",
    "Existing_variation", "PUBMED", "CHR", "POS", "REF", "ALT", "Feature",
]


def _default_row():
    """A TSC1 variant that passes all exclusion filters."""
    return {
        "Import DefGen": None,
        "SYMBOL": "TSC1",
        "REF_TRANS": "NM_000368.5",
        "REGIONS": "exon15",
        "HGVSc": "c.1525C>T",
        "HGVSp": "p.Arg509Trp",
        "USER_CLASS": None,
        "USER_ANNOT": None,
        "USER_COM": None,
        "IMPACT": "HIGH",
        "Consequence": "stop_gained",
        "CLINVAR": None,
        "RNCNT_TOTAL": 0,
        "PJCNT_TOTAL": 0,
        "RNCNT_LOW": 0,
        "PJCNT_LOW": 0,
        "RNCNT_HTZ": 0,
        "PJCNT_HTZ": 0,
        "RNCNT_HOM": 0,
        "PJCNT_HOM": 0,
        "LOVD": None,
        "GNOMAD3_AF": None,
        "GNOMADEX2.1_AF": None,
        "MAX_AF_GNOMADEX2.1": None,
        "SPiP_Interpretation": None,
        "SPiP_InterConfident": None,
        "SpliceAI": None,
        "SpliceAI_num": None,
        "CADD_PHRED": 30.0,
        "REVEL": None,
        "SIFT": None,
        "PolyPhen": None,
        "Existing_variation": None,
        "PUBMED": None,
        "CHR": "chr9",
        "POS": 135786850,
        "REF": "C",
        "ALT": "T",
        "Feature": "NM_000368.5",
        # Sample columns (no prefix — added later)
        "AF_PCT": 48.5,
        "DP": 120,
        "SB": 0.52,
        "CALLNB": 3,
        "CALLAF_mutect2|strelka|gatkHC|deepvariant": "0.49|0.48|0.49|.",
        "CALLAD_mutect2": "61-59",
        "CALLAD_strelka": "62-58",
        "CALLAD_gatkHC": "60-60",
        "CALLAD_deepvariant": ".",
        "CALLFILTER_mutect2": "PASS",
        "CALLFILTER_strelka": "PASS",
        "CALLFILTER_gatkHC": "PASS",
        "CALLFILTER_deepvariant": ".",
    }


@pytest.fixture
def make_xlsx(tmp_path):
    """Factory fixture: build a NiourK-format XLSX from a list of row dicts.

    Usage in tests:
        path = make_xlsx([{...}, {...}])  # returns Path to XLSX
        path = make_xlsx()                # single default row
    """
    def _build(rows=None, sample_id="999999_S1"):
        if rows is None:
            rows = [_default_row()]

        # Build DataFrame with NiourK-style column names
        records = []
        for row in rows:
            base = _default_row()
            base.update(row)
            record = {}
            for col in SHARED_COLUMNS:
                record[col] = base[col]
            for col in SAMPLE_COLUMNS:
                niourk_name = f"{sample_id}\n{col}"
                record[niourk_name] = base[col]
            records.append(record)

        df = pd.DataFrame(records)
        xlsx_path = tmp_path / f"TSC_063_{sample_id}.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name="Variants", index=False)
            # Empty Config and Blanket sheets to mimic NiourK
            pd.DataFrame().to_excel(writer, sheet_name="Config", index=False)
            pd.DataFrame().to_excel(writer, sheet_name="Blanket", index=False)
        return xlsx_path

    return _build


@pytest.fixture
def make_filtered_xlsx(tmp_path):
    """Build a minimal filtered XLSX with multiple sheets."""
    def _build(rows_by_sheet):
        path = tmp_path / "TSC_063_345822_S2_filtered.xlsx"
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            pd.DataFrame({"Key": ["date"], "Value": ["2026-03-24"]}).to_excel(
                writer, sheet_name="_summary", index=False
            )
            for sheet, rows in rows_by_sheet.items():
                pd.DataFrame(rows).to_excel(writer, sheet_name=sheet, index=False)
        return path
    return _build
