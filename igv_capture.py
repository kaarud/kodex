"""igv_capture.py — Generate IGV screenshots for annotated TSC1/TSC2 variants."""

import argparse
import time
from pathlib import Path

import pandas as pd

try:
    import requests
except ImportError:
    requests = None  # type: ignore[assignment]


_ANNOTATION_COLS = ["USER_CLASS", "USER_ANNOT", "USER_COM"]
_DEDUP_KEY = ["CHR", "POS", "REF", "ALT"]


def load_annotated_variants(xlsx_path: Path) -> list[dict]:
    """Read all data sheets (skip _summary), return annotated variants deduplicated."""
    xl = pd.ExcelFile(xlsx_path, engine="openpyxl")
    sheets = [s for s in xl.sheet_names if s != "_summary"]

    frames = []
    for sheet in sheets:
        df = xl.parse(sheet)
        # Keep only rows with at least one annotation field non-empty
        mask = pd.Series(False, index=df.index)
        for col in _ANNOTATION_COLS:
            if col in df.columns:
                mask |= df[col].notna() & (df[col].astype(str).str.strip() != "")
        frames.append(df[mask])

    if not frames:
        return []

    combined = pd.concat(frames, ignore_index=True)

    # Deduplicate by genomic position
    available_key = [c for c in _DEDUP_KEY if c in combined.columns]
    if available_key:
        combined = combined.drop_duplicates(subset=available_key)

    return combined.to_dict(orient="records")
