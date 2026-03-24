"""igv_capture.py — Generate IGV screenshots for annotated TSC1/TSC2 variants."""

import argparse
import time
from pathlib import Path

import pandas as pd

import requests


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


def extract_sample_id(xlsx_path: Path) -> str:
    """Extract sample ID from filename: TSC_063_345822_S2_filtered.xlsx → 345822_S2."""
    stem = xlsx_path.stem  # TSC_063_345822_S2_filtered
    # Remove known suffixes and prefix
    stem = stem.replace("_filtered", "").replace("--beta", "")
    parts = stem.split("_")
    # Drop leading TSC + run number (e.g. TSC, 063)
    # sample_id = last two parts: DEFGEN_id + S_number
    return "_".join(parts[-2:])


def _find_index(bam: Path) -> Path | None:
    """Find BAM index: try {stem}.bai first (primary), then {file}.bam.bai (fallback)."""
    candidate1 = bam.with_suffix(".bai")            # {stem}.bai  — primary
    candidate2 = bam.parent / (bam.name + ".bai")   # {file}.bam.bai — fallback
    if candidate1.exists():
        return candidate1
    if candidate2.exists():
        return candidate2
    return None


def check_igv(port: int) -> bool:
    """Return True if IGV is reachable on the given port."""
    try:
        r = requests.get(f"http://localhost:{port}/", timeout=3)
        return r.status_code < 500
    except Exception:
        return False


def igv_cmd(port: int, command: str, **params) -> str:
    """Send a single command to IGV HTTP API. Returns response text."""
    r = requests.get(f"http://localhost:{port}/{command}", params=params, timeout=10)
    return r.text


def sanitize_filename(s: str) -> str:
    """Replace special characters unsafe for filenames."""
    return (str(s)
            .replace(">", "gt")
            .replace("/", "_")
            .replace("*", "star")
            .replace(":", "_")
            .replace(" ", "_"))


def build_variant_name(variant: dict) -> str:
    """Build a safe filename stem for a variant."""
    symbol = variant.get("SYMBOL") or "UNK"
    hgvsc = variant.get("HGVSc")
    if hgvsc and str(hgvsc).strip() not in ("", "nan", "None"):
        return f"{symbol}_{sanitize_filename(hgvsc)}"
    # Fallback to genomic coordinates
    chr_ = variant.get("CHR", "")
    pos  = variant.get("POS", "")
    ref  = variant.get("REF", "")
    alt  = variant.get("ALT", "")
    return f"{symbol}_{chr_}_{pos}_{ref}_{alt}"


def normalize_chr(chrom: str) -> str:
    """Ensure chromosome name has 'chr' prefix."""
    chrom = str(chrom).strip()
    if not chrom.startswith("chr"):
        return f"chr{chrom}"
    return chrom


def find_bams(sample_id: str, bam_dir: Path) -> dict | None:
    """Locate the 3 BAMs + indices for a sample. Returns None if any missing."""
    specs = {
        "raw":     bam_dir / f"{sample_id}.bam",
        "mutect2": bam_dir / f"{sample_id}.mutect2.bam",
        "chim":    bam_dir / f"{sample_id}.chim.bam",
    }
    result = {}
    for key, bam in specs.items():
        if not bam.exists():
            print(f"[WARN] BAM manquant : {bam}")
            return None
        bai = _find_index(bam)
        if bai is None:
            print(f"[WARN] Index BAI manquant pour : {bam}")
            return None
        result[key] = {"bam": bam, "bai": bai}
    return result
