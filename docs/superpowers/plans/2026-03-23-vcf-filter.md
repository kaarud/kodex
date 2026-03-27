# vcf_filter.py Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Python CLI script that filters NiourK XLSX variant files for TSC1/TSC2 pathogenic candidates, outputting a pre-report and a multi-sheet filtered Excel workbook.

**Architecture:** Single-script CLI (`vcf_filter.py`) reads an XLSX + a YAML config, applies sequential exclusion filters then tier-specific filters, prints a pre-report to stdout, and writes a multi-sheet Excel file. Configuration lives in `config.yaml`. Tests use pytest with a small synthetic XLSX fixture.

**Tech Stack:** Python 3.11, pandas, openpyxl, PyYAML (to install), pytest

**Spec:** `docs/superpowers/specs/2026-03-23-vcf-filter-design.md`

---

## Critical: Column Name Format in NiourK XLSX

NiourK column headers use a **sample-prefix with newline** for per-sample columns:

```
'346769_S1\nAF_PCT'      →  we want 'AF_PCT'
'346769_S1\nDP'           →  we want 'DP'
'346769_S1\nCALLNB'       →  we want 'CALLNB'
'346769_S1\nCALLAF_mutect2|strelka|gatkHC|deepvariant'  →  we want 'CALLAF_mutect2|strelka|gatkHC|deepvariant'
'346769_S1\nCALLFILTER_mutect2'                          →  we want 'CALLFILTER_mutect2'
```

Shared annotation columns have **no prefix**: `SYMBOL`, `CLINVAR`, `GNOMAD3_AF`, etc.

The script must **auto-detect and strip the sample prefix** from per-sample columns after loading.

---

## File Structure

```
07_ngs_reanalysis/
├── vcf_filter.py                    # main CLI script (all logic)
├── config.yaml                      # default configuration
├── tests/
│   ├── conftest.py                  # shared pytest fixtures (synthetic XLSX builder)
│   ├── test_load.py                 # loading + column renaming
│   ├── test_exclusions.py           # exclusion filter logic
│   ├── test_tiers.py                # tier classification logic
│   └── test_cli.py                  # end-to-end CLI test
└── docs/superpowers/
    ├── specs/2026-03-23-vcf-filter-design.md
    └── plans/2026-03-23-vcf-filter.md   # this file
```

`vcf_filter.py` is intentionally a single file. At ~300-400 lines with clear function boundaries it stays readable and avoids premature package structure for a CLI utility.

---

### Task 1: Setup — Install PyYAML and create config.yaml

**Files:**
- Create: `07_ngs_reanalysis/config.yaml`

- [ ] **Step 1: Install PyYAML in the environment**

```bash
micromamba run -n ont_bioinfo pip install pyyaml
```

Expected: `Successfully installed PyYAML-...`

- [ ] **Step 2: Create config.yaml with all spec values**

Create `07_ngs_reanalysis/config.yaml` with this exact content:

```yaml
# vcf_filter config — PRELUDE-TSC
# AF_PCT est en POURCENTAGE dans les XLSX NiourK (ex: 12.5 = 12.5 %)
# MAX_AF_GNOMADEX2.1 est en FRACTION (ex: 0.005 = 0.5 %)

genes: [TSC1, TSC2]

exclusions:
  gnomad_max_af: 0.005           # fraction — hard exclusion > 0.5 %
  rncnt_max_total: 10            # cohorte globale
  pjcnt_max_total: 4             # projet en cours (run TSC_063)
  exclude_hom: true              # RNCNT_HOM > 0 → exclu
  exclude_clinvar_benign: true

  htz_to_mosaic_artifact:
    enabled: true
    min_htz_ratio: 0.80
    min_htz_count: 3
    af_max_pct: 35.0

  mosaic_artifact:
    enabled: true
    rncnt_low_max: 5
    af_max_pct: 35.0

tiers:
  clinvar_lovd:
    enabled: true

  high:
    enabled: true
    spliceai_threshold: 0.50

  moderate:
    enabled: true
    gnomad_max_af: 0.001
    revel_min: 0.75
    cadd_min: 20

  splicing:
    enabled: true
    spliceai_min: 0.20

  mosaic:
    enabled: true
    af_min_pct: 5.0
    af_max_pct: 35.0
    rncnt_low_max: 5
    rncnt_htz_max: 3

output_columns:
  - SYMBOL
  - HGVSc
  - HGVSp
  - CHR
  - POS
  - REF
  - ALT
  - REF_TRANS
  - REGIONS
  - Feature
  - AF_PCT
  - DP
  - SB
  - CALLNB
  - IMPACT
  - Consequence
  - REVEL
  - CADD_PHRED
  - SpliceAI_num
  - SPiP_Interpretation
  - SPiP_InterConfident
  - SIFT
  - PolyPhen
  - GNOMAD3_AF
  - GNOMADEX2.1_AF
  - MAX_AF_GNOMADEX2.1
  - RNCNT_TOTAL
  - PJCNT_TOTAL
  - RNCNT_LOW
  - PJCNT_LOW
  - RNCNT_HTZ
  - PJCNT_HTZ
  - RNCNT_HOM
  - PJCNT_HOM
  - CLINVAR
  - LOVD
  - Existing_variation
  - PUBMED
  - USER_CLASS
  - USER_ANNOT
  - USER_COM
  - CALLAF_mutect2|strelka|gatkHC|deepvariant
  - CALLAD_mutect2
  - CALLAD_strelka
  - CALLAD_gatkHC
  - CALLAD_deepvariant
  - CALLFILTER_mutect2
  - CALLFILTER_strelka
  - CALLFILTER_gatkHC
  - CALLFILTER_deepvariant
```

- [ ] **Step 3: Verify YAML loads**

```bash
micromamba run -n ont_bioinfo python3 -c "import yaml; print(yaml.safe_load(open('config.yaml'))['genes'])"
```

Expected: `['TSC1', 'TSC2']`

- [ ] **Step 4: Commit**

```bash
git add config.yaml
git commit -m "feat: add vcf_filter config.yaml with ACMG/ClinGen thresholds"
```

---

### Task 2: Test fixture — Synthetic XLSX builder

**Files:**
- Create: `07_ngs_reanalysis/tests/conftest.py`

- [ ] **Step 1: Create conftest.py with a `make_xlsx` fixture**

The fixture builds a small XLSX in a temp directory with the exact NiourK column format (sample-prefixed headers). This is the shared test data factory for all subsequent tests.

```python
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
```

- [ ] **Step 2: Run pytest to confirm fixture loads**

```bash
micromamba run -n ont_bioinfo python3 -m pytest tests/conftest.py --co -q
```

Expected: `no tests ran` (fixture file only, no errors)

- [ ] **Step 3: Commit**

```bash
git add tests/conftest.py
git commit -m "test: add synthetic NiourK XLSX fixture for vcf_filter"
```

---

### Task 3: Load XLSX + strip sample prefix from column names

**Files:**
- Create: `07_ngs_reanalysis/tests/test_load.py`
- Create: `07_ngs_reanalysis/vcf_filter.py` (initial skeleton)

- [ ] **Step 1: Write failing test for load + column rename**

Create `tests/test_load.py`:

```python
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


def test_load_filters_genes(make_xlsx):
    path = make_xlsx([
        {"SYMBOL": "TSC1"},
        {"SYMBOL": "TSC2"},
        {"SYMBOL": "BRCA1"},
    ])
    df = load_variants(path, genes=["TSC1", "TSC2"])
    assert set(df["SYMBOL"].unique()) == {"TSC1", "TSC2"}
    assert len(df) == 2
```

- [ ] **Step 2: Run test to verify it fails**

```bash
micromamba run -n ont_bioinfo python3 -m pytest tests/test_load.py -v
```

Expected: FAIL — `ModuleNotFoundError: No module named 'vcf_filter'`

- [ ] **Step 3: Write minimal vcf_filter.py with load_variants**

Create `vcf_filter.py`:

```python
"""vcf_filter.py — Filter NiourK XLSX for TSC1/TSC2 pathogenic candidates."""

import argparse
import hashlib
import re
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import yaml


# ---------------------------------------------------------------------------
# 1. Load & normalise
# ---------------------------------------------------------------------------

def load_variants(xlsx_path: Path, genes: list[str] | None = None) -> pd.DataFrame:
    """Load the 'Variants' sheet and strip NiourK sample-prefix from columns.

    NiourK format: per-sample columns are named '<sample_id>\\n<column_name>'.
    We detect the sample prefix and strip it so downstream code uses clean names.
    """
    df = pd.read_excel(xlsx_path, sheet_name="Variants", engine="openpyxl")

    # Detect sample prefix: look for columns containing '\n'
    renamed = {}
    for col in df.columns:
        if "\n" in str(col):
            # Take everything after the first newline
            clean = str(col).split("\n", 1)[1]
            renamed[col] = clean
    df = df.rename(columns=renamed)

    # Filter to target genes
    if genes:
        df = df[df["SYMBOL"].isin(genes)].reset_index(drop=True)

    return df


def load_config(config_path: Path) -> dict:
    """Load YAML config file."""
    with open(config_path) as f:
        return yaml.safe_load(f)
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
micromamba run -n ont_bioinfo python3 -m pytest tests/test_load.py -v
```

Expected: 2 passed

- [ ] **Step 5: Commit**

```bash
git add vcf_filter.py tests/test_load.py
git commit -m "feat: add load_variants with NiourK column prefix stripping"
```

---

### Task 4: Exclusion filters

**Files:**
- Create: `07_ngs_reanalysis/tests/test_exclusions.py`
- Modify: `07_ngs_reanalysis/vcf_filter.py`

- [ ] **Step 1: Write failing tests for each exclusion rule**

Create `tests/test_exclusions.py`:

```python
import pandas as pd
from vcf_filter import apply_exclusions, load_config


def _cfg():
    """Minimal config dict matching config.yaml structure."""
    return {
        "exclusions": {
            "gnomad_max_af": 0.005,
            "rncnt_max_total": 10,
            "pjcnt_max_total": 4,
            "exclude_hom": True,
            "exclude_clinvar_benign": True,
            "htz_to_mosaic_artifact": {
                "enabled": True,
                "min_htz_ratio": 0.80,
                "min_htz_count": 3,
                "af_max_pct": 35.0,
            },
            "mosaic_artifact": {
                "enabled": True,
                "rncnt_low_max": 5,
                "af_max_pct": 35.0,
            },
        }
    }


def test_exclude_gnomad_too_frequent(make_xlsx):
    path = make_xlsx([
        {"MAX_AF_GNOMADEX2.1": 0.006},  # > 0.005 → excluded
        {"MAX_AF_GNOMADEX2.1": 0.004},  # kept
        {"MAX_AF_GNOMADEX2.1": None},   # novel → kept
    ])
    from vcf_filter import load_variants
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2
    assert report["gnomad"] == 1


def test_exclude_rncnt_total(make_xlsx):
    path = make_xlsx([
        {"RNCNT_TOTAL": 11},  # > 10 → excluded
        {"RNCNT_TOTAL": 10},  # kept (≤ 10)
        {"RNCNT_TOTAL": 0},   # kept
    ])
    from vcf_filter import load_variants
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2
    assert report["rncnt_total"] == 1


def test_exclude_pjcnt_total(make_xlsx):
    path = make_xlsx([
        {"PJCNT_TOTAL": 5},  # > 4 → excluded
        {"PJCNT_TOTAL": 4},  # kept (≤ 4)
    ])
    from vcf_filter import load_variants
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 1
    assert report["pjcnt_total"] == 1


def test_exclude_hom(make_xlsx):
    path = make_xlsx([
        {"RNCNT_HOM": 1},   # > 0 → excluded
        {"RNCNT_HOM": 0},   # kept
        {"RNCNT_HOM": None}, # treated as 0 → kept
    ])
    from vcf_filter import load_variants
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2
    assert report["hom"] == 1


def test_exclude_clinvar_benign(make_xlsx):
    path = make_xlsx([
        {"CLINVAR": "Benign"},
        {"CLINVAR": "Likely_benign"},
        {"CLINVAR": "Benign/Likely_benign"},
        {"CLINVAR": "Pathogenic"},
        {"CLINVAR": None},
    ])
    from vcf_filter import load_variants
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2  # Pathogenic + None
    assert report["clinvar_benign"] == 3


def test_exclude_htz_to_mosaic_artifact(make_xlsx):
    """Variant usually HTZ in cohort but low-AF in sample → artifact."""
    path = make_xlsx([
        {
            "RNCNT_HTZ": 8, "RNCNT_LOW": 1,  # ratio 8/9=0.89 > 0.80 ✓, count ≥ 3 ✓
            "AF_PCT": 15.0,  # < 35% ✓ → excluded
        },
        {
            "RNCNT_HTZ": 8, "RNCNT_LOW": 1,
            "AF_PCT": 48.0,  # > 35% → not artifact (genuinely HTZ in sample)
        },
        {
            "RNCNT_HTZ": 2, "RNCNT_LOW": 0,  # count < 3 → rule doesn't fire
            "AF_PCT": 15.0,
        },
    ])
    from vcf_filter import load_variants
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2
    assert report["htz_to_mosaic_artifact"] == 1


def test_exclude_mosaic_artifact(make_xlsx):
    """Recurrent low-AF variant = sequencing hotspot."""
    path = make_xlsx([
        {"RNCNT_LOW": 6, "AF_PCT": 12.0},  # > 5 & < 35% → excluded
        {"RNCNT_LOW": 6, "AF_PCT": 48.0},  # > 5 but AF > 35% → kept
        {"RNCNT_LOW": 3, "AF_PCT": 12.0},  # ≤ 5 → kept
    ])
    from vcf_filter import load_variants
    df = load_variants(path)
    df, report = apply_exclusions(df, _cfg())
    assert len(df) == 2
    assert report["mosaic_artifact"] == 1
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
micromamba run -n ont_bioinfo python3 -m pytest tests/test_exclusions.py -v
```

Expected: FAIL — `ImportError: cannot import name 'apply_exclusions'`

- [ ] **Step 3: Implement apply_exclusions in vcf_filter.py**

Add to `vcf_filter.py` after `load_config`:

```python
# ---------------------------------------------------------------------------
# 2. Exclusion filters
# ---------------------------------------------------------------------------

def _safe_gt(series: pd.Series, threshold, fill=0.0) -> pd.Series:
    """Compare series > threshold, treating NaN as `fill`."""
    return series.fillna(fill) > threshold


def apply_exclusions(df: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, dict]:
    """Apply sequential exclusion filters. Returns (filtered_df, report_dict).

    report_dict maps filter name → count of rows removed at that step.
    Filters are applied in spec order; counts are cumulative (each step
    operates on the survivors of the previous step).
    """
    exc = cfg["exclusions"]
    report = {}

    # 1. gnomAD > max_af (NaN = novel → keep)
    before = len(df)
    mask = _safe_gt(df["MAX_AF_GNOMADEX2.1"], exc["gnomad_max_af"], fill=0.0)
    df = df[~mask].reset_index(drop=True)
    report["gnomad"] = before - len(df)

    # 2. RNCNT_TOTAL > max
    before = len(df)
    mask = _safe_gt(df["RNCNT_TOTAL"], exc["rncnt_max_total"])
    df = df[~mask].reset_index(drop=True)
    report["rncnt_total"] = before - len(df)

    # 3. PJCNT_TOTAL > max
    before = len(df)
    mask = _safe_gt(df["PJCNT_TOTAL"], exc["pjcnt_max_total"])
    df = df[~mask].reset_index(drop=True)
    report["pjcnt_total"] = before - len(df)

    # 4. RNCNT_HOM > 0 (NaN → 0 → keep)
    if exc.get("exclude_hom", True):
        before = len(df)
        mask = _safe_gt(df["RNCNT_HOM"], 0)
        df = df[~mask].reset_index(drop=True)
        report["hom"] = before - len(df)

    # 5. ClinVar Benign / Likely_benign
    if exc.get("exclude_clinvar_benign", True):
        before = len(df)
        clinvar = df["CLINVAR"].fillna("").astype(str)
        mask = clinvar.str.contains("Benign|Likely_benign", case=False, regex=True)
        df = df[~mask].reset_index(drop=True)
        report["clinvar_benign"] = before - len(df)

    # 6. Artifact: HTZ → low-AF (false mosaic)
    art_htz = exc.get("htz_to_mosaic_artifact", {})
    if art_htz.get("enabled", False):
        before = len(df)
        htz = df["RNCNT_HTZ"].fillna(0)
        low = df["RNCNT_LOW"].fillna(0)
        total = htz + low
        ratio = htz / total.replace(0, float("nan"))  # avoid div/0
        mask = (
            (ratio > art_htz["min_htz_ratio"])
            & (htz >= art_htz["min_htz_count"])
            & (df["AF_PCT"].fillna(0) < art_htz["af_max_pct"])
        )
        df = df[~mask].reset_index(drop=True)
        report["htz_to_mosaic_artifact"] = before - len(df)

    # 7. Artifact: recurrent mosaic hotspot
    art_mos = exc.get("mosaic_artifact", {})
    if art_mos.get("enabled", False):
        before = len(df)
        mask = (
            (_safe_gt(df["RNCNT_LOW"], art_mos["rncnt_low_max"]))
            & (df["AF_PCT"].fillna(0) < art_mos["af_max_pct"])
        )
        df = df[~mask].reset_index(drop=True)
        report["mosaic_artifact"] = before - len(df)

    return df, report
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
micromamba run -n ont_bioinfo python3 -m pytest tests/test_exclusions.py -v
```

Expected: 7 passed

- [ ] **Step 5: Commit**

```bash
git add vcf_filter.py tests/test_exclusions.py
git commit -m "feat: add sequential exclusion filters with artifact detection"
```

---

### Task 5: Tier classification

**Files:**
- Create: `07_ngs_reanalysis/tests/test_tiers.py`
- Modify: `07_ngs_reanalysis/vcf_filter.py`

- [ ] **Step 1: Write failing tests for each tier**

Create `tests/test_tiers.py`:

```python
from vcf_filter import classify_tiers


def _cfg():
    return {
        "tiers": {
            "clinvar_lovd": {"enabled": True},
            "high": {"enabled": True, "spliceai_threshold": 0.50},
            "moderate": {"enabled": True, "gnomad_max_af": 0.001, "revel_min": 0.75, "cadd_min": 20},
            "splicing": {"enabled": True, "spliceai_min": 0.20},
            "mosaic": {"enabled": True, "af_min_pct": 5.0, "af_max_pct": 35.0, "rncnt_low_max": 5, "rncnt_htz_max": 3},
        }
    }


def test_tier_clinvar_pathogenic(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"CLINVAR": "Pathogenic", "IMPACT": "MODIFIER"},
        {"CLINVAR": None, "LOVD": "some_entry", "IMPACT": "MODIFIER"},
        {"CLINVAR": None, "LOVD": None, "IMPACT": "MODIFIER"},
    ])
    df = load_variants(path)
    tiers = classify_tiers(df, _cfg())
    assert len(tiers["CLINVAR_LOVD"]) == 2


def test_tier_high_impact(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"IMPACT": "HIGH", "Consequence": "stop_gained"},
        {"IMPACT": "MODERATE", "SpliceAI_num": 0.60},       # SpliceAI ≥ 0.50
        {"IMPACT": "MODERATE", "SPiP_Interpretation": "Altered", "SPiP_InterConfident": "Yes"},
        {"IMPACT": "MODERATE", "SpliceAI_num": 0.10},       # not HIGH
    ])
    df = load_variants(path)
    tiers = classify_tiers(df, _cfg())
    assert len(tiers["HIGH"]) == 3


def test_tier_moderate_missense_revel(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"IMPACT": "MODERATE", "Consequence": "missense_variant", "REVEL": 0.80, "MAX_AF_GNOMADEX2.1": 0.0005},
        {"IMPACT": "MODERATE", "Consequence": "missense_variant", "REVEL": 0.50, "MAX_AF_GNOMADEX2.1": 0.0005},  # REVEL too low
        {"IMPACT": "MODERATE", "Consequence": "missense_variant", "REVEL": None, "CADD_PHRED": 25, "MAX_AF_GNOMADEX2.1": 0.0005},  # fallback CADD
        {"IMPACT": "MODERATE", "Consequence": "missense_variant", "REVEL": None, "CADD_PHRED": 15, "MAX_AF_GNOMADEX2.1": 0.0005},  # both fail
        {"IMPACT": "MODERATE", "Consequence": "inframe_insertion", "REVEL": None, "CADD_PHRED": 22, "MAX_AF_GNOMADEX2.1": 0.0005},  # non-missense → CADD
    ])
    df = load_variants(path)
    tiers = classify_tiers(df, _cfg())
    assert len(tiers["MODERATE"]) == 3  # rows 0, 2, 4


def test_tier_moderate_gnomad_strict(make_xlsx):
    """MODERATE tier has stricter gnomAD threshold (0.001)."""
    from vcf_filter import load_variants
    path = make_xlsx([
        {"IMPACT": "MODERATE", "Consequence": "missense_variant", "REVEL": 0.80, "MAX_AF_GNOMADEX2.1": 0.002},  # > 0.001 → excluded from MODERATE
    ])
    df = load_variants(path)
    tiers = classify_tiers(df, _cfg())
    assert len(tiers["MODERATE"]) == 0


def test_tier_splicing(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"SpliceAI_num": 0.25, "IMPACT": "LOW", "Consequence": "synonymous_variant"},
        {"SpliceAI_num": 0.10, "SPiP_Interpretation": "Altered"},  # SPiP alone
        {"SpliceAI_num": 0.05, "SPiP_Interpretation": None},  # nothing
    ])
    df = load_variants(path)
    tiers = classify_tiers(df, _cfg())
    assert len(tiers["SPLICING"]) == 2


def test_tier_mosaic(make_xlsx):
    from vcf_filter import load_variants
    path = make_xlsx([
        {"AF_PCT": 15.0, "RNCNT_LOW": 2, "RNCNT_HTZ": 1},   # ✓ mosaic candidate
        {"AF_PCT": 48.0, "RNCNT_LOW": 0, "RNCNT_HTZ": 0},   # AF too high
        {"AF_PCT": 3.0, "RNCNT_LOW": 0, "RNCNT_HTZ": 0},    # AF < 5% → noise
        {"AF_PCT": 15.0, "RNCNT_LOW": 6, "RNCNT_HTZ": 0},   # RNCNT_LOW > 5 → hotspot
        {"AF_PCT": 15.0, "RNCNT_LOW": 0, "RNCNT_HTZ": 4},   # RNCNT_HTZ ≥ 3 → usually HTZ
    ])
    df = load_variants(path)
    tiers = classify_tiers(df, _cfg())
    assert len(tiers["MOSAIC"]) == 1  # only row 0
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
micromamba run -n ont_bioinfo python3 -m pytest tests/test_tiers.py -v
```

Expected: FAIL — `ImportError: cannot import name 'classify_tiers'`

- [ ] **Step 3: Implement classify_tiers in vcf_filter.py**

Add to `vcf_filter.py`:

```python
# ---------------------------------------------------------------------------
# 3. Tier classification
# ---------------------------------------------------------------------------

def classify_tiers(df: pd.DataFrame, cfg: dict) -> dict[str, pd.DataFrame]:
    """Classify surviving variants into named tiers. Returns dict of tier_name → DataFrame."""
    tcfg = cfg["tiers"]
    tiers = {}

    # Tier 1: CLINVAR_LOVD
    if tcfg.get("clinvar_lovd", {}).get("enabled", False):
        clinvar = df["CLINVAR"].fillna("").astype(str)
        lovd = df["LOVD"].fillna("").astype(str)
        mask = (
            clinvar.str.contains("Pathogenic|Likely_pathogenic", case=False, regex=True)
            | ((lovd != "") & (lovd != "."))
        )
        tiers["CLINVAR_LOVD"] = df[mask].copy()

    # Tier 2: HIGH
    if tcfg.get("high", {}).get("enabled", False):
        hcfg = tcfg["high"]
        spliceai = df["SpliceAI_num"].fillna(0)
        spip_int = df["SPiP_Interpretation"].fillna("").astype(str)
        spip_conf = df["SPiP_InterConfident"].fillna("").astype(str)
        mask = (
            (df["IMPACT"] == "HIGH")
            | (spliceai >= hcfg["spliceai_threshold"])
            | ((spip_int == "Altered") & (spip_conf == "Yes"))
        )
        tiers["HIGH"] = df[mask].copy()

    # Tier 3: MODERATE
    if tcfg.get("moderate", {}).get("enabled", False):
        mcfg = tcfg["moderate"]
        is_moderate = df["IMPACT"] == "MODERATE"
        gnomad_ok = df["MAX_AF_GNOMADEX2.1"].fillna(0) < mcfg["gnomad_max_af"]
        consequence = df["Consequence"].fillna("").astype(str)
        is_missense = consequence.str.contains("missense", case=False)
        revel = df["REVEL"]
        cadd = df["CADD_PHRED"]

        revel_pass = revel.notna() & (revel >= mcfg["revel_min"])
        cadd_pass = cadd.notna() & (cadd >= mcfg["cadd_min"])

        # Missense: require REVEL. Non-missense or REVEL absent: accept CADD.
        score_ok = (is_missense & revel_pass) | (~is_missense & cadd_pass) | (is_missense & revel.isna() & cadd_pass)

        mask = is_moderate & gnomad_ok & score_ok
        tiers["MODERATE"] = df[mask].copy()

    # Tier 4: SPLICING
    if tcfg.get("splicing", {}).get("enabled", False):
        scfg = tcfg["splicing"]
        spliceai = df["SpliceAI_num"].fillna(0)
        spip_int = df["SPiP_Interpretation"].fillna("").astype(str)
        mask = (spliceai >= scfg["spliceai_min"]) | (spip_int == "Altered")
        tiers["SPLICING"] = df[mask].copy()

    # Tier 5: MOSAIC
    if tcfg.get("mosaic", {}).get("enabled", False):
        mocfg = tcfg["mosaic"]
        af = df["AF_PCT"].fillna(0)
        rncnt_low = df["RNCNT_LOW"].fillna(0)
        rncnt_htz = df["RNCNT_HTZ"].fillna(0)
        mask = (
            (af >= mocfg["af_min_pct"])
            & (af < mocfg["af_max_pct"])
            & (rncnt_low <= mocfg["rncnt_low_max"])
            & (rncnt_htz < mocfg["rncnt_htz_max"])
        )
        tiers["MOSAIC"] = df[mask].copy()

    return tiers
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
micromamba run -n ont_bioinfo python3 -m pytest tests/test_tiers.py -v
```

Expected: 6 passed

- [ ] **Step 5: Commit**

```bash
git add vcf_filter.py tests/test_tiers.py
git commit -m "feat: add tier classification (CLINVAR_LOVD, HIGH, MODERATE, SPLICING, MOSAIC)"
```

---

### Task 6: Pre-report + Excel output + CLI

**Files:**
- Create: `07_ngs_reanalysis/tests/test_cli.py`
- Modify: `07_ngs_reanalysis/vcf_filter.py`

- [ ] **Step 1: Write failing CLI end-to-end test**

Create `tests/test_cli.py`:

```python
import subprocess
import sys
from pathlib import Path

import pandas as pd


def test_cli_produces_filtered_xlsx(make_xlsx, tmp_path):
    """End-to-end: run vcf_filter.py on a synthetic XLSX, check output file."""
    # Build input with one variant per tier
    rows = [
        # CLINVAR_LOVD tier
        {"CLINVAR": "Pathogenic", "IMPACT": "MODIFIER", "Consequence": "intron_variant",
         "CHR": "chr9", "POS": 100, "REF": "A", "ALT": "G"},
        # HIGH tier
        {"IMPACT": "HIGH", "Consequence": "stop_gained",
         "CHR": "chr9", "POS": 200, "REF": "C", "ALT": "T"},
        # MODERATE tier (missense + REVEL)
        {"IMPACT": "MODERATE", "Consequence": "missense_variant", "REVEL": 0.80,
         "MAX_AF_GNOMADEX2.1": 0.0005,
         "CHR": "chr9", "POS": 300, "REF": "G", "ALT": "A"},
        # SPLICING tier
        {"IMPACT": "LOW", "Consequence": "synonymous_variant", "SpliceAI_num": 0.30,
         "CHR": "chr9", "POS": 400, "REF": "T", "ALT": "C"},
        # MOSAIC tier
        {"AF_PCT": 15.0, "RNCNT_LOW": 1, "RNCNT_HTZ": 0, "IMPACT": "MODERATE",
         "Consequence": "missense_variant", "REVEL": 0.80, "MAX_AF_GNOMADEX2.1": 0.0005,
         "CHR": "chr9", "POS": 500, "REF": "A", "ALT": "T"},
        # Excluded: gnomAD too high
        {"MAX_AF_GNOMADEX2.1": 0.01, "IMPACT": "HIGH", "Consequence": "stop_gained",
         "CHR": "chr9", "POS": 600, "REF": "C", "ALT": "A"},
    ]
    xlsx_path = make_xlsx(rows)

    # Write a config to tmp_path
    config_src = Path(__file__).resolve().parent.parent / "config.yaml"

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
    assert "CLINVAR_LOVD" in sheets
    assert "HIGH" in sheets
    assert "MODERATE" in sheets
    assert "SPLICING" in sheets
    assert "MOSAIC" in sheets


def test_cli_pre_report_counts(make_xlsx, tmp_path, capsys):
    """Verify pre-report shows correct exclusion counts."""
    rows = [
        {"MAX_AF_GNOMADEX2.1": 0.01},   # excluded by gnomAD
        {"RNCNT_TOTAL": 15},             # excluded by RNCNT
        {"IMPACT": "HIGH", "Consequence": "stop_gained"},  # survives
    ]
    xlsx_path = make_xlsx(rows)
    config_src = Path(__file__).resolve().parent.parent / "config.yaml"

    result = subprocess.run(
        [sys.executable, "vcf_filter.py", str(xlsx_path), "--config", str(config_src)],
        capture_output=True, text=True,
        cwd=str(Path(__file__).resolve().parent.parent),
    )
    assert "gnomAD" in result.stdout
    assert "RNCNT_TOTAL" in result.stdout
```

- [ ] **Step 2: Run test to verify it fails**

```bash
micromamba run -n ont_bioinfo python3 -m pytest tests/test_cli.py -v
```

Expected: FAIL (no `main` function / CLI entry point yet)

- [ ] **Step 3: Implement pre-report, Excel writer, and CLI main**

Add to `vcf_filter.py`:

```python
# ---------------------------------------------------------------------------
# 4. Pre-report
# ---------------------------------------------------------------------------

def build_pre_report(
    total_rows: int,
    unique_variants: int,
    exclusion_report: dict,
    tiers: dict[str, pd.DataFrame],
    config_path: str,
    input_name: str,
) -> str:
    """Build the pre-report string for terminal output."""
    lines = []
    today = date.today().isoformat()
    lines.append(f"=== Pre-report {input_name} — {today} ===")
    lines.append(f"config : {config_path}")
    lines.append("")
    lines.append(f"Variants TSC1/TSC2 bruts (lignes)         : {total_rows}")
    lines.append(f"  dont variants uniques (CHR+POS+REF+ALT) : {unique_variants}")
    lines.append("")
    lines.append("--- Exclusions communes (séquentielles, sur lignes) ---")

    labels = {
        "gnomad": "gnomAD > 0.5 %",
        "rncnt_total": "RNCNT_TOTAL > seuil",
        "pjcnt_total": "PJCNT_TOTAL > seuil",
        "hom": "RNCNT_HOM > 0",
        "clinvar_benign": "ClinVar Benign",
        "htz_to_mosaic_artifact": "Artefact HTZ→low-AF",
        "mosaic_artifact": "Artefact mosaïque récurrente",
    }

    remaining = total_rows  # Start from all rows (not unique) to match apply_exclusions counts
    for key, label in labels.items():
        count = exclusion_report.get(key, 0)
        remaining -= count
        lines.append(f"  {label:<36s}: -{count:>4d}  → {remaining:>4d} restants")

    lines.append("")
    lines.append(f"\n--- Tiers (sur {remaining} lignes restantes) ---")
    for name, tier_df in tiers.items():
        n_rows = len(tier_df)
        n_unique = len(tier_df.drop_duplicates(subset=["CHR", "POS", "REF", "ALT"]))
        lines.append(f"  {name:<16s}: {n_unique:>3d} variants ({n_rows} lignes)")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 5. Excel output
# ---------------------------------------------------------------------------

def write_output(
    output_path: Path,
    tiers: dict[str, pd.DataFrame],
    pre_report: str,
    config_path: str,
    output_columns: list[str],
) -> None:
    """Write filtered tiers to a multi-sheet Excel file."""
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        # _summary sheet
        summary_df = pd.DataFrame({"Pre-report": pre_report.split("\n")})
        # Add config hash for traceability
        cfg_hash = hashlib.md5(Path(config_path).read_bytes()).hexdigest()[:8]
        meta = pd.DataFrame({
            "Key": ["config_path", "config_md5", "date"],
            "Value": [config_path, cfg_hash, date.today().isoformat()],
        })
        meta.to_excel(writer, sheet_name="_summary", index=False, startrow=0)
        summary_df.to_excel(writer, sheet_name="_summary", index=False, startrow=len(meta) + 2)

        # Tier sheets
        for name, tier_df in tiers.items():
            # Select output columns that exist in the DataFrame
            cols = [c for c in output_columns if c in tier_df.columns]
            tier_df[cols].to_excel(writer, sheet_name=name, index=False)


# ---------------------------------------------------------------------------
# 6. CLI
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Filter NiourK XLSX for TSC1/TSC2 pathogenic candidates.")
    parser.add_argument("xlsx", type=Path, help="Input NiourK XLSX file")
    parser.add_argument("--config", type=Path, default=Path("config.yaml"), help="YAML config file")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)

    # Load + gene filter
    df = load_variants(args.xlsx, genes=cfg.get("genes"))
    total_rows = len(df)
    unique_variants = len(df.drop_duplicates(subset=["CHR", "POS", "REF", "ALT"]))

    # Exclusions
    df_filtered, exc_report = apply_exclusions(df, cfg)

    # Tiers
    tiers = classify_tiers(df_filtered, cfg)

    # Pre-report
    input_name = args.xlsx.stem
    report = build_pre_report(total_rows, unique_variants, exc_report, tiers, str(args.config), input_name)
    print(report)

    # Write output
    output_path = args.xlsx.with_name(f"{args.xlsx.stem}_filtered.xlsx")
    write_output(output_path, tiers, report, str(args.config), cfg.get("output_columns", []))
    print(f"\n→ Excel écrit : {output_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run all tests**

```bash
micromamba run -n ont_bioinfo python3 -m pytest tests/ -v
```

Expected: All tests pass (test_load: 2, test_exclusions: 7, test_tiers: 6, test_cli: 2 = 17 total)

- [ ] **Step 5: Commit**

```bash
git add vcf_filter.py tests/test_cli.py
git commit -m "feat: add pre-report, Excel output, and CLI entry point"
```

---

### Task 7: Smoke test on real S1 data

**Files:** none modified (validation only)

- [ ] **Step 1: Run vcf_filter.py on the real S1 XLSX**

```bash
cd /Users/duraak/projets/PRELUDE-TSC/07_ngs_reanalysis
micromamba run -n ont_bioinfo python3 vcf_filter.py /Users/duraak/data/TSC_063/xlsx/TSC_063_346769_S1.xlsx --config config.yaml
```

Expected: Pre-report printed to terminal. Verify the numbers make sense (some variants in each tier, reasonable exclusion counts).

- [ ] **Step 2: Inspect the output Excel**

```bash
micromamba run -n ont_bioinfo python3 -c "
import pandas as pd
f = '/Users/duraak/data/TSC_063/xlsx/TSC_063_346769_S1_filtered.xlsx'
for sheet in pd.ExcelFile(f).sheet_names:
    df = pd.read_excel(f, sheet_name=sheet)
    print(f'{sheet}: {len(df)} rows, {len(df.columns)} cols')
"
```

- [ ] **Step 3: Commit final state (if any fixes were needed)**

```bash
git add vcf_filter.py tests/ config.yaml
git commit -m "fix: adjustments after smoke test on real S1 data"
```
