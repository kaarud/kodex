# IGV Capture Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Créer `igv_capture.py` — script autonome qui lit un XLSX filtré, trouve les variants annotés par l'utilisateur, et génère 2 captures IGV par variant (strand bias + soft clips) avec les 3 BAMs empilés via l'API HTTP port d'IGV Desktop.

**Architecture:** Script unique `igv_capture.py` avec 5 fonctions pures testables + CLI argparse. Toutes les fonctions qui appellent IGV reçoivent le port en paramètre pour être mockables. TDD sur toutes les fonctions sauf les appels HTTP réels à IGV (mockés avec `unittest.mock`).

**Tech Stack:** Python 3.11, pandas, openpyxl, requests, pathlib, argparse, unittest.mock (tests)

---

## File Structure

```
07_ngs_reanalysis/
  igv_capture.py                      ← script principal (créer)
  tests/
    test_igv_capture.py               ← tests unitaires (créer)
```

---

### Task 1 : Chargement des variants annotés depuis le XLSX

**Files:**
- Create: `igv_capture.py` (squelette + fonction `load_annotated_variants`)
- Create: `tests/test_igv_capture.py`

- [ ] **Step 1 : Écrire le test**

Ajouter le fixture dans `tests/conftest.py` (déjà existant) pour partager avec tous les tests :

```python
# Ajouter à la fin de tests/conftest.py
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
    assert len(variants) == 1  # dédupliqué


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
```

- [ ] **Step 2 : Vérifier que les tests échouent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py::test_load_annotated_variants_single_sheet -v
```
Expected: `ImportError: cannot import name 'load_annotated_variants'`

- [ ] **Step 3 : Implémenter `load_annotated_variants`**

```python
# igv_capture.py
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
```

- [ ] **Step 4 : Vérifier que les tests passent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -k "load_annotated" -v
```
Expected: 3 PASSED

- [ ] **Step 5 : Commit**

```bash
git add igv_capture.py tests/test_igv_capture.py
git commit -m "feat(igv): add load_annotated_variants — reads all sheets, deduplicates by position"
```

---

### Task 2 : Extraction du sample_id et découverte des BAMs

**Files:**
- Modify: `igv_capture.py`
- Modify: `tests/test_igv_capture.py`

- [ ] **Step 1 : Écrire les tests**

```python
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
    """Accept {file}.bai if {file}.bam.bai absent."""
    from igv_capture import find_bams
    sid = "345822_S2"
    (tmp_path / f"{sid}.bam").touch()
    (tmp_path / f"{sid}.bai").touch()          # .bai not .bam.bai
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
```

- [ ] **Step 2 : Vérifier que les tests échouent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -k "sample_id or bam" -v
```
Expected: ImportError

- [ ] **Step 3 : Implémenter `extract_sample_id` et `find_bams`**

```python
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


def find_bams(sample_id: str, bam_dir: Path) -> dict | None:
    """Locate the 3 BAMs + indices for a sample. Returns None if any missing."""
    specs = {
        "raw":    bam_dir / f"{sample_id}.bam",
        "mutect2": bam_dir / f"{sample_id}.mutect2.bam",
        "chim":   bam_dir / f"{sample_id}.chim.bam",
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
```

- [ ] **Step 4 : Vérifier que les tests passent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -k "sample_id or bam" -v
```
Expected: 4 PASSED

- [ ] **Step 5 : Commit**

```bash
git add igv_capture.py tests/test_igv_capture.py
git commit -m "feat(igv): add extract_sample_id and find_bams with BAI fallback"
```

---

### Task 3 : Client IGV (API HTTP)

**Files:**
- Modify: `igv_capture.py`
- Modify: `tests/test_igv_capture.py`

- [ ] **Step 1 : Écrire les tests (avec mock HTTP)**

```python
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
```

- [ ] **Step 2 : Vérifier que les tests échouent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -k "igv" -v
```
Expected: ImportError

- [ ] **Step 3 : Implémenter `check_igv` et `igv_cmd`**

```python
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
```

- [ ] **Step 4 : Vérifier que les tests passent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -k "igv" -v
```
Expected: 4 PASSED

- [ ] **Step 5 : Commit**

```bash
git add igv_capture.py tests/test_igv_capture.py
git commit -m "feat(igv): add IGV HTTP client (check_igv, igv_cmd)"
```

---

### Task 4 : Utilitaires (sanitisation + normalisation CHR)

**Files:**
- Modify: `igv_capture.py`
- Modify: `tests/test_igv_capture.py`

- [ ] **Step 1 : Écrire les tests**

```python
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
```

- [ ] **Step 2 : Vérifier que les tests échouent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -k "sanitize or variant_name or normalize" -v
```
Expected: ImportError

- [ ] **Step 3 : Implémenter les utilitaires**

```python
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
```

- [ ] **Step 4 : Vérifier que les tests passent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -k "sanitize or variant_name or normalize" -v
```
Expected: 4 PASSED

- [ ] **Step 5 : Commit**

```bash
git add igv_capture.py tests/test_igv_capture.py
git commit -m "feat(igv): add filename sanitization and CHR normalization utilities"
```

---

### Task 5 : Initialisation IGV et capture par variant

**Files:**
- Modify: `igv_capture.py`
- Modify: `tests/test_igv_capture.py`

- [ ] **Step 1 : Écrire les tests**

```python
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
```

- [ ] **Step 2 : Vérifier que les tests échouent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -k "init_igv or capture_variant" -v
```
Expected: ImportError

- [ ] **Step 3 : Implémenter `init_igv` et `capture_variant`**

```python
_HALF_WINDOW = 145   # ±145 pb → 290 pb total


def init_igv(port: int, bams: dict, genome: str = "hg38") -> None:
    """Reset IGV session, load genome and 3 BAMs."""
    igv_cmd(port, "new")
    igv_cmd(port, "genome", name=genome)
    labels = {"raw": "RAW", "mutect2": "MUTECT2", "chim": "CHIM"}
    for key, paths in bams.items():
        igv_cmd(port, "load",
                file=str(paths["bam"].resolve()),
                index=str(paths["bai"].resolve()),
                name=labels[key])


def capture_variant(
    port: int,
    variant: dict,
    out_dir: Path,
    delay: float = 2.0,
) -> dict:
    """Generate strand-bias and soft-clip screenshots for one variant."""
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    chrom = normalize_chr(str(variant.get("CHR", "")))
    pos   = int(variant["POS"])
    dp    = int(variant.get("DP") or 100)
    name  = build_variant_name(variant)

    panel_height = dp * 2 + 150
    locus = f"{chrom}:{pos - _HALF_WINDOW}-{pos + _HALF_WINDOW}"

    # Set display mode and panel height
    igv_cmd(port, "setPreference", name="SAM.DISPLAY_MODE",   value="SQUISHED")
    igv_cmd(port, "setPreference", name="SAM.SAMPLING_COUNT", value=str(dp))
    igv_cmd(port, "maxPanelHeight", value=str(panel_height))
    igv_cmd(port, "goto", locus=locus)
    time.sleep(delay)

    # Capture 1 — strand bias
    strand_png = out_dir / f"{name}_strand.png"
    igv_cmd(port, "setPreference", name="SAM.COLOR_BY", value="READ_STRAND")
    igv_cmd(port, "snapshot", filename=str(strand_png))

    # Capture 2 — soft clips (strand ALWAYS before soft clips — reset depends on this order)
    softclip_png = out_dir / f"{name}_softclip.png"
    igv_cmd(port, "setPreference", name="SAM.SHOW_SOFT_CLIPPED", value="true")
    igv_cmd(port, "setPreference", name="SAM.COLOR_BY",          value="NO_COLORING")
    igv_cmd(port, "snapshot", filename=str(softclip_png))

    # Reset for next variant
    igv_cmd(port, "setPreference", name="SAM.SHOW_SOFT_CLIPPED", value="false")

    return {"strand": strand_png, "softclip": softclip_png}
```

- [ ] **Step 4 : Vérifier que les tests passent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -k "init_igv or capture_variant" -v
```
Expected: 2 PASSED

- [ ] **Step 5 : Commit**

```bash
git add igv_capture.py tests/test_igv_capture.py
git commit -m "feat(igv): add init_igv and capture_variant with dynamic panel height"
```

---

### Task 6 : CLI et orchestration principale

**Files:**
- Modify: `igv_capture.py`
- Modify: `tests/test_igv_capture.py`

- [ ] **Step 1 : Écrire le test CLI**

```python
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
```

- [ ] **Step 2 : Vérifier que les tests échouent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -k "main" -v
```
Expected: ImportError sur `main`

- [ ] **Step 3 : Implémenter `main`**

```python
def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Génère des captures IGV pour les variants annotés d'un XLSX filtré."
    )
    parser.add_argument("xlsx", type=Path, help="Fichier *_filtered.xlsx")
    parser.add_argument("--port",    type=int,  default=60151)
    parser.add_argument("--bam-dir", type=Path,
                        default=Path("~/data/TSC_063/bam").expanduser())
    parser.add_argument("--out",     type=Path,
                        default=Path("~/data/TSC_063/igv_captures").expanduser())
    parser.add_argument("--delay",   type=float, default=2.0)
    parser.add_argument("--genome",  type=str,   default="hg38")
    args = parser.parse_args(argv)

    xlsx_path = args.xlsx.expanduser().resolve()
    bam_dir   = args.bam_dir.expanduser().resolve()
    out_dir   = args.out.expanduser().resolve()

    # 1. Charger les variants annotés
    variants = load_annotated_variants(xlsx_path)
    if not variants:
        print("Aucun variant annoté trouvé (USER_CLASS / USER_ANNOT / USER_COM vides).")
        raise SystemExit(0)
    print(f"{len(variants)} variant(s) annoté(s) à capturer.")

    # 2. Vérifier IGV
    if not check_igv(args.port):
        print(f"[ERREUR] IGV non joignable sur le port {args.port}. "
              "Ouvrir IGV Desktop et réessayer.")
        raise SystemExit(1)

    # 3. Trouver les BAMs
    sample_id = extract_sample_id(xlsx_path)
    bams = find_bams(sample_id, bam_dir)
    if bams is None:
        print(f"[ERREUR] BAMs manquants pour {sample_id} dans {bam_dir}. Abandon.")
        raise SystemExit(1)

    # 4. Initialiser IGV
    sample_out = out_dir / sample_id
    print(f"Initialisation IGV (genome={args.genome}, 3 BAMs chargés)...")
    init_igv(args.port, bams, genome=args.genome)

    # 5. Capturer chaque variant
    for i, variant in enumerate(variants, 1):
        chrom = variant.get("CHR", "")
        pos   = variant.get("POS", "")
        if not chrom or not pos:
            print(f"[WARN] Variant {i} sans CHR/POS — ignoré.")
            continue
        name = build_variant_name(variant)
        print(f"  [{i}/{len(variants)}] {name} @ {chrom}:{pos}")
        result = capture_variant(args.port, variant, sample_out, delay=args.delay)
        print(f"    → {result['strand'].name}")
        print(f"    → {result['softclip'].name}")

    print(f"\nCaptures terminées → {sample_out}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4 : Vérifier que tous les tests passent**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/test_igv_capture.py -v
```
Expected: tous PASSED

- [ ] **Step 5 : Commit**

```bash
git add igv_capture.py tests/test_igv_capture.py
git commit -m "feat(igv): add CLI orchestration — complete igv_capture.py"
```

---

### Task 7 : Alias bash et vérification end-to-end

**Files:**
- Modify: `~/.bashrc`

- [ ] **Step 1 : Ajouter l'alias**

Ajouter dans `~/.bashrc` après le bloc `vcf_filter` :

```bash
# === IGV CAPTURE (PRELUDE-TSC) ===
# igv_capture <*_filtered.xlsx> [--port 60151] [--delay 2.0]
alias igv_capture="~/bin/micromamba run -n ont_bioinfo python ~/projets/PRELUDE-TSC/07_ngs_reanalysis/igv_capture.py"
```

- [ ] **Step 2 : Recharger et vérifier**

```bash
source ~/.bashrc
igv_capture --help
```
Expected: affichage de l'aide argparse

- [ ] **Step 3 : Lancer la suite de tests complète**

```bash
eval "$(~/bin/micromamba shell hook -s bash)" && micromamba activate ont_bioinfo && \
python -m pytest tests/ -v
```
Expected: 24 tests existants + nouveaux tests igv_capture = tous PASSED

- [ ] **Step 4 : Commit final**

Note: `~/.bashrc` est hors du dépôt git — ne pas l'ajouter au commit.

```bash
git add igv_capture.py tests/test_igv_capture.py tests/conftest.py
git commit -m "chore: igv_capture complete — all tests pass"
```
