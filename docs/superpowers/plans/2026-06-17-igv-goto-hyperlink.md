# IGV go-to-position hyperlink + CLI — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let a clinician jump from a variant row in a `*_filtered.xlsx` (or from a CLI command) straight to that position in IGV Desktop, loading the patient's 3 BAMs plus a fixed reference ("normal") BAM for artefact comparison.

**Architecture:** A new pure-logic module `igv_link.py` builds the IGV file set (with reference resolution + dedup), the locus window, and the `http://localhost:PORT/load?...` URL — no I/O, fully unit-testable. `igv_goto.py` is a thin CLI that uses `igv_link` plus the existing socket helpers in `igv_capture.py` to drive IGV. `vcf_filter.py`'s `write_output` calls `igv_link` to add a clickable `IGV` column to each data sheet. All four BAM tracks are encoded in a single GET, so one click loads them and navigates.

**Tech Stack:** Python 3.11, openpyxl, pandas, pytest, `urllib.parse`. Reuses `igv_capture.py` (`check_igv`, `igv_cmd`, `find_bams`, `extract_sample_id`, `normalize_chr`).

---

## File Structure

- **Create `igv_link.py`** — pure functions: locus window, reference-BAM resolution, ordered BAM set assembly, IGV load-URL builder. No sockets, no IGV calls.
- **Create `igv_goto.py`** — CLI. Resolves target (sample+locus, or xlsx row), checks IGV, sends `new`/`genome`/`load`/`goto` over the socket via `igv_capture` helpers.
- **Modify `vcf_filter.py`** — `write_output` (around line 1379) injects an `IGV` hyperlink column per data sheet; `main`/`write_output` thread through the `igv:` config block and the output sample id.
- **Modify `config.yaml`** — add the `igv:` block.
- **Create `tests/test_igv_link.py`** — unit tests for the pure logic.
- **Create `tests/test_igv_goto.py`** — CLI tests (sockets mocked).
- **Modify `tests/test_cli.py`** — assert the `IGV` column appears in a generated filtered file.
- **Modify `README.md`** — document the hyperlink + CLI.

Reference resolution rule (used everywhere): given `bam_dir` and `reference_lane` (e.g. `S2`), the reference is the sample whose id ends with `_{lane}`, loading only its raw `.bam`. If that sample equals the current sample, the reference is skipped (no duplicate track).

---

## Task 1: `igv_link.build_locus` — center a window on a position

**Files:**
- Create: `igv_link.py`
- Test: `tests/test_igv_link.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_igv_link.py
def test_build_locus_centers_window():
    from igv_link import build_locus
    assert build_locus("9", 1000, half_window=145) == "chr9:855-1145"

def test_build_locus_already_prefixed():
    from igv_link import build_locus
    assert build_locus("chr16", 2098100, half_window=10) == "chr16:2098090-2098110"

def test_build_locus_clamps_start_to_one():
    from igv_link import build_locus
    assert build_locus("chrX", 50, half_window=145) == "chrX:1-195"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `micromamba run -n kodex pytest tests/test_igv_link.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'igv_link'`

- [ ] **Step 3: Write minimal implementation**

```python
# igv_link.py
"""igv_link.py — KODEX : logique pure pour piloter IGV (URL + locus + BAM set)."""

from pathlib import Path
from urllib.parse import quote

from igv_capture import normalize_chr


def build_locus(chrom: str, pos: int, half_window: int = 145) -> str:
    """Return an IGV locus string centered on pos, e.g. 'chr9:855-1145'."""
    chrom = normalize_chr(chrom)
    pos = int(pos)
    start = max(1, pos - half_window)
    end = pos + half_window
    return f"{chrom}:{start}-{end}"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `micromamba run -n kodex pytest tests/test_igv_link.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add igv_link.py tests/test_igv_link.py
git commit -m "feat: add igv_link.build_locus for centered IGV windows"
```

---

## Task 2: `igv_link.resolve_reference_bam` — find the normal track

**Files:**
- Modify: `igv_link.py`
- Test: `tests/test_igv_link.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_igv_link.py — append
def _touch_run(tmp_path):
    """Create a fake bam_dir with two samples (S1, S2), each with 3 BAMs."""
    for sid in ("100008_S1", "999111_S2"):
        for name in (f"{sid}.bam", f"{sid}.bam.bai",
                     f"{sid}.chim.bam", f"{sid}.chim.bam.bai",
                     f"{sid}.mutect2.bam", f"{sid}.mutect2.bam.bai"):
            (tmp_path / name).touch()
    return tmp_path

def test_resolve_reference_bam_by_lane(tmp_path):
    from igv_link import resolve_reference_bam
    bam_dir = _touch_run(tmp_path)
    ref = resolve_reference_bam(bam_dir, "S2", current_sample_id="100008_S1")
    assert ref.name == "999111_S2.bam"

def test_resolve_reference_bam_skips_when_current_is_reference(tmp_path):
    from igv_link import resolve_reference_bam
    bam_dir = _touch_run(tmp_path)
    ref = resolve_reference_bam(bam_dir, "S2", current_sample_id="999111_S2")
    assert ref is None

def test_resolve_reference_bam_missing_lane_returns_none(tmp_path):
    from igv_link import resolve_reference_bam
    bam_dir = _touch_run(tmp_path)
    ref = resolve_reference_bam(bam_dir, "S9", current_sample_id="100008_S1")
    assert ref is None

def test_resolve_reference_bam_ignores_chim_and_mutect2(tmp_path):
    """The raw .bam is returned, never the .chim.bam / .mutect2.bam of the lane."""
    from igv_link import resolve_reference_bam
    bam_dir = _touch_run(tmp_path)
    ref = resolve_reference_bam(bam_dir, "S2", current_sample_id="100008_S1")
    assert ref.name == "999111_S2.bam"
    assert ".chim." not in ref.name and ".mutect2." not in ref.name
```

- [ ] **Step 2: Run test to verify it fails**

Run: `micromamba run -n kodex pytest tests/test_igv_link.py -k reference -v`
Expected: FAIL with `ImportError: cannot import name 'resolve_reference_bam'`

- [ ] **Step 3: Write minimal implementation**

```python
# igv_link.py — append
def resolve_reference_bam(
    bam_dir: Path,
    reference_lane: str,
    current_sample_id: str,
) -> Path | None:
    """Return the raw .bam of the sample whose id ends with _{lane}.

    Returns None if no such sample exists, or if it is the current sample
    (avoid loading the same track twice). Only the raw BAM is considered —
    .chim.bam and .mutect2.bam are ignored.
    """
    bam_dir = Path(bam_dir)
    suffix = f"_{reference_lane}"
    for bam in sorted(bam_dir.glob("*.bam")):
        stem = bam.name[:-len(".bam")]            # strip trailing .bam
        if stem.endswith(".chim") or stem.endswith(".mutect2"):
            continue
        if not stem.endswith(suffix):
            continue
        if stem == current_sample_id:
            return None
        return bam
    return None
```

- [ ] **Step 4: Run test to verify it passes**

Run: `micromamba run -n kodex pytest tests/test_igv_link.py -k reference -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add igv_link.py tests/test_igv_link.py
git commit -m "feat: resolve reference BAM by lane with self-dedup"
```

---

## Task 3: `igv_link.collect_bam_tracks` — ordered (path, label) list

**Files:**
- Modify: `igv_link.py`
- Test: `tests/test_igv_link.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_igv_link.py — append
def test_collect_bam_tracks_patient_plus_reference(tmp_path):
    from igv_link import collect_bam_tracks
    bam_dir = _touch_run(tmp_path)   # S1 + S2 present
    tracks = collect_bam_tracks("100008_S1", bam_dir, reference_lane="S2")
    labels = [label for _, label in tracks]
    names = [p.name for p, _ in tracks]
    assert labels == ["RAW", "CHIM", "MUTECT2", "REF_S2"]
    assert names == ["100008_S1.bam", "100008_S1.chim.bam",
                     "100008_S1.mutect2.bam", "999111_S2.bam"]

def test_collect_bam_tracks_no_reference_when_self(tmp_path):
    from igv_link import collect_bam_tracks
    bam_dir = _touch_run(tmp_path)
    tracks = collect_bam_tracks("999111_S2", bam_dir, reference_lane="S2")
    labels = [label for _, label in tracks]
    assert labels == ["RAW", "CHIM", "MUTECT2"]   # reference == self → skipped

def test_collect_bam_tracks_returns_none_when_patient_bams_missing(tmp_path):
    from igv_link import collect_bam_tracks
    # empty dir → find_bams returns None
    assert collect_bam_tracks("100008_S1", tmp_path, reference_lane="S2") is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `micromamba run -n kodex pytest tests/test_igv_link.py -k collect -v`
Expected: FAIL with `ImportError: cannot import name 'collect_bam_tracks'`

- [ ] **Step 3: Write minimal implementation**

```python
# igv_link.py — add import at top alongside the existing one
from igv_capture import normalize_chr, find_bams


# igv_link.py — append
def collect_bam_tracks(
    sample_id: str,
    bam_dir: Path,
    reference_lane: str,
) -> list[tuple[Path, str]] | None:
    """Return ordered (bam_path, label) tuples to load into IGV.

    Order: patient RAW, CHIM, MUTECT2, then the reference raw BAM (REF_{lane}).
    Returns None if the patient's 3 BAMs cannot all be located.
    The reference track is omitted if missing or if it is the current sample.
    """
    bam_dir = Path(bam_dir)
    bams = find_bams(sample_id, bam_dir)
    if bams is None:
        return None
    tracks: list[tuple[Path, str]] = [
        (bams["raw"]["bam"], "RAW"),
        (bams["chim"]["bam"], "CHIM"),
        (bams["mutect2"]["bam"], "MUTECT2"),
    ]
    ref = resolve_reference_bam(bam_dir, reference_lane, sample_id)
    if ref is not None:
        tracks.append((ref, f"REF_{reference_lane}"))
    return tracks
```

- [ ] **Step 4: Run test to verify it passes**

Run: `micromamba run -n kodex pytest tests/test_igv_link.py -k collect -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add igv_link.py tests/test_igv_link.py
git commit -m "feat: assemble ordered IGV BAM track set (patient + reference)"
```

---

## Task 4: `igv_link.build_load_url` — single-GET IGV load URL

**Files:**
- Modify: `igv_link.py`
- Test: `tests/test_igv_link.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_igv_link.py — append
from pathlib import Path as _P

def test_build_load_url_full():
    from igv_link import build_load_url
    tracks = [
        (_P("/data/100008_S1.bam"), "RAW"),
        (_P("/data/100008_S1.chim.bam"), "CHIM"),
        (_P("/data/100008_S1.mutect2.bam"), "MUTECT2"),
        (_P("/data/999111_S2.bam"), "REF_S2"),
    ]
    url = build_load_url(tracks, locus="chr16:2098000-2098290",
                         port=60151, genome="hg38")
    assert url == (
        "http://localhost:60151/load?"
        "file=/data/100008_S1.bam,/data/100008_S1.chim.bam,"
        "/data/100008_S1.mutect2.bam,/data/999111_S2.bam"
        "&genome=hg38&merge=false"
        "&locus=chr16%3A2098000-2098290"
        "&name=RAW,CHIM,MUTECT2,REF_S2"
    )

def test_build_load_url_encodes_spaces_in_path():
    from igv_link import build_load_url
    tracks = [(_P("/da ta/s.bam"), "RAW")]
    url = build_load_url(tracks, locus="chr9:1-2", port=60151, genome="hg38")
    assert "/da%20ta/s.bam" in url
    assert "file=/da%20ta/s.bam&" in url
```

- [ ] **Step 2: Run test to verify it fails**

Run: `micromamba run -n kodex pytest tests/test_igv_link.py -k build_load_url -v`
Expected: FAIL with `ImportError: cannot import name 'build_load_url'`

- [ ] **Step 3: Write minimal implementation**

```python
# igv_link.py — append
def build_load_url(
    tracks: list[tuple[Path, str]],
    locus: str,
    port: int = 60151,
    genome: str = "hg38",
) -> str:
    """Build a single IGV HTTP GET that loads all tracks and navigates to locus.

    merge=false clears previously loaded tracks on each click. File paths are
    URL-encoded (slashes kept). Commas separate files/names and are left literal.
    """
    files = ",".join(quote(str(p.resolve()), safe="/") for p, _ in tracks)
    names = ",".join(label for _, label in tracks)
    locus_q = quote(locus, safe="-")          # keep '-', encode ':'
    return (
        f"http://localhost:{port}/load?"
        f"file={files}"
        f"&genome={genome}&merge=false"
        f"&locus={locus_q}"
        f"&name={names}"
    )
```

> Note: `test_build_load_url_full` uses non-existent `/data/...` paths; `Path.resolve()` on a non-existent path still returns an absolute path without touching the filesystem, and these are already absolute, so the strings are unchanged.

- [ ] **Step 4: Run test to verify it passes**

Run: `micromamba run -n kodex pytest tests/test_igv_link.py -k build_load_url -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add igv_link.py tests/test_igv_link.py
git commit -m "feat: build single-GET IGV load URL from track set"
```

---

## Task 5: `igv_goto.py` — CLI that drives IGV over the socket

**Files:**
- Create: `igv_goto.py`
- Test: `tests/test_igv_goto.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_igv_goto.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `micromamba run -n kodex pytest tests/test_igv_goto.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'igv_goto'`

- [ ] **Step 3: Write minimal implementation**

```python
# igv_goto.py
"""igv_goto.py — KODEX : charge les BAMs d'un patient (+ référence) et saute sur une position dans IGV."""

import argparse
from pathlib import Path

import pandas as pd
import yaml

from igv_capture import check_igv, igv_cmd, extract_sample_id
from igv_link import build_locus, collect_bam_tracks


def _resolve_target(args) -> tuple[str, str, int]:
    """Return (sample_id, chrom, pos) from either the sample+locus or xlsx+row form."""
    target = Path(args.target)
    if target.suffix.lower() == ".xlsx":
        if args.row is None:
            raise SystemExit("[ERREUR] --row est requis avec un fichier .xlsx")
        xl = pd.ExcelFile(target, engine="openpyxl")
        sheet = args.sheet or next(s for s in xl.sheet_names if s != "_summary")
        df = xl.parse(sheet)
        row = df.iloc[int(args.row)]
        return extract_sample_id(target), str(row["CHR"]), int(row["POS"])
    # sample + locus form: locus is "chrN:pos" or "N:pos"
    if not args.locus:
        raise SystemExit("[ERREUR] un locus 'chrN:pos' est requis avec un sample id")
    chrom, _, pos = args.locus.partition(":")
    return args.target, chrom, int(pos)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="KODEX — Charge les BAMs (patient + référence) et navigue dans IGV."
    )
    parser.add_argument("target", help="Sample id (ex: 100008_S1) ou fichier *_filtered.xlsx")
    parser.add_argument("locus", nargs="?", help="Position 'chrN:pos' (si target = sample id)")
    parser.add_argument("--sheet", default=None, help="Feuille XLSX (défaut: 1re feuille de données)")
    parser.add_argument("--row", type=int, default=None, help="Index de ligne 0-based (mode xlsx)")
    parser.add_argument("--port", type=int, default=None)
    parser.add_argument("--bam-dir", type=Path, default=None)
    parser.add_argument("--genome", default=None)
    parser.add_argument("--reference-lane", default=None)
    parser.add_argument("--half-window", type=int, default=None)
    parser.add_argument("--config", type=Path, default=Path(__file__).parent / "config.yaml")
    args = parser.parse_args(argv)

    # Config defaults (igv: block)
    cfg = {}
    if args.config.exists():
        with open(args.config) as f:
            cfg = (yaml.safe_load(f) or {}).get("igv", {}) or {}
    port = args.port or cfg.get("port", 60151)
    genome = args.genome or cfg.get("genome", "hg38")
    bam_dir = (args.bam_dir or Path(cfg.get("bam_dir", "<run_dir>/bam"))).expanduser()
    reference_lane = args.reference_lane or cfg.get("reference_lane", "S2")
    half_window = args.half_window if args.half_window is not None else cfg.get("half_window", 145)

    sample_id, chrom, pos = _resolve_target(args)

    if not check_igv(port):
        print(f"[ERREUR] IGV non joignable sur le port {port}. Ouvrir IGV Desktop et réessayer.")
        raise SystemExit(1)

    tracks = collect_bam_tracks(sample_id, bam_dir, reference_lane)
    if tracks is None:
        print(f"[ERREUR] BAMs manquants pour {sample_id} dans {bam_dir}. Abandon.")
        raise SystemExit(1)

    locus = build_locus(chrom, pos, half_window=half_window)

    igv_cmd(port, "new")
    igv_cmd(port, "genome", name=genome)
    for path, label in tracks:
        igv_cmd(port, "load", file=str(path.resolve()), name=label)
    igv_cmd(port, "goto", locus=locus)
    print(f"{len(tracks)} piste(s) chargée(s) — navigation sur {locus}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `micromamba run -n kodex pytest tests/test_igv_goto.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add igv_goto.py tests/test_igv_goto.py
git commit -m "feat: add igv_goto.py CLI for one-shot IGV load + navigate"
```

---

## Task 6: Hyperlink column injection in `vcf_filter.py`

**Files:**
- Modify: `vcf_filter.py` (`write_output` near line 1379; `main` near line 1470)
- Test: `tests/test_igv_link.py` (injection helper, openpyxl)

- [ ] **Step 1: Write the failing test**

```python
# tests/test_igv_link.py — append
import openpyxl

def test_inject_igv_column_adds_hyperlinks(tmp_path):
    from igv_link import inject_igv_column
    bam_dir = _touch_run(tmp_path)
    # Build a tiny sheet with CHR/POS headers + 2 data rows
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["SYMBOL", "CHR", "POS"])
    ws.append(["TSC1", "chr9", 1000])
    ws.append(["TSC2", "chr16", 2098100])

    inject_igv_column(ws, sample_id="100008_S1", bam_dir=bam_dir,
                      reference_lane="S2", port=60151, genome="hg38",
                      half_window=145)

    assert ws.cell(row=1, column=1).value == "IGV"      # new first column header
    link1 = ws.cell(row=2, column=1).hyperlink.target
    assert link1.startswith("http://localhost:60151/load?file=")
    assert "locus=chr9%3A855-1145" in link1
    link2 = ws.cell(row=3, column=1).hyperlink.target
    assert "locus=chr16%3A2097955-2098245" in link2

def test_inject_igv_column_noop_when_bams_missing(tmp_path):
    from igv_link import inject_igv_column
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(["SYMBOL", "CHR", "POS"])
    ws.append(["TSC1", "chr9", 1000])
    # empty bam_dir → collect_bam_tracks returns None → no column added
    added = inject_igv_column(ws, sample_id="100008_S1", bam_dir=tmp_path,
                              reference_lane="S2", port=60151, genome="hg38",
                              half_window=145)
    assert added is False
    assert ws.cell(row=1, column=1).value == "SYMBOL"   # unchanged
```

- [ ] **Step 2: Run test to verify it fails**

Run: `micromamba run -n kodex pytest tests/test_igv_link.py -k inject -v`
Expected: FAIL with `ImportError: cannot import name 'inject_igv_column'`

- [ ] **Step 3: Write minimal implementation**

```python
# igv_link.py — add to imports at top
from openpyxl.styles import Font
from openpyxl.worksheet.worksheet import Worksheet


# igv_link.py — append
def inject_igv_column(
    ws: "Worksheet",
    sample_id: str,
    bam_dir: Path,
    reference_lane: str,
    port: int = 60151,
    genome: str = "hg38",
    half_window: int = 145,
) -> bool:
    """Insert a leftmost 'IGV' column of per-row hyperlinks. Returns False (no-op)
    if the patient's BAMs are missing — the sheet is left untouched.

    All rows share the same BAM set (one sample per file); only locus varies.
    """
    tracks = collect_bam_tracks(sample_id, bam_dir, reference_lane)
    if tracks is None:
        return False

    # Locate CHR / POS columns by header
    headers = {ws.cell(row=1, column=c).value: c for c in range(1, ws.max_column + 1)}
    chr_col, pos_col = headers.get("CHR"), headers.get("POS")
    if chr_col is None or pos_col is None:
        return False

    ws.insert_cols(1)
    ws.cell(row=1, column=1, value="IGV")
    for r in range(2, ws.max_row + 1):
        chrom = ws.cell(row=r, column=chr_col + 1).value   # +1: column shifted right
        pos = ws.cell(row=r, column=pos_col + 1).value
        if chrom in (None, "") or pos in (None, ""):
            continue
        locus = build_locus(str(chrom), int(pos), half_window=half_window)
        url = build_load_url(tracks, locus, port=port, genome=genome)
        cell = ws.cell(row=r, column=1, value="▶")     # ▶
        cell.hyperlink = url
        cell.font = Font(color="0563C1", underline="single")
    return True
```

- [ ] **Step 4: Run test to verify it passes**

Run: `micromamba run -n kodex pytest tests/test_igv_link.py -k inject -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Wire into `write_output` and `main`**

In `vcf_filter.py`, add the import near the top (with the other local imports):

```python
from igv_link import inject_igv_column
from igv_capture import extract_sample_id
```

Change the `write_output` signature (line 1354) to accept an optional `igv_cfg` and `sample_id`:

```python
def write_output(
    output_path: Path,
    tiers: dict[str, pd.DataFrame],
    pre_report: str,
    config_path: str,
    output_columns: list[str],
    scoring: dict | None = None,
    discordance: dict | None = None,
    af_column: str = "MAX_AF_GNOMADEX2.1",
    igv_cfg: dict | None = None,
    sample_id: str | None = None,
) -> None:
```

Add a **separate loop after** the cross-sheet duplicate/formula block (after line ~1391), so the leftmost-column insert can never interfere with the HGVSc column lookups in `_collect_hgvsc_across_sheets` / `_mark_duplicate_hgvsc` / `_inject_usercom_formulas`:

```python
        # IGV hyperlink column — last, so column inserts don't shift earlier lookups
        if igv_cfg and igv_cfg.get("enabled") and sample_id:
            bam_dir = Path(igv_cfg.get("bam_dir", "<run_dir>/bam")).expanduser()
            added_any = False
            for name in tiers:
                if inject_igv_column(
                    writer.sheets[name],
                    sample_id=sample_id,
                    bam_dir=bam_dir,
                    reference_lane=igv_cfg.get("reference_lane", "S2"),
                    port=igv_cfg.get("port", 60151),
                    genome=igv_cfg.get("genome", "hg38"),
                    half_window=igv_cfg.get("half_window", 145),
                ):
                    added_any = True
            if not added_any:
                print(f"[IGV] Colonne omise — BAMs introuvables pour {sample_id} dans {bam_dir}")
```

In `main` (the `write_output(...)` call at line 1470), pass the new args:

```python
    write_output(output_path, sheets, report, str(config_path), cfg.get("output_columns", []),
                 scoring=cfg.get("scoring"),
                 discordance=cfg.get("af_discordance"),
                 af_column=af_col,
                 igv_cfg=cfg.get("igv"),
                 sample_id=extract_sample_id(output_path))
```

> The injection loop runs strictly after the cross-sheet block (lines 1387-1391), so the leftmost-column insert cannot shift the columns those functions read. Step 7 still runs the full suite to confirm no regression.

- [ ] **Step 6: Commit**

```bash
git add igv_link.py vcf_filter.py tests/test_igv_link.py
git commit -m "feat: inject per-row IGV hyperlink column into filtered sheets"
```

- [ ] **Step 7: Regression check on cross-sheet formatting**

Run: `micromamba run -n kodex pytest tests/ -v`
Expected: PASS (all existing tests, including `test_tiers.py` cross-sheet/HGVSc tests). If any HGVSc-column test fails because it assumed a fixed column index, fix that test to look up the column by header name (matching how the code already resolves columns).

---

## Task 7: Config block + end-to-end column assertion

**Files:**
- Modify: `config.yaml`
- Modify: `tests/test_cli.py`

- [ ] **Step 1: Add the `igv:` block to `config.yaml`**

Append after the `igv_exclusion_terms:` block (line ~199):

```yaml
igv:
  enabled: true
  port: 60151
  genome: hg38
  bam_dir: <run_dir>/bam      # emplacement des BAMs locaux
  reference_lane: S2               # piste « normale » toujours chargée (suffixe _S2)
  half_window: 145
```

Also add the same block (with `enabled: true`) to `config_legacy.yaml` so the legacy (run 1) path gets hyperlinks too.

- [ ] **Step 2: Write the failing end-to-end test**

```python
# tests/test_cli.py — append
def test_cli_adds_igv_column_when_bams_present(make_xlsx, tmp_path):
    """With a config pointing bam_dir at a folder containing the sample's BAMs,
    the filtered sheets gain an IGV hyperlink column."""
    import openpyxl
    import yaml as _yaml
    from pathlib import Path

    # make_xlsx defaults to sample_id="999999_S1" → file run_1_999999_S1.xlsx
    # (conftest.py:105,123). So extract_sample_id() → "999999_S1".
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

    cfg_path = tmp_path / "config.yaml"
    base = _yaml.safe_load(Path("config.yaml").read_text())
    base["igv"] = {"enabled": True, "port": 60151, "genome": "hg38",
                   "bam_dir": str(bam_dir), "reference_lane": "S2", "half_window": 145}
    cfg_path.write_text(_yaml.safe_dump(base))

    from vcf_filter import main as filter_main
    filter_main([str(xlsx_in), "--config", str(cfg_path)])

    out = xlsx_in.with_name(f"{xlsx_in.stem}_filtered.xlsx")
    wb = openpyxl.load_workbook(out)
    ws = wb["HIGH_IMPACT"]
    assert ws.cell(row=1, column=1).value == "IGV"
    target = ws.cell(row=2, column=1).hyperlink.target
    assert target.startswith("http://localhost:60151/load?")
    assert "name=RAW,CHIM,MUTECT2,REF_S2" in target
```

- [ ] **Step 3: Run test to verify it fails, then passes**

Run: `micromamba run -n kodex pytest tests/test_cli.py::test_cli_adds_igv_column_when_bams_present -v`
Expected: first FAIL if the `igv:` block / wiring is incomplete, then PASS after Step 1 + Task 6 wiring are in place.

- [ ] **Step 4: Commit**

```bash
git add config.yaml config_legacy.yaml tests/test_cli.py
git commit -m "feat: enable IGV hyperlink column via config (069 + legacy)"
```

---

## Task 8: Documentation

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Document the feature**

Add a subsection under the existing IGV usage area (near the "Captures IGV" / "Bookmark IGV" sections, ~line 97-108):

````markdown
### Lien IGV par position (go-to)

Chaque feuille du `*_filtered.xlsx` reçoit une colonne **IGV** : un clic charge dans
IGV Desktop les 3 BAMs du patient (RAW + CHIM + MUTECT2) plus le BAM brut d'un
échantillon de référence (« normal », par défaut lane `S2`) et navigue sur la
position du variant. Un softclip présent à la fois sur le patient et la référence
est probablement un artefact.

IGV doit être ouvert (port 60151). Les BAMs sont lus en local depuis `igv.bam_dir`.

Équivalent en ligne de commande :

```bash
# par sample + position
micromamba run -n kodex python igv_goto.py 100008_S1 chr16:2098100

# depuis une ligne d'un fichier filtré (index 0-based)
micromamba run -n kodex python igv_goto.py sample_filtered.xlsx --sheet HIGH_IMPACT --row 7
```

Configuration (`config.yaml`, bloc `igv:`) : `port`, `genome`, `bam_dir`,
`reference_lane`, `half_window`.
````

Also add `igv_goto.py` and `igv_link.py` to the file-list section (near lines 153-154, alongside `igv_capture.py` / `igv_bookmark.py`).

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "docs: document IGV go-to hyperlink and igv_goto CLI"
```

---

## Final verification

- [ ] Run the full suite: `micromamba run -n kodex pytest tests/ -v` — all green.
- [ ] Manual smoke (optional, needs IGV + real BAMs): open IGV, run
  `python igv_goto.py <real_sample> <chr:pos>`, confirm 3+1 tracks load and the
  view centers on the position. **Do not run KODEX filtering on a real input whose
  `_filtered.xlsx` you keep — it overwrites silently** (see project memory).
```
