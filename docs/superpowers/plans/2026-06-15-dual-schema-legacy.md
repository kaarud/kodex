# KODEX Dual-Schema (063/069) via `--legacy` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `vcf_filter.py` process both the new NiourK run 2 schema (default) and the old run 1 schema (via `--legacy`), with all schema differences expressed in YAML config rather than code branches.

**Architecture:** Internal column vocabulary stays as today's 063 names. At load time, `column_aliases` renames the raw 069 columns onto those internal names, `synthesize_user_columns` fabricates the absent USER_* columns, and `required_columns` validation produces a clean schema-mismatch error instead of a `KeyError`. A single `af_filter_column` config key replaces the six hard-coded `MAX_AF_GNOMADEX2.1` references so frequency filtering, rarity scoring, and AF cell-coloring all read the schema-appropriate gnomAD column. `config.yaml` carries 069 defaults; a new `config_legacy.yaml` carries the 063 values; `--legacy` is shorthand for `--config config_legacy.yaml`.

**Tech Stack:** Python 3.11, pandas, openpyxl, PyYAML, argparse, pytest. Run tests with `python -m pytest` (bare `python` is not on PATH).

**Reference spec:** `docs/superpowers/specs/2026-06-15-dual-schema-legacy-design.md`

---

## File Structure

- **Modify** `vcf_filter.py`
  - `load_variants(...)` (`:20`) — new kwargs `column_aliases`, `synthesize_user_columns`, `required_columns`; alias rename + USER_* synthesis + validation.
  - New module-level `SchemaMismatchError` + `_USER_COLUMNS` constant.
  - `apply_exclusions(df, cfg)` (`:69`) — read AF column from `cfg["af_filter_column"]`.
  - `classify_sheets(df, cfg)` (`:130`) — same, three sites.
  - `_cell_fill(...)` (`:327`) + `_apply_cell_colors(...)` (`:1251`) — thread `af_column`.
  - `_score_variant(...)` (`:532`) + `_apply_variant_scores(...)` (`:900`) — thread `af_column`.
  - `write_output(...)` (`:1314`) — accept + forward `af_column`.
  - `main(...)` (`:1357`) — `--legacy` flag, config-path resolution, AF-column wiring, clean `SchemaMismatchError` exit. (The `is_file()` guard already present here gets committed as part of this task.)
- **Create** `config_legacy.yaml` — copy of today's `config.yaml` + three legacy-valued keys.
- **Modify** `config.yaml` — add the three keys with 069 values, rewrite `output_columns` for 069.
- **Modify** `tests/conftest.py` — add `make_xlsx_069` factory fixture + `_USER_COLS`/`_ALIAS_069` helpers.
- **Create** `tests/test_dual_schema.py` — unit + integration + error-path tests.

---

## Task 1: `load_variants` — aliases, USER_* synthesis, schema validation

**Files:**
- Modify: `vcf_filter.py:16-51` (region around `load_variants`)
- Modify: `tests/conftest.py` (add 069 factory)
- Test: `tests/test_dual_schema.py` (new)

- [ ] **Step 1: Add the 069 factory fixture to `tests/conftest.py`**

Append these module-level constants after the existing `SHARED_COLUMNS` block (after line 33) and the fixture after the existing `make_xlsx` fixture (after line 125):

```python
# --- 069 schema helpers -----------------------------------------------------
# 069 renames these Variants columns and drops the USER_* columns entirely
# (they live in the Valideur*/Technicien sheets instead). gnomAD v4 is added.
_ALIAS_069 = {"SYMBOL": "GENE", "Consequence": "CSQ", "REGIONS": "REGION"}
_USER_COLS = ["USER_CLASS", "USER_ANNOT", "USER_COM"]


@pytest.fixture
def make_xlsx_069(tmp_path):
    """Factory fixture: build a NiourK run 2-format XLSX from row dicts.

    Columns use 069 raw names (GENE/CSQ/REGION), USER_* are absent, and a
    GNOMAD4_AF column is present. Set GNOMAD4_AF per-row to drive filtering;
    if unset it mirrors MAX_AF_GNOMADEX2.1.
    """
    def _build(rows=None, sample_id="100005_S1"):
        if rows is None:
            rows = [{}]

        records = []
        for row in rows:
            base = _default_row()
            base.update(row)
            gnomad4 = base.get("GNOMAD4_AF", base["MAX_AF_GNOMADEX2.1"])
            record = {}
            for col in SHARED_COLUMNS:
                if col in _USER_COLS:
                    continue
                record[_ALIAS_069.get(col, col)] = base[col]
            record["GNOMAD4_AF"] = gnomad4
            for col in SAMPLE_COLUMNS:
                record[f"{sample_id}\n{col}"] = base[col]
            records.append(record)

        df = pd.DataFrame(records)
        xlsx_path = tmp_path / f"run_2_Singularity_{sample_id}.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            df.to_excel(writer, sheet_name="Variants", index=False)
            for sheet in ("Config", "Blanket", "Valideur1",
                          "Valideur2", "Technicien", "AIDE"):
                pd.DataFrame().to_excel(writer, sheet_name=sheet, index=False)
        return xlsx_path

    return _build
```

- [ ] **Step 2: Write the failing unit tests**

Create `tests/test_dual_schema.py`:

```python
import pytest

from vcf_filter import SchemaMismatchError, load_variants

ALIASES_069 = {"GENE": "SYMBOL", "CSQ": "Consequence", "REGION": "REGIONS"}


def test_load_069_applies_aliases(make_xlsx_069):
    path = make_xlsx_069([{"SYMBOL": "TSC2"}])  # SYMBOL→stored as GENE by factory
    df = load_variants(
        path,
        genes=["TSC1", "TSC2"],
        column_aliases=ALIASES_069,
        synthesize_user_columns=True,
        required_columns=["SYMBOL", "GNOMAD4_AF"],
    )
    assert "SYMBOL" in df.columns
    assert "Consequence" in df.columns
    assert "REGIONS" in df.columns
    assert "GENE" not in df.columns
    assert df["SYMBOL"].iloc[0] == "TSC2"


def test_load_069_synthesizes_empty_user_columns(make_xlsx_069):
    path = make_xlsx_069()
    df = load_variants(
        path,
        genes=["TSC1", "TSC2"],
        column_aliases=ALIASES_069,
        synthesize_user_columns=True,
        required_columns=["SYMBOL", "GNOMAD4_AF"],
    )
    for col in ("USER_CLASS", "USER_ANNOT", "USER_COM"):
        assert col in df.columns
        assert df[col].isna().all()


def test_load_069_missing_required_column_raises(make_xlsx_069):
    # No aliases passed → GENE never renamed to SYMBOL → mismatch
    path = make_xlsx_069()
    with pytest.raises(SchemaMismatchError) as exc:
        load_variants(
            path,
            genes=["TSC1", "TSC2"],
            column_aliases=None,
            synthesize_user_columns=False,
            required_columns=["SYMBOL"],
        )
    assert "SYMBOL" in str(exc.value)
    assert "--legacy" in str(exc.value)


def test_load_063_unchanged_under_legacy_defaults(make_xlsx):
    # Legacy path: no aliases, no synthesis — behaves exactly as before.
    path = make_xlsx([{"SYMBOL": "TSC1"}, {"SYMBOL": "TSC2"}])
    df = load_variants(
        path,
        genes=["TSC1", "TSC2"],
        column_aliases={},
        synthesize_user_columns=False,
        required_columns=["SYMBOL", "MAX_AF_GNOMADEX2.1"],
    )
    assert set(df["SYMBOL"].unique()) == {"TSC1", "TSC2"}
    assert "USER_CLASS" in df.columns  # 063 already has it
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `python -m pytest tests/test_dual_schema.py -q`
Expected: FAIL — `ImportError: cannot import name 'SchemaMismatchError'`.

- [ ] **Step 4: Implement in `vcf_filter.py`**

Add the exception class and constant immediately after the imports block (after line 13):

```python
class SchemaMismatchError(Exception):
    """Raised when input columns don't match the active schema profile."""


_USER_COLUMNS = ["USER_CLASS", "USER_ANNOT", "USER_COM"]
```

Replace the `load_variants` signature and body (lines 20-51) with:

```python
def load_variants(
    xlsx_path: Path,
    genes: list[str] | None = None,
    column_aliases: dict[str, str] | None = None,
    synthesize_user_columns: bool = False,
    required_columns: list[str] | None = None,
) -> pd.DataFrame:
    """Load the 'Variants' sheet and normalise it to the internal vocabulary.

    NiourK format: per-sample columns are named '<sample_id>\\n<column_name>'.
    After stripping that prefix we (1) apply `column_aliases` to rename the
    raw schema columns onto the internal names, (2) optionally synthesize the
    USER_* columns when the schema lacks them, then (3) validate that the
    columns the rest of the pipeline depends on are present.
    """
    df = pd.read_excel(xlsx_path, sheet_name="Variants", engine="openpyxl")

    # Strip sample prefix: columns containing '\n'
    renamed = {}
    for col in df.columns:
        if "\n" in str(col):
            renamed[col] = str(col).split("\n", 1)[1]
    df = df.rename(columns=renamed)

    # Map raw schema columns onto the internal vocabulary
    if column_aliases:
        present = {raw: internal for raw, internal in column_aliases.items()
                   if raw in df.columns}
        df = df.rename(columns=present)

    # Fabricate USER_* columns the 069 Variants sheet doesn't carry
    if synthesize_user_columns:
        for col in _USER_COLUMNS:
            if col not in df.columns:
                df[col] = pd.NA

    # Validate the schema profile matches this file
    for col in (required_columns or []):
        if col not in df.columns:
            raise SchemaMismatchError(
                f"Erreur : colonne attendue '{col}' absente du fichier — "
                "mauvais profil de schéma ? Essayez --legacy."
            )

    # Filter to target genes
    if genes:
        df = df[df["SYMBOL"].isin(genes)].reset_index(drop=True)

    # Deduplicate by REF_TRANS: for each position, keep only the canonical
    # transcript (REF_TRANS == "VRAI") when one exists; keep all rows otherwise.
    pos_key = ["CHR", "POS", "REF", "ALT"]
    if "REF_TRANS" in df.columns and all(c in df.columns for c in pos_key):
        is_ref = df["REF_TRANS"].astype(str).str.upper().isin(["VRAI", "TRUE", "1"])
        has_ref = df.groupby(pos_key)["REF_TRANS"].transform(
            lambda s: s.astype(str).str.upper().isin(["VRAI", "TRUE", "1"]).any()
        ).astype(bool)
        df = df[~has_ref | is_ref].reset_index(drop=True)

    return df
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `python -m pytest tests/test_dual_schema.py -q`
Expected: PASS (4 passed).

- [ ] **Step 6: Run the existing suite to confirm no regression**

Run: `python -m pytest tests/test_load.py -q`
Expected: PASS — existing calls use only `genes=`, so new kwargs default to no-ops.

- [ ] **Step 7: Commit**

```bash
git add vcf_filter.py tests/conftest.py tests/test_dual_schema.py
git commit -m "feat: schema aliasing, USER_* synthesis and validation in load_variants"
```

---

## Task 2: `af_filter_column` in the frequency filters

**Files:**
- Modify: `vcf_filter.py:69-123` (`apply_exclusions`), `vcf_filter.py:130-190` (`classify_sheets`)
- Test: `tests/test_dual_schema.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_dual_schema.py`:

```python
from vcf_filter import apply_exclusions, classify_sheets


def _base_cfg():
    return {
        "exclusions": {
            "gnomad_max_af": 0.005,
            "rncnt_max_total": 6,
            "pjcnt_max_total": 10,
            "exclude_hom": True,
            "exclude_clinvar_benign": True,
            "mosaic_artifact": {"enabled": False},
        },
        "af_filter_column": "GNOMAD4_AF",
    }


def test_exclusions_use_af_filter_column(make_xlsx_069):
    import pandas as pd
    # Variant common in GNOMAD4_AF but rare in popmax → must be excluded
    df = pd.DataFrame([
        {"GNOMAD4_AF": 0.02, "MAX_AF_GNOMADEX2.1": 0.0001,
         "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0, "RNCNT_HOM": 0,
         "RNCNT_LOW": 0, "AF_PCT": 48.0, "CLINVAR": None},
        {"GNOMAD4_AF": 0.0001, "MAX_AF_GNOMADEX2.1": 0.0001,
         "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0, "RNCNT_HOM": 0,
         "RNCNT_LOW": 0, "AF_PCT": 48.0, "CLINVAR": None},
    ])
    out, report = apply_exclusions(df, _base_cfg())
    assert len(out) == 1
    assert report["gnomad"] == 1


def test_classify_sheets_use_af_filter_column():
    import pandas as pd
    cfg = {
        "af_filter_column": "GNOMAD4_AF",
        "sheets": {"fq_moderate": {"enabled": True, "gnomad_max_af": 0.001}},
    }
    df = pd.DataFrame([
        # rare in GNOMAD4_AF + MODERATE → in FQ_MODERATE
        {"GNOMAD4_AF": 0.0001, "MAX_AF_GNOMADEX2.1": 0.5, "IMPACT": "MODERATE",
         "CLINVAR": None},
        # common in GNOMAD4_AF → excluded from FQ_MODERATE
        {"GNOMAD4_AF": 0.5, "MAX_AF_GNOMADEX2.1": 0.0001, "IMPACT": "MODERATE",
         "CLINVAR": None},
    ])
    sheets = classify_sheets(df, cfg)
    assert len(sheets["FQ_MODERATE"]) == 1
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_dual_schema.py -k af_filter_column -q`
Expected: FAIL — `KeyError: 'GNOMAD4_AF'` (code still reads the literal `MAX_AF_GNOMADEX2.1`).

- [ ] **Step 3: Implement — `apply_exclusions`**

In `apply_exclusions`, after `exc = cfg["exclusions"]` (line 76) add:

```python
    af_col = cfg.get("af_filter_column", "MAX_AF_GNOMADEX2.1")
```

Change line 81 from:

```python
    mask = _safe_gt(df["MAX_AF_GNOMADEX2.1"], exc["gnomad_max_af"], fill=0.0)
```

to:

```python
    mask = _safe_gt(df[af_col], exc["gnomad_max_af"], fill=0.0)
```

- [ ] **Step 4: Implement — `classify_sheets`**

In `classify_sheets`, after `scfg = cfg["sheets"]` (line 136) add:

```python
    af_col = cfg.get("af_filter_column", "MAX_AF_GNOMADEX2.1")
```

Change each of the three lines 171, 178, 185 from:

```python
        gnomad_ok = df["MAX_AF_GNOMADEX2.1"].fillna(0) < fq_cfg["gnomad_max_af"]
```

to:

```python
        gnomad_ok = df[af_col].fillna(0) < fq_cfg["gnomad_max_af"]
```

(`replace_all` is safe here — all three lines are identical.)

- [ ] **Step 5: Run to verify pass**

Run: `python -m pytest tests/test_dual_schema.py -k af_filter_column -q`
Expected: PASS (2 passed).

- [ ] **Step 6: Run exclusions + tiers suites for regression**

Run: `python -m pytest tests/test_exclusions.py tests/test_tiers.py -q`
Expected: PASS — these configs have no `af_filter_column` key, so `af_col` defaults to `MAX_AF_GNOMADEX2.1`, unchanged.

- [ ] **Step 7: Commit**

```bash
git add vcf_filter.py tests/test_dual_schema.py
git commit -m "feat: drive frequency filtering from af_filter_column config key"
```

---

## Task 3: `af_column` in rarity scoring and AF cell-coloring

**Files:**
- Modify: `vcf_filter.py:327` (`_cell_fill`), `:368`, `:1251` (`_apply_cell_colors`)
- Modify: `vcf_filter.py:532` (`_score_variant`), `:738`, `:900` (`_apply_variant_scores`)
- Modify: `vcf_filter.py:1314` (`write_output`)
- Test: `tests/test_dual_schema.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_dual_schema.py`:

```python
from vcf_filter import _cell_fill, _score_variant


def test_cell_fill_honours_af_column():
    # With af_column=GNOMAD4_AF, that column gets the AF coloring...
    assert _cell_fill("GNOMAD4_AF", None, af_column="GNOMAD4_AF") is not None
    # ...and the default popmax column no longer matches.
    assert _cell_fill("MAX_AF_GNOMADEX2.1", 0.0, af_column="GNOMAD4_AF") is None


def test_score_variant_reads_af_column():
    row = {
        "SB": 0.50, "CALLNB": 4, "RNCNT_TOTAL": 0, "PJCNT_TOTAL": 0,
        "DP": 200, "AF_PCT": 48.0, "CLINVAR": "Pathogenic", "IMPACT": "HIGH",
        "Consequence": "stop_gained", "CADD_PHRED": 35.0, "REVEL": 0.9,
        "GNOMAD4_AF": None, "MAX_AF_GNOMADEX2.1": 0.5,  # popmax common, v4 novel
    }
    *_, comment = _score_variant(row, af_column="GNOMAD4_AF")
    assert "Absent gnomAD (novel)" in comment
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_dual_schema.py -k "af_column" -q`
Expected: FAIL — `_cell_fill()`/`_score_variant()` got an unexpected keyword argument `af_column`.

- [ ] **Step 3: Implement — `_cell_fill`**

Change the signature at line 327 from:

```python
def _cell_fill(col: str, value) -> "PatternFill | None":
```

to:

```python
def _cell_fill(col: str, value, af_column: str = "MAX_AF_GNOMADEX2.1") -> "PatternFill | None":
```

Change line 368 from:

```python
    elif col == "MAX_AF_GNOMADEX2.1":
```

to:

```python
    elif col == af_column:
```

- [ ] **Step 4: Implement — `_apply_cell_colors`**

Change the signature at line 1251 from:

```python
def _apply_cell_colors(ws) -> None:
```

to:

```python
def _apply_cell_colors(ws, af_column: str = "MAX_AF_GNOMADEX2.1") -> None:
```

Change line 1259 from:

```python
            fill = _cell_fill(col_name, cell.value)
```

to:

```python
            fill = _cell_fill(col_name, cell.value, af_column=af_column)
```

- [ ] **Step 5: Implement — `_score_variant`**

Change the signature at line 532 from:

```python
def _score_variant(row: dict, scoring: dict | None = None) -> tuple[float, float, float, float, str]:
```

to:

```python
def _score_variant(row: dict, scoring: dict | None = None, af_column: str = "MAX_AF_GNOMADEX2.1") -> tuple[float, float, float, float, str]:
```

Change line 738 from:

```python
    gnomad = _num(row.get("MAX_AF_GNOMADEX2.1"))
```

to:

```python
    gnomad = _num(row.get(af_column))
```

- [ ] **Step 6: Implement — `_apply_variant_scores`**

Change the signature at line 900 from:

```python
def _apply_variant_scores(ws, scoring: dict | None = None) -> None:
```

to:

```python
def _apply_variant_scores(ws, scoring: dict | None = None, af_column: str = "MAX_AF_GNOMADEX2.1") -> None:
```

Find the `_score_variant(` call inside this function and pass `af_column=af_column`. Read the function body first to locate the exact call:

Run: `grep -n "_score_variant(" vcf_filter.py`

For the call site inside `_apply_variant_scores` (between lines 900 and 928), add the kwarg, e.g. change `_score_variant(row_dict, scoring)` to `_score_variant(row_dict, scoring, af_column=af_column)` (match the actual argument names present).

- [ ] **Step 7: Implement — `write_output` threads `af_column`**

Change the `write_output` signature (lines 1314-1322) to add the parameter after `output_columns`:

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
) -> None:
```

In the data-sheet loop (lines 1342-1343) change:

```python
            _apply_cell_colors(writer.sheets[name])
            _apply_variant_scores(writer.sheets[name], scoring)
```

to:

```python
            _apply_cell_colors(writer.sheets[name], af_column=af_column)
            _apply_variant_scores(writer.sheets[name], scoring, af_column=af_column)
```

- [ ] **Step 8: Run to verify pass**

Run: `python -m pytest tests/test_dual_schema.py -k "af_column" -q`
Expected: PASS (2 passed).

- [ ] **Step 9: Run the full CLI suite for regression**

Run: `python -m pytest tests/test_cli.py -q`
Expected: PASS — every new param defaults to `MAX_AF_GNOMADEX2.1`, so 063 behaviour is unchanged.

- [ ] **Step 10: Commit**

```bash
git add vcf_filter.py tests/test_dual_schema.py
git commit -m "feat: thread af_column through rarity scoring and AF cell-coloring"
```

---

## Task 4: `main()` — `--legacy`, config resolution, clean mismatch exit

**Files:**
- Modify: `vcf_filter.py:1357-1408` (`main`)
- Test: `tests/test_dual_schema.py`

- [ ] **Step 1: Write the failing integration + error-path tests**

Append to `tests/test_dual_schema.py`:

```python
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent


def _run(args):
    return subprocess.run(
        [sys.executable, "vcf_filter.py", *args],
        capture_output=True, text=True, cwd=str(REPO),
    )


def test_cli_069_default_config_produces_output(make_xlsx_069):
    rows = [
        {"IMPACT": "HIGH", "Consequence": "stop_gained",
         "GNOMAD4_AF": 0.0, "MAX_AF_GNOMADEX2.1": 0.0,
         "CHR": "chr9", "POS": 200, "REF": "C", "ALT": "T"},
    ]
    path = make_xlsx_069(rows)
    result = _run([str(path)])  # no --config, no --legacy → 069 default
    assert result.returncode == 0, f"STDERR: {result.stderr}"
    out = path.with_name(path.stem + "_filtered.xlsx")
    assert out.exists()
    sheets = pd.ExcelFile(out, engine="openpyxl").sheet_names
    assert "HIGH_IMPACT" in sheets


def test_cli_063_legacy_flag_still_works(make_xlsx):
    rows = [
        {"IMPACT": "HIGH", "Consequence": "stop_gained",
         "CHR": "chr9", "POS": 200, "REF": "C", "ALT": "T"},
    ]
    path = make_xlsx(rows)
    result = _run(["--legacy", str(path)])
    assert result.returncode == 0, f"STDERR: {result.stderr}"
    out = path.with_name(path.stem + "_filtered.xlsx")
    assert out.exists()


def test_cli_063_under_default_config_exits_clean(make_xlsx):
    # 063 file lacks GNOMAD4_AF → default (069) config must fail cleanly.
    path = make_xlsx([{"SYMBOL": "TSC1"}])
    result = _run([str(path)])
    assert result.returncode != 0
    assert "GNOMAD4_AF" in (result.stdout + result.stderr)
    assert "--legacy" in (result.stdout + result.stderr)
    assert "Traceback" not in result.stderr


def test_cli_069_under_legacy_exits_clean(make_xlsx_069):
    # 069 file under legacy config: GENE never aliased → SYMBOL missing.
    path = make_xlsx_069()
    result = _run(["--legacy", str(path)])
    assert result.returncode != 0
    assert "SYMBOL" in (result.stdout + result.stderr)
    assert "Traceback" not in result.stderr
```

- [ ] **Step 2: Run to verify failure**

Run: `python -m pytest tests/test_dual_schema.py -k cli -q`
Expected: FAIL — `config_legacy.yaml` does not exist yet and/or `--legacy` is an unrecognized argument. (These pass after Task 5; this task wires the CLI, Task 5 adds the configs. Expect them red until Task 5 is done — that is acceptable; re-run at the end of Task 5.)

- [ ] **Step 3: Implement — replace `main()` argument parsing + load**

Replace lines 1358-1374 (from the `parser = ...` line through the `df = load_variants(...)` call) with:

```python
    parser = argparse.ArgumentParser(description="KODEX — Déchiffre les variants NiourK (TSC1/TSC2).")
    parser.add_argument("xlsx", type=Path, help="Input NiourK XLSX file")
    parser.add_argument("--config", type=Path, default=None,
                        help="YAML config file (overrides --legacy)")
    parser.add_argument("--legacy", action="store_true",
                        help="Use the run 1 schema profile (config_legacy.yaml)")
    parser.add_argument("--beta", action="store_true",
                        help="Skip pre-exclusions — apply sheet filters on raw variants only")
    args = parser.parse_args(argv)

    if not args.xlsx.is_file():
        sys.exit(
            f"Erreur : fichier introuvable : '{args.xlsx}'\n"
            "Vérifiez le chemin (un motif comme '*01.xlsx' qui ne correspond à aucun "
            "fichier est transmis tel quel par le shell)."
        )

    here = Path(__file__).parent
    if args.config is not None:
        config_path = args.config
    elif args.legacy:
        config_path = here / "config_legacy.yaml"
    else:
        config_path = here / "config.yaml"

    cfg = load_config(config_path)
    af_col = cfg.get("af_filter_column", "MAX_AF_GNOMADEX2.1")
    aliases = cfg.get("column_aliases")
    synth = cfg.get("synthesize_user_columns", False)

    # Load + gene filter + REF_TRANS dedup
    try:
        df = load_variants(
            args.xlsx, genes=cfg.get("genes"),
            column_aliases=aliases, synthesize_user_columns=synth,
            required_columns=["SYMBOL", af_col],
        )
    except SchemaMismatchError as e:
        sys.exit(str(e))
```

- [ ] **Step 4: Implement — update the remaining `main()` body**

The rest of `main()` still references `args.config`. Replace each remaining `args.config` with `config_path` (three sites: the `build_pre_report(... str(args.config) ...)` call ~line 1399, and the `write_output(... str(args.config) ...)` call ~line 1404). Update the extra-genes load (~line 1392) to forward the new kwargs, and add `af_column=af_col` to the `write_output` call. The end of `main()` becomes:

```python
    # Feuille fourre-tout pour les gènes extra (RHEB, PKD1, …) — aucun filtre qualité
    extra_genes = cfg.get("sheets", {}).get("extra_genes", [])
    if extra_genes:
        df_extra = load_variants(
            args.xlsx, genes=extra_genes,
            column_aliases=aliases, synthesize_user_columns=synth,
            required_columns=["SYMBOL", af_col],
        )
        if not df_extra.empty:
            sheet_name = "_".join(extra_genes)
            sheets[sheet_name] = df_extra

    # Pre-report
    input_name = args.xlsx.stem
    report = build_pre_report(total_rows, unique_variants, exc_report, sheets, str(config_path), input_name)
    print(report)

    # Write output
    output_path = args.xlsx.with_name(f"{args.xlsx.stem}_filtered.xlsx")
    write_output(output_path, sheets, report, str(config_path), cfg.get("output_columns", []),
                 scoring=cfg.get("scoring"),
                 discordance=cfg.get("af_discordance"),
                 af_column=af_col)
    print(f"\n→ Excel écrit : {output_path}")
```

- [ ] **Step 5: Defer test run to Task 5**

The CLI tests need `config_legacy.yaml` and the updated `config.yaml`. Proceed to Task 5, then run them.

- [ ] **Step 6: Commit**

```bash
git add vcf_filter.py tests/test_dual_schema.py
git commit -m "feat: add --legacy flag, config-path resolution and clean schema-mismatch exit"
```

---

## Task 5: Config files — 069 default + 063 legacy

**Files:**
- Create: `config_legacy.yaml`
- Modify: `config.yaml`

- [ ] **Step 1: Create `config_legacy.yaml` from the current `config.yaml`**

```bash
cp config.yaml config_legacy.yaml
```

- [ ] **Step 2: Add the three legacy-valued keys to `config_legacy.yaml`**

Insert immediately after the `genes: [TSC1, TSC2]` line near the top:

```yaml
# --- Schema profile (run 1 / legacy) ---
column_aliases: {}                       # 063 columns already match internal names
af_filter_column: MAX_AF_GNOMADEX2.1     # popmax gnomAD ex2.1
synthesize_user_columns: false           # 063 Variants sheet already has USER_*
```

(`output_columns` in `config_legacy.yaml` stays exactly as the current file — no change.)

- [ ] **Step 3: Add the three 069 keys to `config.yaml`**

Insert immediately after the `genes: [TSC1, TSC2]` line:

```yaml
# --- Schema profile (run 2 / default) ---
column_aliases:                          # raw 069 column -> internal name
  GENE: SYMBOL
  CSQ: Consequence
  REGION: REGIONS
af_filter_column: GNOMAD4_AF             # gnomAD v4 overall AF (see spec Open question)
synthesize_user_columns: true            # 069 USER_* live in Valideur*/Technicien sheets
```

- [ ] **Step 4: Rewrite `output_columns` in `config.yaml` for 069**

Replace the entire `output_columns:` block (current lines 192-273, from `output_columns:` to end of file) with:

```yaml
output_columns:
  # --- Identifiant & localisation ---
  - Import DefGen
  - SYMBOL
  - REF_TRANS
  - REGIONS
  - HGVSc
  - HGVSp
  # --- Données échantillon ---
  - AF_PCT
  - DP
  - SB
  # --- Annotations utilisateur (synthétisées vides en 069) ---
  - USER_CLASS
  - USER_ANNOT
  - USER_COM
  # --- Callers ---
  - CALLNB
  # --- Impact fonctionnel ---
  - IMPACT
  - Consequence
  - CLINVAR
  # --- Compteurs réseau ---
  - RNCNT_TOTAL
  - PJCNT_TOTAL
  - RNCNT_LOW
  - PJCNT_LOW
  - RNCNT_HTZ
  - PJCNT_HTZ
  - RNCNT_HOM
  - PJCNT_HOM
  # --- Bases de données ---
  - LOVD
  - VARSOME
  - FILTER
  # --- Fréquences populationnelles (gnomAD v4 + popmax ex2.1) ---
  - MAX_AF_GNOMADEX2.1
  - GNOMAD4_AF
  - GNOMAD4_AC
  - GNOMAD4_AN
  - GNOMAD4_nhomalt
  # --- Scores prédictifs ---
  - SPiP_Interpretation
  - SPiP_InterConfident
  - SpliceAI
  - SpliceAI_num
  - CADD_PHRED
  - REVEL
  - SIFT
  - PolyPhen
  # --- Références ---
  - Existing_variation
  - PUBMED
  # --- Coordonnées génomiques ---
  - CHR
  - POS
  - REF
  - ALT
  - Feature
  # --- Détail callers ---
  - "CALLAF_mutect2|strelka|gatkHC|deepvariant"
  - CALLAD_mutect2
  - CALLAD_strelka
  - CALLAD_gatkHC
  - CALLAD_deepvariant
  - CALLFILTER_mutect2
  - CALLFILTER_strelka
  - CALLFILTER_gatkHC
  - CALLFILTER_deepvariant
```

- [ ] **Step 5: Sanity-check both YAML files parse**

Run: `python -c "import yaml; yaml.safe_load(open('config.yaml')); yaml.safe_load(open('config_legacy.yaml')); print('OK')"`
Expected: `OK`.

- [ ] **Step 6: Run the deferred CLI tests from Task 4**

Run: `python -m pytest tests/test_dual_schema.py -k cli -q`
Expected: PASS (4 passed) — 069 default works, 063 `--legacy` works, both mismatch paths exit cleanly.

- [ ] **Step 7: Run the full test suite**

Run: `python -m pytest -q`
Expected: PASS — existing 063 tests use `--config config.yaml` explicitly... **verify this assumption** (see note below) before relying on it.

> **Important regression note:** `tests/test_cli.py` passes `--config <repo>/config.yaml` explicitly. After this task `config.yaml` is the **069** profile, so those 063-format `make_xlsx` fixtures will now hit a schema mismatch (missing `GNOMAD4_AF`). If `python -m pytest tests/test_cli.py -q` fails for that reason, update those tests to use `--config config_legacy.yaml` (they exercise 063-format fixtures). Make that edit, then re-run. Do **not** weaken the schema validation to make them pass.

- [ ] **Step 8: Commit**

```bash
git add config.yaml config_legacy.yaml tests/test_cli.py
git commit -m "feat: 069 default config + config_legacy.yaml for 063 schema"
```

---

## Task 6: Manual end-to-end verification on real files

**Files:** none (verification only)

- [ ] **Step 1: Run KODEX on a real 069 file (default profile)**

```bash
python vcf_filter.py <run_dir>/xlsx/run_2_Singularity_100005_S1.xlsx
```
Expected: pre-report prints, `…_filtered.xlsx` written next to the input, no traceback.

- [ ] **Step 2: Run KODEX on a real 063 file with `--legacy`**

```bash
python vcf_filter.py --legacy <run_dir>/xlsx/run_1_100003_S7.xlsx
```
Expected: pre-report prints, `…_filtered.xlsx` written, no traceback.

- [ ] **Step 3: Confirm the mismatch guidance on a real 063 file under default profile**

```bash
python vcf_filter.py <run_dir>/xlsx/run_1_100003_S7.xlsx
```
Expected: clean one-line error naming `GNOMAD4_AF` and suggesting `--legacy`; exit code non-zero; no Python traceback.

- [ ] **Step 4: Open both output workbooks and spot-check**

Confirm in the 069 output: `SYMBOL`/`Consequence`/`REGIONS` headers present (aliased), `USER_CLASS/ANNOT/COM` columns present but empty, `GNOMAD4_AF` column coloured green/yellow, scoring comments populated. Confirm the 063 `--legacy` output is visually identical to pre-change runs.

---

## Self-Review

**Spec coverage:**
- Schema diff / aliasing → Task 1 (`column_aliases`). ✓
- USER_* synthesis → Task 1 (`synthesize_user_columns`). ✓
- `af_filter_column` replacing all 6 references → Task 2 (filters `:81,:171,:178,:185`) + Task 3 (scoring `:738`, coloring `:368`). ✓
- `--legacy` shorthand, explicit `--config` wins, no autodetection → Task 4. ✓
- Clean schema-mismatch error (no traceback) → Task 1 (`SchemaMismatchError`) + Task 4 (`sys.exit`). ✓
- `config.yaml` 069 defaults + `config_legacy.yaml` → Task 5. ✓
- `output_columns` for 069 → Task 5 Step 4. ✓
- Testing plan (unit load, unit legacy unchanged, integration both schemas, error paths, fixtures) → Tasks 1-5 tests + Task 6 manual. ✓
- gnomAD threshold recalibration → explicitly out of scope (spec Open question); thresholds kept as-is. ✓

**Placeholder scan:** No TBD/TODO; every code step shows full code. One deliberate dependency: Task 4's CLI tests are written then run after Task 5 (configs exist) — flagged explicitly in both tasks, not a placeholder.

**Type consistency:** `af_column` (str) threaded identically through `_cell_fill`, `_apply_cell_colors`, `_score_variant`, `_apply_variant_scores`, `write_output`; default `"MAX_AF_GNOMADEX2.1"` everywhere. `load_variants` kwargs (`column_aliases: dict`, `synthesize_user_columns: bool`, `required_columns: list[str]`) match their call sites in `main()`. `SchemaMismatchError` raised in `load_variants`, imported in tests, caught in `main()`. Config keys (`column_aliases`, `af_filter_column`, `synthesize_user_columns`) match between YAML and `cfg.get(...)` reads.
