# KODEX — Dual-schema support (NiourK 063 ↔ 069) via `--legacy`

**Date:** 2026-06-15
**Status:** Approved (design), pending implementation plan

## Problem

KODEX (`vcf_filter.py`) was built against the NiourK run 1 XLSX schema. The newer
run 2 NiourK export uses a different, restructured schema. Running KODEX on a 069
file currently crashes with `KeyError: 'SYMBOL'`.

We want one tool that processes **both** schemas: the new (069) schema as the default,
and the old (063) schema via a `--legacy` flag. Schema differences should live in YAML
config, with minimal branching in code.

## Schema diff (063 → 069)

Confirmed by diffing raw inputs
`run 1/xlsx/run_1_100003_S7.xlsx` vs `run_2_xlsx/run_2_Singularity_100005_S1.xlsx`.

**Sheets:** 063 = `Variants, Config, Blanket`. 069 adds `Valideur1, Valideur2,
Technicien, AIDE`.

**Variants columns:** 68 (063) → 54 (069).

Column-level changes that matter to KODEX:

| Role | 063 name | 069 name | Used in code at |
|---|---|---|---|
| Gene symbol | `SYMBOL` | `GENE` | gene filter `:39`; LOVD URL `:1300`, `:1306` |
| Consequence | `Consequence` | `CSQ` | scoring `:644`; downstream-gene coloring `:1118-1145` |
| Region | `REGIONS` | `REGION` | output only (cosmetic) |
| Pop. AF (filter) | `MAX_AF_GNOMADEX2.1` (present in BOTH) | `GNOMAD4_AF` (069 only) | exclusions/scoring `:81,:171,:178,:185,:368,:738` |
| User annotation | `USER_CLASS/ANNOT/COM` | absent from Variants (moved to Valideur*/Technicien sheets) | coloring `:393`; cross-sheet formulas `:1016-1086` |

**069-only additions:** `GNOMAD4_AC/AN/AF/nhomalt`, `FILTER`, `VARSOME`.
**063-only (dropped in 069):** all `GNOMAD3_AF*` and gnomADex2.1 per-population AF columns.

Note: `MAX_AF_GNOMADEX2.1` exists in BOTH schemas; the `CALLAF_*` caller-order column
name varies between individual 069 samples.

## Decisions (from brainstorming)

1. **Schema selection:** separate config files, no autodetection.
   `config.yaml` = 069 default; `config_legacy.yaml` = 063.
2. **USER_* columns in 069:** create them empty, preserving the 063 coloring +
   cross-sheet-formula workflow for the analyst to fill in.
3. **gnomAD source for 069:** switch frequency filtering to `GNOMAD4_AF`
   (a deliberate scientific change — see Open question).

## Design (Approach 1: config-driven aliases, internal vocabulary unchanged)

The internal column vocabulary stays as today's names (the 063 names already used
throughout the code). The new schema is mapped **onto** those names at load time. This
keeps code churn small (~10 call sites) and the diff reviewable for a clinical tool.

### A. Schema selection / CLI
- `config.yaml` is the default config = **069 schema**.
- `--legacy` is shorthand for `--config config_legacy.yaml`.
- If both `--legacy` and an explicit `--config X` are passed, **explicit `--config` wins**.
- No autodetection. A schema/file mismatch (e.g. 063 file under default config, or 069
  file under `--legacy`) must fail with a **clean message**, not a traceback — reusing
  the `is_file()`/clean-exit pattern already added in `main()`. Message names the missing
  expected column and suggests checking `--legacy`, e.g.:
  `Erreur : colonne attendue 'SYMBOL' absente du fichier — mauvais profil de schéma ? Essayez --legacy.`

### B. New config keys (consumed in `load_variants` + frequency filters)
- `column_aliases:` — a `{raw_name: internal_name}` rename map applied in
  `load_variants` immediately after the sample-prefix strip.
  - `config.yaml`: `{GENE: SYMBOL, CSQ: Consequence, REGION: REGIONS}`
  - `config_legacy.yaml`: empty (063 already matches internal names)
- `af_filter_column:` — the column driving frequency filtering. Replaces the 6 hardcoded
  `MAX_AF_GNOMADEX2.1` references.
  - `config.yaml`: `GNOMAD4_AF`
  - `config_legacy.yaml`: `MAX_AF_GNOMADEX2.1`
- `synthesize_user_columns: true` — when `USER_CLASS/USER_ANNOT/USER_COM` are absent
  after aliasing, `load_variants` creates them as empty columns so coloring +
  cross-sheet formula injection keep working.
  - `config.yaml`: `true`; `config_legacy.yaml`: `false` (063 already has them)

### C. output_columns
- `config.yaml`: 069 column set — add `GNOMAD4_AC/AN/AF/nhomalt`, `FILTER`, `VARSOME`;
  drop the gnomADex2.1 per-population list and `GNOMAD3_*`. Internal names (`SYMBOL`,
  `Consequence`, `REGIONS`) are used because aliasing already renamed the raw columns.
- `config_legacy.yaml`: today's `output_columns` list verbatim.
- Output selection already self-filters to present columns (`:1339`), so any listed
  column missing from a given file is silently skipped — no crash.

### D. Code changes (summary)
- `load_variants(...)`: accept the schema config; after prefix-strip, (1) apply
  `column_aliases` rename, (2) if `synthesize_user_columns`, add any missing USER_* as
  empty columns. Validate that required internal columns (`SYMBOL`) are present; if not,
  raise the clean schema-mismatch error.
- Frequency-filter sites (`:81,:171,:178,:185,:368,:738`): read the column name from
  `af_filter_column` instead of the literal `MAX_AF_GNOMADEX2.1`.
- `main()`: add `--legacy`; resolve config path (`--config` > `--legacy` > default).
- New file `config_legacy.yaml`: copy of current `config.yaml` content + the three new
  keys set to legacy values.
- Current `config.yaml`: add the three new keys (069 values) + update `output_columns`.

## Testing
- Unit: `load_variants` on a 069 fixture yields internal names (`SYMBOL`, `Consequence`,
  `REGIONS`) and empty USER_* columns.
- Unit: `load_variants` on a 063 fixture under legacy config is unchanged.
- Integration: full run on one 069 file (`run_2_Singularity_100005_S1.xlsx`) produces a
  `_filtered.xlsx` without error; full run on one 063 file under `--legacy` still works.
- Error path: 063 file under default config (and 069 under `--legacy`) exits cleanly with
  the schema-mismatch message, no traceback.
- Fixtures: small sliced copies of one 063 and one 069 Variants sheet in `tests/`.

## Open question (defer to user, not blocking)
- `GNOMAD4_AF` is overall AF, not popmax; the `0.005` / `0.001` thresholds were
  calibrated against `MAX_AF_GNOMADEX2.1` popmax. The plan keeps the current threshold
  numbers as the starting 069 default; recalibration is a separate scientific task.

## Out of scope
- Batch / multi-file processing (single file per invocation, as today).
- Reading actual USER_* values from the 069 Valideur*/Technicien sheets (we create them
  empty instead).
- gnomAD threshold recalibration.
