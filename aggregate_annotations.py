"""aggregate_annotations.py — KODEX : agrège les variants annotés de tous les *_filtered.xlsx."""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vcf_filter import _format_sheet, _apply_cell_colors, _fix_lovd_url
from igv_capture import extract_sample_id

_ANNOTATION_COLS = ["USER_CLASS", "USER_ANNOT", "USER_COM"]
_DEDUP_KEY = ["CHR", "POS", "REF", "ALT"]

# Colonnes affichées en priorité dans le tableau de synthèse
_PRIORITY_COLS = [
    "SAMPLE", "SOURCE_SHEET",
    "USER_CLASS", "USER_ANNOT", "USER_COM",
    "SYMBOL", "HGVSc", "HGVSp",
    "AF_PCT", "DP", "SB", "CALLNB",
    "IMPACT", "Consequence",
    "CLINVAR", "LOVD",
    "SpliceAI_num", "CADD_PHRED", "REVEL",
    "SPiP_Interpretation", "SPiP_InterConfident",
    "RNCNT_TOTAL", "PJCNT_TOTAL",
    "RNCNT_LOW", "PJCNT_LOW",
    "MAX_AF_GNOMADEX2.1",
    "CHR", "POS", "REF", "ALT",
    "Existing_variation", "PUBMED",
]


def _collect_annotated(xlsx_path: Path) -> list[dict]:
    """Retourne les variants annotés d'un fichier filtré, avec la colonne SOURCE_SHEET.

    - Parcourt toutes les feuilles sauf _summary
    - Garde les lignes ayant au moins un champ USER_* non vide
    - Déduplique par CHR+POS+REF+ALT : si un variant apparaît dans plusieurs feuilles,
      les noms de feuilles sont concaténés dans SOURCE_SHEET
    """
    xl = pd.ExcelFile(xlsx_path, engine="openpyxl")
    data_sheets = [s for s in xl.sheet_names if s != "_summary"]

    seen: dict[tuple, dict] = {}  # (chr, pos, ref, alt) → {data, sheets}

    for sheet in data_sheets:
        df = xl.parse(sheet)
        mask = pd.Series(False, index=df.index)
        for col in _ANNOTATION_COLS:
            if col in df.columns:
                vals = df[col].astype(str).str.strip()
                mask |= df[col].notna() & (vals != "") & (~vals.str.startswith("="))

        for _, row in df[mask].iterrows():
            key = tuple(str(row.get(k, "")) for k in _DEDUP_KEY)
            if key not in seen:
                seen[key] = {"data": row.to_dict(), "sheets": [sheet]}
            elif sheet not in seen[key]["sheets"]:
                seen[key]["sheets"].append(sheet)

    result = []
    for info in seen.values():
        d = info["data"]
        d["SOURCE_SHEET"] = " | ".join(info["sheets"])
        result.append(d)

    return result


def aggregate_all(xlsx_dir: Path, files: list[Path] | None = None) -> pd.DataFrame:
    """Agrège les variants annotés de tous les *_filtered.xlsx (ou *.xlsx) du répertoire."""
    if files is None:
        files = sorted(xlsx_dir.glob("*_filtered.xlsx"))
        if not files:
            # Fallback : tous les xlsx sauf le fichier de sortie et les lock files
            files = sorted(
                f for f in xlsx_dir.glob("*.xlsx")
                if not f.name.startswith("~$")
                and "annotations" not in f.name.lower()
            )
    if not files:
        print(f"Aucun fichier xlsx trouvé dans {xlsx_dir}")
        return pd.DataFrame()

    all_rows: list[dict] = []

    for xlsx_path in files:
        try:
            sample_id = extract_sample_id(xlsx_path)
        except Exception:
            sample_id = xlsx_path.stem.replace("_filtered", "")

        try:
            rows = _collect_annotated(xlsx_path)
        except Exception as e:
            print(f"  {sample_id:<18s} [SKIP] {e}")
            continue

        if rows:
            try:
                _write_annotation_sheet(rows, xlsx_path)
            except Exception as e:
                print(f"  {sample_id:<18s} [WARN feuille locale] {e}")

        for row in rows:
            row["SAMPLE"] = sample_id
        all_rows.extend(rows)

        status = f"{len(rows)} variant(s)" if rows else "—"
        print(f"  {sample_id:<18s} {status}")

    if not all_rows:
        return pd.DataFrame()

    df = pd.DataFrame(all_rows)

    if "POS" in df.columns:
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce")

    df = df.sort_values(
        ["SAMPLE", "CHR", "POS"], na_position="last"
    ).reset_index(drop=True)

    return df


def _write_annotation_sheet(rows: list[dict], xlsx_path: Path) -> None:
    """Écrit (ou remplace) la feuille _annotations dans un *_filtered.xlsx existant."""
    # Colonnes prioritaires sans SAMPLE (redondant dans un fichier mono-échantillon)
    priority = [c for c in _PRIORITY_COLS if c != "SAMPLE"]

    df = pd.DataFrame(rows)
    if "POS" in df.columns:
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce")

    remaining = [c for c in df.columns if c not in priority]
    ordered = [c for c in priority if c in df.columns] + remaining
    df_out = _fix_lovd_url(df[ordered])

    with pd.ExcelWriter(
        xlsx_path, engine="openpyxl", mode="a", if_sheet_exists="replace"
    ) as writer:
        df_out.to_excel(writer, sheet_name="_annotations", index=False)
        ws = writer.sheets["_annotations"]
        _format_sheet(ws)
        _apply_cell_colors(ws)


def _write_summary(df: pd.DataFrame, output_path: Path) -> None:
    """Écrit le tableau de synthèse dans un XLSX formaté."""
    remaining = [c for c in df.columns if c not in _PRIORITY_COLS]
    ordered = [c for c in _PRIORITY_COLS if c in df.columns] + remaining

    df_out = _fix_lovd_url(df[ordered])

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df_out.to_excel(writer, sheet_name="ANNOTATIONS", index=False)
        ws = writer.sheets["ANNOTATIONS"]
        _format_sheet(ws)
        _apply_cell_colors(ws)

    print(f"\n→ {output_path}")
    print(f"  {len(df)} variant(s) annotés — {df['SAMPLE'].nunique()} sample(s)")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="KODEX — Agrège les variants annotés (USER_CLASS/ANNOT/COM) de tous les *_filtered.xlsx."
    )
    parser.add_argument(
        "targets", type=Path, nargs="*",
        help="Répertoire OU fichier(s) xlsx (ex: TSC_063*_filtered.xlsx). "
             "Défaut: ~/data/TSC_063/xlsx",
    )
    parser.add_argument(
        "--out", type=Path,
        help="Fichier de sortie (défaut: <répertoire>/TSC_063_annotations.xlsx)",
    )
    args = parser.parse_args(argv)

    targets = [t.expanduser().resolve() for t in args.targets] if args.targets else []

    # Déterminer la liste de fichiers et le répertoire de sortie
    if not targets:
        xlsx_dir = Path("~/data/TSC_063/xlsx").expanduser().resolve()
        files = None  # aggregate_all fera le glob
    elif len(targets) == 1 and targets[0].is_dir():
        xlsx_dir = targets[0]
        files = None
    else:
        # Un ou plusieurs fichiers passés directement
        files = [f for f in targets if f.suffix == ".xlsx" and f.exists()]
        invalid = [f for f in targets if not f.exists()]
        if invalid:
            for f in invalid:
                print(f"[WARN] Fichier introuvable : {f}")
        if not files:
            print("[ERREUR] Aucun fichier xlsx valide trouvé.")
            raise SystemExit(1)
        xlsx_dir = files[0].parent

    out_path = (args.out or xlsx_dir / "TSC_063_annotations.xlsx").expanduser().resolve()

    if files is not None:
        print(f"Lecture de {len(files)} fichier(s) ...")
        df = aggregate_all(xlsx_dir, files=files)
    else:
        print(f"Lecture de {xlsx_dir} ...")
        df = aggregate_all(xlsx_dir)

    if df.empty:
        print("Aucun variant annoté trouvé.")
        raise SystemExit(0)

    _write_summary(df, out_path)


if __name__ == "__main__":
    main()
