"""aggregate_annotations.py — KODEX : agrège les variants annotés de tous les *_filtered.xlsx."""

import argparse
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vcf_filter import (
    _apply_cell_colors,
    _apply_variant_scores,
    _fix_lovd_url,
    _format_sheet,
    add_score_column,
    load_config,
)
from igv_capture import extract_sample_id
from igv_link import inject_aggregate_igv_columns, resolve_bam_dir

_ANNOTATION_COLS = ["USER_CLASS", "USER_ANNOT"]
_DEDUP_KEY = ["CHR", "POS", "REF", "ALT"]


def _run_prefix(xlsx_dir: Path) -> str:
    """Run identifier déduit des noms de fichiers du répertoire.

    Évite le hardcode : le fichier de synthèse hérite du préfixe <lettres>_<chiffres>
    du run réellement présent (run_2_..._S1.xlsx → run_2). Fallback "RUN".
    """
    import re
    for f in sorted(xlsx_dir.glob("*.xlsx")):
        m = re.match(r"([A-Za-z]+_\d+)_", f.name)
        if m:
            return m.group(1)
    return "RUN"

# Colonnes affichées en priorité dans le tableau de synthèse
_PRIORITY_COLS = [
    "SCORE",
    "RETENU",                 # mode --all : le variant portait-il une annotation
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


def _annotation_mask(df: pd.DataFrame) -> "pd.Series":
    """Lignes portant une annotation de relecture — le critère de « variant retenu ».

    Une ligne est retenue dès qu'un champ USER_* est rempli (hors formules `=`).
    Score convention : un "0" nu dans USER_CLASS = artefact confirmé → pas une
    trouvaille en soi. "0?" (ou toute valeur finissant par "?") = doute → gardé.
    Une ligne avec un USER_ANNOT rempli ressort via l'autre branche du OR, même
    si USER_CLASS vaut "0".
    """
    mask = pd.Series(False, index=df.index)
    for col in _ANNOTATION_COLS:
        if col not in df.columns:
            continue
        vals = df[col].astype(str).str.strip()
        # Une colonne USER_CLASS ne contenant que des chiffres est relue en
        # float64 : "0" revient "0.0" et la règle du zéro nu ne s'appliquait
        # jamais — les artefacts confirmés remontaient dans la synthèse.
        vals = vals.str.replace(r"^(\d+)\.0$", r"\1", regex=True)
        valid = df[col].notna() & (vals != "") & (~vals.str.startswith("="))
        if col == "USER_CLASS":
            valid &= vals != "0"
        mask |= valid
    return mask


def _collect_all(
    xlsx_path: Path, extra_genes: list[str] | None = None
) -> list[dict]:
    """Retourne *tous* les variants d'un fichier filtré, retenus ou non.

    Contrepartie de `_collect_annotated` pour le mode --all : sert à contrôler le
    scoring, en vérifiant où les variants retenus à la main tombent dans un
    classement de l'ensemble. La colonne RETENU ("oui" / "") marque ceux que
    `_collect_annotated` aurait gardés — sans elle il faudrait croiser deux
    fichiers à la main.

    Les gènes de `extra_genes` (RHEB, PKD1) sont exclus : leur feuille n'applique
    aucun filtre qualité, ils gonfleraient le classement sans être comparables.
    Ils restent dans le classeur par patient, pour ce qu'ils sont — un cas TSC
    publié avec variant RHEB, et les délétions contiguës TSC2-PKD1.

    Déduplication par CHR+POS+REF+ALT, noms de feuilles concaténés, comme
    `_collect_annotated`.
    """
    excluded = {g.upper() for g in (extra_genes or [])}
    xl = pd.ExcelFile(xlsx_path, engine="openpyxl")
    data_sheets = [s for s in xl.sheet_names if not s.startswith("_")]

    seen: dict[tuple, dict] = {}

    for sheet in data_sheets:
        df = xl.parse(sheet)
        if df.empty:
            continue
        retained = _annotation_mask(df)

        for idx, row in df.iterrows():
            if excluded and str(row.get("SYMBOL", "")).upper() in excluded:
                continue
            key = tuple(str(row.get(k, "")) for k in _DEDUP_KEY)
            if key not in seen:
                data = row.to_dict()
                data["RETENU"] = "oui" if retained.loc[idx] else ""
                seen[key] = {"data": data, "sheets": [sheet]}
            else:
                # Un variant vu dans plusieurs feuilles est retenu s'il l'est
                # dans l'une d'elles — l'annotation vit sur une feuille à la fois.
                if retained.loc[idx]:
                    seen[key]["data"]["RETENU"] = "oui"
                if sheet not in seen[key]["sheets"]:
                    seen[key]["sheets"].append(sheet)

    result = []
    for info in seen.values():
        d = info["data"]
        d["SOURCE_SHEET"] = " | ".join(info["sheets"])
        result.append(d)

    return result


def _collect_annotated(xlsx_path: Path) -> list[dict]:
    """Retourne les variants annotés d'un fichier filtré, avec la colonne SOURCE_SHEET.

    - Parcourt toutes les feuilles de données (ignore les feuilles méta _*)
    - Garde les lignes ayant au moins un champ USER_* non vide (hors formules =)
    - USER_CLASS "0" nu = artefact confirmé → ignoré ; "0?" / "<x>?" = doute,
      gardé ; un nombre ou un "." = sûr, gardé. Un USER_ANNOT rempli garde la
      ligne même si USER_CLASS vaut "0".
    - Déduplique par CHR+POS+REF+ALT : si un variant apparaît dans plusieurs feuilles,
      les noms de feuilles sont concaténés dans SOURCE_SHEET
    """
    xl = pd.ExcelFile(xlsx_path, engine="openpyxl")
    data_sheets = [s for s in xl.sheet_names if not s.startswith("_")]

    seen: dict[tuple, dict] = {}  # (chr, pos, ref, alt) → {data, sheets}

    for sheet in data_sheets:
        df = xl.parse(sheet)
        for _, row in df[_annotation_mask(df)].iterrows():
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


def aggregate_all(
    xlsx_dir: Path,
    files: list[Path] | None = None,
    scoring: dict | None = None,
    every_variant: bool = False,
    extra_genes: list[str] | None = None,
    af_column: str = "MAX_AF_GNOMADEX2.1",
) -> pd.DataFrame:
    """Agrège les variants de tous les *_filtered.xlsx (ou *.xlsx) du répertoire.

    Par défaut, seuls les variants annotés. Avec `every_variant`, tous les
    variants sauf `extra_genes`, classés par score sur l'ensemble du run —
    de quoi vérifier où les variants retenus tombent dans le classement.
    """
    if files is None:
        files = sorted(
            f for f in xlsx_dir.glob("*_filtered.xlsx")
            if not f.name.startswith("~$")
        )
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
            rows = (
                _collect_all(xlsx_path, extra_genes=extra_genes)
                if every_variant
                else _collect_annotated(xlsx_path)
            )
        except Exception as e:
            print(f"  {sample_id:<18s} [SKIP] {e}")
            continue

        # La feuille _annotations locale reste le reflet de la relecture : en
        # mode --all on ne la réécrit pas avec l'intégralité des variants.
        if rows and not every_variant:
            try:
                _write_annotation_sheet(rows, xlsx_path, scoring=scoring,
                                        af_column=af_column)
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

    # Le SCORE est reconstruit ici : les feuilles sources peuvent venir d'une
    # version antérieure du scoring, et sans lui _sort_by_score retomberait sur
    # l'ordre génomique — ce qui viderait --all de son objet.
    #
    # Reconstruit aussi quand la colonne est *incomplète* : un run mélange des
    # classeurs produits par plusieurs versions de KODEX, les uns portant SCORE
    # et les autres non. La colonne existe alors mais reste NaN pour la majorité
    # des lignes, que le tri renvoie en fin de fichier — enterrant précisément
    # les variants en cours de relecture.
    if "SCORE" not in df.columns or df["SCORE"].isna().any():
        df = df.drop(columns=["SCORE"], errors="ignore")
        df = add_score_column(df, scoring, af_column=af_column)

    return _sort_by_score(df, group_by_sample=not every_variant).reset_index(drop=True)


def _sort_by_score(df: pd.DataFrame, group_by_sample: bool = True) -> pd.DataFrame:
    """Rank variants best score first, patients grouped or not.

    Patients stay grouped by default because the '▶ BAM' link reloads IGV tracks
    on every change of sample. The --all ranking passes group_by_sample=False:
    it exists to see where the retained variants land across the whole run, so
    the best variant leads whoever it belongs to. Falls back to genomic order
    for files produced before the SCORE column existed.
    """
    keys = ["SAMPLE", "SCORE"] if group_by_sample else ["SCORE"]
    if "SCORE" in df.columns:
        ascending = [True, False] if group_by_sample else False
        return df.sort_values(keys, ascending=ascending, na_position="last")
    fallback = ["SAMPLE", "CHR", "POS"] if group_by_sample else ["CHR", "POS"]
    return df.sort_values(fallback, na_position="last")


def _write_annotation_sheet(
    rows: list[dict], xlsx_path: Path, scoring: dict | None = None,
    af_column: str = "MAX_AF_GNOMADEX2.1",
) -> None:
    """Écrit (ou remplace) la feuille _annotations dans un *_filtered.xlsx existant."""
    # Colonnes prioritaires sans SAMPLE (redondant dans un fichier mono-échantillon)
    priority = [c for c in _PRIORITY_COLS if c != "SAMPLE"]

    df = pd.DataFrame(rows)
    if "POS" in df.columns:
        df["POS"] = pd.to_numeric(df["POS"], errors="coerce")

    if "SCORE" in df.columns:
        df = df.sort_values("SCORE", ascending=False, na_position="last")

    remaining = [c for c in df.columns if c not in priority]
    ordered = [c for c in priority if c in df.columns] + remaining
    df_out = _fix_lovd_url(df[ordered])

    with pd.ExcelWriter(
        xlsx_path, engine="openpyxl", mode="a", if_sheet_exists="replace"
    ) as writer:
        df_out.to_excel(writer, sheet_name="_annotations", index=False)
        ws = writer.sheets["_annotations"]
        _format_sheet(ws)
        _apply_cell_colors(ws, af_column=af_column)
        _apply_variant_scores(ws, scoring, af_column=af_column)


def _write_summary(
    df: pd.DataFrame,
    output_path: Path,
    igv_cfg: dict | None = None,
    scoring: dict | None = None,
    af_column: str = "MAX_AF_GNOMADEX2.1",
) -> None:
    """Écrit le tableau de synthèse dans un XLSX formaté."""
    remaining = [c for c in df.columns if c not in _PRIORITY_COLS]
    ordered = [c for c in _PRIORITY_COLS if c in df.columns] + remaining

    df_out = _fix_lovd_url(df[ordered])

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df_out.to_excel(writer, sheet_name="ANNOTATIONS", index=False)
        ws = writer.sheets["ANNOTATIONS"]
        _format_sheet(ws)
        _apply_cell_colors(ws, af_column=af_column)
        _apply_variant_scores(ws, scoring, af_column=af_column)

        # IGV links — last, so the column inserts don't shift earlier lookups.
        if igv_cfg and igv_cfg.get("enabled"):
            bam_dir = resolve_bam_dir(igv_cfg.get("bam_dir", "../bam"),
                                      output_path.parent)
            n = inject_aggregate_igv_columns(
                ws,
                bam_dir=bam_dir,
                reference_lane=igv_cfg.get("reference_lane", "S2"),
                port=igv_cfg.get("port", 60151),
                genome=igv_cfg.get("genome", "hg38"),
                half_window=igv_cfg.get("half_window", 145),
            )
            if n == 0:
                print(f"[IGV] Liens omis — BAMs introuvables dans {bam_dir}")

    label = "variant(s)" if "RETENU" in df.columns else "variant(s) annotés"
    print(f"\n→ {output_path}")
    print(f"  {len(df)} {label} — {df['SAMPLE'].nunique()} sample(s)")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="KODEX — Agrège les variants annotés (USER_CLASS/ANNOT/COM) de tous les *_filtered.xlsx."
    )
    parser.add_argument(
        "targets", type=Path, nargs="*",
        help="Répertoire OU fichier(s) xlsx (ex: run_1*_filtered.xlsx). "
             "Défaut: répertoire courant",
    )
    parser.add_argument(
        "--out", type=Path,
        help="Fichier de sortie (défaut: <répertoire>/<RUN>_annotations_<AAAAMMJJ_HHMM>.xlsx ; "
             "un chemin donné ici est utilisé tel quel et écrase s'il existe)",
    )
    parser.add_argument("--config", type=Path, default=None,
                        help="YAML config file (overrides --legacy)")
    parser.add_argument("--legacy", action="store_true",
                        help="Use the run 1 schema profile (config_legacy.yaml)")
    parser.add_argument(
        "--all", dest="every_variant", action="store_true",
        help="Agrège TOUS les variants (pas seulement les annotés), classés par "
             "score sur l'ensemble du run. Colonne RETENU = variant annoté. "
             "RHEB/PKD1 exclus. Sortie : <RUN>_all_variants_<date>.xlsx",
    )
    args = parser.parse_args(argv)

    here = Path(__file__).resolve().parent
    if args.config is not None:
        config_path = args.config
    elif args.legacy:
        config_path = here / "config_legacy.yaml"
    else:
        config_path = here / "config.yaml"
    cfg = load_config(config_path)
    igv_cfg = cfg.get("igv")
    scoring = cfg.get("scoring")

    targets = [t.expanduser().resolve() for t in args.targets] if args.targets else []

    # Déterminer la liste de fichiers et le répertoire de sortie
    if not targets:
        xlsx_dir = Path.cwd()
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

    # Horodaté par défaut : _write_summary ouvre le classeur en mode écriture, donc
    # un nom fixe recréait <RUN>_annotations.xlsx de zéro à chaque passage et
    # emportait toute édition faite dans ce fichier. --out reste littéral : le nom
    # donné explicitement est respecté tel quel, écrasement compris.
    stem = "all_variants" if args.every_variant else "annotations"
    out_path = (
        args.out
        or xlsx_dir / f"{_run_prefix(xlsx_dir)}_{stem}_{datetime.now():%Y%m%d_%H%M}.xlsx"
    ).expanduser().resolve()

    extra_genes = cfg.get("sheets", {}).get("extra_genes", [])
    # Profil-dépendant : GNOMAD4_AF en 069, MAX_AF_GNOMADEX2.1 en legacy. Les
    # deux colonnes coexistent dans les feuilles, donc se tromper est silencieux
    # — l'agrégat noterait la rareté autrement que le fichier patient d'origine.
    af_column = cfg.get("af_filter_column", "MAX_AF_GNOMADEX2.1")

    if files is not None:
        print(f"Lecture de {len(files)} fichier(s) ...")
    else:
        print(f"Lecture de {xlsx_dir} ...")
    df = aggregate_all(
        xlsx_dir, files=files, scoring=scoring,
        every_variant=args.every_variant, extra_genes=extra_genes,
        af_column=af_column,
    )

    if df.empty:
        print("Aucun variant trouvé." if args.every_variant
              else "Aucun variant annoté trouvé.")
        raise SystemExit(0)

    if args.every_variant and "RETENU" in df.columns:
        n_ret = int((df["RETENU"] == "oui").sum())
        if n_ret:
            ranks = df.index[df["RETENU"] == "oui"] + 1
            print(f"\n{len(df)} variants classés — {n_ret} retenu(s), "
                  f"rangs {ranks.min()} à {ranks.max()} "
                  f"(médiane {int(ranks.to_series().median())})")

    _write_summary(df, out_path, igv_cfg=igv_cfg, scoring=scoring,
                   af_column=af_column)


if __name__ == "__main__":
    main()
