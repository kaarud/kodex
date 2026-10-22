"""bg_filter.py — KODEX : retrouver les variants qui ressemblent au temoin positif.

`run_background.py` ecrit sa mesure dans un commentaire de cellule, lisible a
l'oeil mais pas triable : sur 726 lignes on ne repere pas a la main les quelques
sites ou le patient porte des reads alternatifs que personne d'autre n'a. Ce
script relit ces commentaires, en extrait les chiffres, et ne garde que la
signature du variant qui a motive la mesure — TSC2 c.4090dup chez 100015_S14 :
AF 2,6 %, bruit du run nul, contraste x47, reads primaires intacts.

Ce qui discrimine n'est pas la VAF basse (les artefacts aussi) mais la
**conjonction** : bruit median nul ET tres peu d'echantillons avec le moindre
read alternatif ET contraste eleve ET rapport primaire proche de 1. Chaque
critere pris seul laisse passer une classe d'artefacts entiere.

Usage :
    python bg_filter.py fichier_bg.xlsx [autres...] [-o sortie.xlsx]
"""

import re
from pathlib import Path

# Les chiffres du commentaire, dans l'ordre ou format_comment les ecrit.
_RE_PATIENT = re.compile(r"->\s*alt\s+([\d.]+|—)\s*%")
_RE_DP = re.compile(r"DP=(\d+)")
_RE_STATS = re.compile(
    r"mediane\s+([\d.]+|—)\s*%\s*\|\s*max\s+([\d.]+|—)\s*%"
    r"\s*\|\s*avec alt\s*:\s*(\d+)/(\d+)")
_RE_CONTRAST = re.compile(r"contraste\s+x([\d.]+)")
_RE_RATIO = re.compile(r"rapport\s+([\d.]+)")
_RE_CARRIERS = re.compile(r"porteurs exclus du bruit\s*:\s*(.+)")


def parse_comment(text: str) -> dict:
    """Chiffres du commentaire DP. Valeur absente -> None, jamais 0.

    La distinction vaut la peine : « pas mesure » et « mesure a zero » menent a
    des conclusions opposees sur un site.
    """
    num = lambda m, g=1: (None if m is None or m.group(g) == "—"
                          else float(m.group(g)))
    stats = _RE_STATS.search(text)
    carriers = _RE_CARRIERS.search(text)
    return {
        "AF_BAM_PCT": num(_RE_PATIENT.search(text)),
        "DP_BAM": int(_RE_DP.search(text).group(1)) if _RE_DP.search(text) else 0,
        "BG_MEDIAN_PCT": num(stats, 1) if stats else None,
        "BG_MAX_PCT": num(stats, 2) if stats else None,
        "BG_N_WITH_ALT": int(stats.group(3)) if stats else None,
        "BG_N_MEASURED": int(stats.group(4)) if stats else None,
        "BG_CONTRAST": num(_RE_CONTRAST.search(text)),
        "BG_PRIMARY_RATIO": num(_RE_RATIO.search(text)),
        "BG_CARRIERS": carriers.group(1).strip() if carriers else "",
    }


def keeps(m: dict, max_af=30.0, min_af=0.5, min_contrast=5.0,
          max_median=0.0, max_with_alt_frac=0.35, min_ratio=0.8,
          min_dp=100) -> bool:
    """Le variant a-t-il la signature du temoin positif ?

    `max_median = 0` est volontairement strict : c'est le critere qui a fait
    remarquer c.4090dup — 0 read alternatif sur 22 echantillons. L'assouplir
    (0,1 %) ouvre aux sites legerement bruites, au prix du volume.

    `max_with_alt_frac` est cale sur le temoin positif du run 2 :
    c.2220+1_2220+2del (0,69 %, contraste x7,3) a 5 echantillons sur 18 portant
    un ou deux reads de trace, soit 0,28. Un seuil a 0,15 le rejetait. A basse
    VAF, quelques reads parasites chez les autres sont la norme, pas un signe
    de bruit partage — c'est la mediane et le contraste qui portent ce verdict.

    `min_ratio` accepte un rapport absent (indel court non mesurable) plutot que
    de rejeter : l'absence de mesure n'est pas une charge contre le variant.
    """
    if (m.get("DP_BAM") or 0) < min_dp:
        # Sous une centaine de reads, le plancher 1/profondeur gonfle le
        # contraste sans qu'aucune fraction ne soit reellement mesuree.
        return False
    af, contrast = m["AF_BAM_PCT"], m["BG_CONTRAST"]
    med, n_alt, n_meas = m["BG_MEDIAN_PCT"], m["BG_N_WITH_ALT"], m["BG_N_MEASURED"]
    if af is None or not (min_af <= af <= max_af):
        return False
    if contrast is None or contrast < min_contrast:
        return False
    if med is None or med > max_median:
        return False
    if not n_meas or n_alt / n_meas > max_with_alt_frac:
        return False
    ratio = m["BG_PRIMARY_RATIO"]
    return ratio is None or ratio >= min_ratio


_KEEP_COLS = ("SAMPLE", "SYMBOL", "HGVSc", "HGVSp", "AF_PCT", "DP", "SCORE",
              "IMPACT", "Consequence", "CLINVAR", "LOVD", "CALLNB",
              "SOURCE_SHEET", "CHR", "POS", "REF", "ALT",
              "MAX_AF_GNOMADEX2.1", "GNOMAD4_AF", "CADD_PHRED", "REVEL",
              "SpliceAI_num", "RNCNT_TOTAL", "PJCNT_TOTAL")


def scan(path: Path, impacts=None, **thresholds) -> list[dict]:
    """Toutes les lignes commentees d'un classeur qui passent le filtre."""
    from openpyxl import load_workbook

    wb = load_workbook(path)
    rows = []
    for ws in wb.worksheets:
        if ws.title.startswith("_"):
            continue
        idx = {c.value: c.column for c in ws[1]}
        if "DP" not in idx:
            continue
        for r in range(2, ws.max_row + 1):
            comment = ws.cell(r, idx["DP"]).comment
            if not comment:
                continue
            m = parse_comment(comment.text)
            if not keeps(m, **thresholds):
                continue
            if impacts and "IMPACT" in idx:
                if str(ws.cell(r, idx["IMPACT"]).value).upper() not in impacts:
                    continue
            row = {"FILE": path.name, "SHEET": ws.title, "ROW": r}
            row.update({c: ws.cell(r, idx[c]).value for c in _KEEP_COLS
                        if c in idx})
            row.update(m)
            rows.append(row)
    return rows


def main(argv=None):
    import argparse
    import pandas as pd

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("xlsx", type=Path, nargs="+", help="classeurs *_bg.xlsx")
    ap.add_argument("-o", "--out", type=Path, default=Path("bg_candidats.xlsx"))
    ap.add_argument("--min-af", type=float, default=0.5)
    ap.add_argument("--max-af", type=float, default=30.0)
    ap.add_argument("--min-contrast", type=float, default=5.0)
    ap.add_argument("--max-median", type=float, default=0.0,
                    help="bruit median tolere, en %% (defaut 0 = aucun)")
    ap.add_argument("--max-with-alt-frac", type=float, default=0.35)
    ap.add_argument("--min-ratio", type=float, default=0.8)
    ap.add_argument("--min-dp", type=int, default=100,
                    help="profondeur minimale du patient dans le BAM")
    ap.add_argument("--impact", nargs="*", default=None,
                    help="ne garder que ces IMPACT (ex: HIGH MODERATE LOW)")
    a = ap.parse_args(argv)

    th = dict(min_af=a.min_af, max_af=a.max_af, min_contrast=a.min_contrast,
              max_median=a.max_median, max_with_alt_frac=a.max_with_alt_frac,
              min_ratio=a.min_ratio, min_dp=a.min_dp)
    impacts = {i.upper() for i in a.impact} if a.impact else None
    rows = [r for p in a.xlsx for r in scan(p, impacts=impacts, **th)]
    if not rows:
        print("Aucun variant ne passe le filtre.")
        return
    df = pd.DataFrame(rows).sort_values("BG_CONTRAST", ascending=False)
    df.to_excel(a.out, index=False)
    print(f"{len(df)} variants -> {a.out}")
    cols = [c for c in ("SAMPLE", "SYMBOL", "HGVSc", "AF_BAM_PCT", "DP_BAM",
                        "BG_CONTRAST", "BG_N_WITH_ALT", "BG_N_MEASURED",
                        "BG_PRIMARY_RATIO", "IMPACT") if c in df]
    print(df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
