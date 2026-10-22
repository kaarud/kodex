"""minibam.py — KODEX : tranches BAM légères autour des variants annotés.

Lit la feuille ANNOTATIONS produite par aggregate_annotations.py et découpe,
pour chaque échantillon, les BAM du run (RAW + CHIM + MUTECT2 + comparateur)
sur les seules régions des variants retenus (± flanc). Résultat : quelques Mo
transférables au lieu de plusieurs Go, ouvrables directement dans IGV.
"""

import argparse
import shutil
import subprocess
import sys
import tarfile
from datetime import datetime
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from igv_capture import normalize_chr
from igv_link import collect_bam_tracks, resolve_bam_dir
from vcf_filter import load_config

SHEET = "ANNOTATIONS"
REQUIRED_COLS = ("SAMPLE", "CHR", "POS")


# ---------------------------------------------------------------------------
# Logique pure — testable sans samtools ni BAM
# ---------------------------------------------------------------------------

def merge_regions(positions: list[tuple[str, int]], flank: int = 500) -> list[str]:
    """Retourne les régions samtools ± flanc, dédupliquées et fusionnées.

    Deux variants distants de moins de 2 × flanc partagent une seule tranche —
    inutile d'extraire deux fois les mêmes reads.
    """
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, pos in positions:
        chrom = normalize_chr(str(chrom))
        start = max(1, int(pos) - flank)
        by_chrom.setdefault(chrom, []).append((start, int(pos) + flank))

    regions = []
    for chrom in sorted(by_chrom):
        merged: list[list[int]] = []
        for start, end in sorted(by_chrom[chrom]):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        regions.extend(f"{chrom}:{s}-{e}" for s, e in merged)
    return regions


def slice_name(bam_path: Path, label: str = "RAW") -> str:
    """'X_S7.chim.bam' → 'X_S7.chim.mini.bam' (le type de piste est conservé).

    La tranche de comparaison est préfixée REF_ : elle porte le nom d'un autre
    patient et se retrouve dans le dossier du patient courant.
    """
    name = Path(bam_path).name[: -len(".bam")] + ".mini.bam"
    return f"REF_{name}" if label.startswith("REF") else name


def build_view_cmd(bam: Path, regions: list[str], out: Path) -> list[str]:
    """Commande samtools extrayant toutes les régions en une passe."""
    return ["samtools", "view", "-b", "-h", "-o", str(out), str(bam), *regions]


def write_manifest(entries: list[dict], path: Path) -> None:
    """Tableau récapitulatif — quel fichier vient de quel BAM, quelle taille."""
    cols = ["sample", "track", "file", "source_bam", "n_variants", "n_regions", "bytes"]
    pd.DataFrame(entries, columns=cols).to_csv(path, sep="\t", index=False)


def read_annotations(xlsx_path: Path) -> pd.DataFrame:
    """Charge la feuille ANNOTATIONS et vérifie les colonnes indispensables."""
    sheets = pd.ExcelFile(xlsx_path, engine="openpyxl").sheet_names
    if SHEET not in sheets:
        raise SystemExit(
            f"[ERREUR] Feuille '{SHEET}' absente de {xlsx_path.name}. "
            "Attendu : le fichier produit par aggregate_annotations.py."
        )
    df = pd.read_excel(xlsx_path, sheet_name=SHEET, engine="openpyxl")
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise SystemExit(f"[ERREUR] Colonnes manquantes dans {SHEET} : {', '.join(missing)}")
    df = df.dropna(subset=list(REQUIRED_COLS))
    if df.empty:
        raise SystemExit(f"[ERREUR] Aucun variant exploitable dans {SHEET}.")
    return df


# ---------------------------------------------------------------------------
# Extraction
# ---------------------------------------------------------------------------

def run_extraction(
    df: pd.DataFrame,
    bam_dir: Path,
    out_dir: Path,
    reference_lane: str = "S2",
    flank: int = 500,
    dry_run: bool = False,
) -> list[dict]:
    """Découpe les BAM de chaque sample. Retourne les entrées du manifest."""
    entries: list[dict] = []

    for sample, rows in df.groupby("SAMPLE", sort=True):
        sample = str(sample)
        tracks = collect_bam_tracks(sample, bam_dir, reference_lane)
        if tracks is None:
            print(f"  {sample:<18s} [SKIP] aucun BAM dans {bam_dir}")
            continue

        regions = merge_regions(list(zip(rows["CHR"], rows["POS"])), flank=flank)
        sample_dir = out_dir / sample
        if not dry_run:
            sample_dir.mkdir(parents=True, exist_ok=True)

        print(f"  {sample:<18s} {len(rows)} variant(s) → {len(regions)} région(s), "
              f"{len(tracks)} piste(s)")

        for bam, label in tracks:
            out = sample_dir / slice_name(bam, label)
            cmd = build_view_cmd(bam, regions, out)
            if dry_run:
                print("    " + " ".join(cmd))
                continue
            subprocess.run(cmd, check=True)
            subprocess.run(["samtools", "index", str(out)], check=True)
            entries.append({
                "sample": sample,
                "track": label,
                "file": str(out.relative_to(out_dir)),
                "source_bam": str(bam),
                "n_variants": len(rows),
                "n_regions": len(regions),
                "bytes": out.stat().st_size,
            })

    return entries


_README = """\
Mini-BAM KODEX — run {run}
Généré le {date}

Contenu : pour chaque échantillon, les BAM du pipeline NiourK tronqués aux
seules régions des variants annotés, avec {flank} pb de part et d'autre.

  <SAMPLE>/<SAMPLE>.mini.bam           lecture brute (RAW)
  <SAMPLE>/<SAMPLE>.chim.mini.bam      reads chimériques
  <SAMPLE>/<SAMPLE>.mutect2.mini.bam   réalignement Mutect2 (si présent)
  <SAMPLE>/REF_*.mini.bam              échantillon de comparaison du run

Ouverture dans IGV : génome hg38, puis glisser-déposer les .bam voulus
(les .bai doivent rester à côté). En dehors des régions listées dans
manifest.tsv les pistes sont vides — c'est normal, les reads n'ont pas été
extraits.
"""


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="KODEX — Découpe des mini-BAM autour des variants de la feuille d'agrégation."
    )
    parser.add_argument("xlsx", type=Path, help="Fichier <RUN>_annotations.xlsx")
    parser.add_argument("--bam-dir", type=Path, default=None,
                        help="Répertoire des BAM (défaut : igv.bam_dir du config)")
    parser.add_argument("--out", type=Path, default=None,
                        help="Répertoire de sortie (défaut : <dir xlsx>/minibam/<RUN>_<horodatage>)")
    parser.add_argument("--flank", type=int, default=500,
                        help="Paires de bases de part et d'autre du variant (défaut : 500)")
    parser.add_argument("--config", type=Path, default=None,
                        help="YAML config file (overrides --legacy)")
    parser.add_argument("--legacy", action="store_true",
                        help="Use the run 1 schema profile (config_legacy.yaml)")
    parser.add_argument("--archive", action="store_true",
                        help="Produire aussi une archive .tar.gz prête à envoyer")
    parser.add_argument("--dry-run", action="store_true",
                        help="Afficher les commandes samtools sans rien écrire")
    args = parser.parse_args(argv)

    if shutil.which("samtools") is None:
        raise SystemExit("[ERREUR] samtools introuvable dans le PATH")

    here = Path(__file__).resolve().parent
    if args.config is not None:
        config_path = args.config
    elif args.legacy:
        config_path = here / "config_legacy.yaml"
    else:
        config_path = here / "config.yaml"
    igv_cfg = load_config(config_path).get("igv") or {}

    xlsx_path = args.xlsx.expanduser().resolve()
    df = read_annotations(xlsx_path)

    bam_dir = (args.bam_dir.expanduser().resolve() if args.bam_dir
               else resolve_bam_dir(igv_cfg.get("bam_dir", "../bam"), xlsx_path.parent))
    if not bam_dir.is_dir():
        raise SystemExit(f"[ERREUR] Répertoire BAM introuvable : {bam_dir}")

    run = xlsx_path.stem.replace("_annotations", "")
    stamp = datetime.now().strftime("%Y%m%d-%H%M")
    # Sortie horodatée : jamais d'écrasement d'un envoi déjà préparé.
    out_dir = (args.out.expanduser().resolve() if args.out
               else xlsx_path.parent / "minibam" / f"{run}_{stamp}")

    print(f"BAM   : {bam_dir}")
    print(f"Sortie: {out_dir}{'  (dry-run)' if args.dry_run else ''}")
    print(f"Flanc : ±{args.flank} pb\n")

    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)

    entries = run_extraction(
        df, bam_dir, out_dir,
        reference_lane=igv_cfg.get("reference_lane", "S2"),
        flank=args.flank,
        dry_run=args.dry_run,
    )

    if args.dry_run:
        return
    if not entries:
        print("\nAucun mini-BAM produit.")
        raise SystemExit(1)

    write_manifest(entries, out_dir / "manifest.tsv")
    (out_dir / "README.txt").write_text(
        _README.format(run=run, date=datetime.now().strftime("%d/%m/%Y %H:%M"), flank=args.flank),
        encoding="utf-8",
    )

    total = sum(e["bytes"] for e in entries)
    print(f"\n→ {out_dir}")
    print(f"  {len(entries)} fichier(s) — {total / 1e6:.1f} Mo")

    if args.archive:
        archive = out_dir.with_suffix(".tar.gz")
        with tarfile.open(archive, "w:gz") as tar:
            tar.add(out_dir, arcname=out_dir.name)
        print(f"  archive : {archive} ({archive.stat().st_size / 1e6:.1f} Mo)")


if __name__ == "__main__":
    main()
