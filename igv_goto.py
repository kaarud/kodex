"""igv_goto.py — KODEX : charge les BAMs d'un patient (+ référence) et saute sur une position dans IGV."""

import argparse
from pathlib import Path

import pandas as pd
import yaml

from igv_capture import check_igv, igv_cmd, extract_sample_id
from igv_link import build_locus, collect_bam_tracks, resolve_bam_dir


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
    # Relative bam_dir is anchored on the xlsx dir (mode xlsx), else on cwd.
    target = Path(args.target)
    anchor = target.parent if target.suffix.lower() == ".xlsx" else Path.cwd()
    bam_dir = resolve_bam_dir(args.bam_dir or cfg.get("bam_dir", "../bam"), anchor)
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
