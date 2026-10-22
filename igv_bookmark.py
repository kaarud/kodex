"""igv_bookmark.py — KODEX : injecte les variants annotés comme ROI dans IGV."""

import argparse
from pathlib import Path

import yaml

from igv_capture import (
    load_annotated_variants,
    check_igv,
    normalize_chr,
    build_variant_name,
    igv_batch,
)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="KODEX — Injecte les variants annotés d'un XLSX filtré comme ROI dans IGV."
    )
    parser.add_argument("xlsx", type=Path, help="Fichier *_filtered.xlsx")
    parser.add_argument("--port", type=int, default=60151)
    parser.add_argument("--config", type=Path,
                        default=Path(__file__).parent / "config.yaml")
    args = parser.parse_args(argv)

    xlsx_path = args.xlsx.expanduser().resolve()

    # Load exclusion terms from config
    exclusion_terms = []
    if args.config.exists():
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}
        exclusion_terms = cfg.get("igv_exclusion_terms", [])

    variants = load_annotated_variants(xlsx_path, exclusion_terms=exclusion_terms)
    if not variants:
        print("Aucun variant annoté trouvé (USER_CLASS / USER_ANNOT / USER_COM vides).")
        raise SystemExit(0)
    print(f"{len(variants)} variant(s) annoté(s) à marquer.")

    if not check_igv(args.port):
        print(f"[ERREUR] IGV non joignable sur le port {args.port}. Ouvrir IGV Desktop et réessayer.")
        raise SystemExit(1)

    # Build ROI entries
    roi_entries = []
    for v in variants:
        chrom = normalize_chr(str(v.get("CHR", "")))
        pos = v.get("POS")
        if not chrom or not pos:
            continue
        pos = int(pos)
        name = build_variant_name(v)

        # Build description from user annotations only (skip formulas)
        desc_parts = []
        for col in ("USER_CLASS", "USER_ANNOT", "USER_COM"):
            val = str(v.get(col, "") or "").strip()
            if val and val not in ("", "nan", "None") and not val.startswith("="):
                desc_parts.append(val)
        desc = " | ".join(desc_parts)
        full_desc = f"{name}  {desc}".rstrip() if desc else name

        roi_entries.append((chrom, pos, full_desc))

    # Write a BED file for IGV Regions > Import Regions (tab-separated → spaces OK)
    bed_path = xlsx_path.with_suffix(".roi.bed")
    with open(bed_path, "w") as fh:
        for chrom, pos, desc in roi_entries:
            fh.write(f"{chrom}\t{pos}\t{pos}\t{desc}\n")

    # Also inject via batch region commands (description underscored for IGV parser)
    commands = []
    for chrom, pos, desc in roi_entries:
        safe_desc = desc.replace(" ", "_")
        commands.append(f"region {chrom} {pos} {pos} {safe_desc}")

    igv_batch(args.port, commands)
    print(f"{len(commands)} région(s) injectée(s).")
    print(f"BED exporté → {bed_path}")
    print("→ IGV : Regions > Import Regions pour charger le BED (descriptions complètes)")
    print("→ IGV : Regions > Region Navigator pour naviguer entre les variants.")


if __name__ == "__main__":
    main()
