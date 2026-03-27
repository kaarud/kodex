"""igv_capture.py — Generate IGV screenshots for annotated TSC1/TSC2 variants."""

import argparse
import socket
import time
from pathlib import Path

import pandas as pd
import yaml

import requests


_ANNOTATION_COLS = ["USER_CLASS", "USER_ANNOT", "USER_COM"]
_DEDUP_KEY = ["CHR", "POS", "REF", "ALT"]


def _matches_exclusion(variant: dict, terms: list[str]) -> bool:
    """Return True if any annotation column contains an exclusion term (case-insensitive)."""
    for col in _ANNOTATION_COLS:
        val = str(variant.get(col, "") or "").strip()
        if not val or val in ("nan", "None") or val.startswith("="):
            continue
        val_lower = val.lower()
        for term in terms:
            if term.lower() in val_lower:
                return True
    return False


def load_annotated_variants(
    xlsx_path: Path,
    exclusion_terms: list[str] | None = None,
) -> list[dict]:
    """Read all data sheets (skip _summary), return annotated variants deduplicated.

    If exclusion_terms is provided, rows whose USER_CLASS / USER_ANNOT / USER_COM
    contain any of the terms are excluded (case-insensitive substring match).
    """
    xl = pd.ExcelFile(xlsx_path, engine="openpyxl")
    sheets = [s for s in xl.sheet_names if s != "_summary"]

    frames = []
    for sheet in sheets:
        df = xl.parse(sheet)
        # Keep only rows with at least one annotation field non-empty
        # Ignore formula cells (start with "=") — these are cross-sheet references
        mask = pd.Series(False, index=df.index)
        for col in _ANNOTATION_COLS:
            if col in df.columns:
                vals = df[col].astype(str).str.strip()
                mask |= df[col].notna() & (vals != "") & (~vals.str.startswith("="))
        frames.append(df[mask])

    if not frames:
        return []

    combined = pd.concat(frames, ignore_index=True)

    # Deduplicate by genomic position
    available_key = [c for c in _DEDUP_KEY if c in combined.columns]
    if available_key:
        combined = combined.drop_duplicates(subset=available_key)

    variants = combined.to_dict(orient="records")

    # Apply exclusion terms
    if exclusion_terms:
        before = len(variants)
        variants = [v for v in variants if not _matches_exclusion(v, exclusion_terms)]
        excluded = before - len(variants)
        if excluded:
            print(f"  [{excluded} variant(s) exclu(s) par termes IGV : {', '.join(exclusion_terms)}]")

    return variants


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


def check_igv(port: int) -> bool:
    """Return True if IGV is reachable on the given port."""
    try:
        r = requests.get(f"http://localhost:{port}/", timeout=3)
        return r.status_code < 500
    except Exception:
        return False


def _build_cmd_str(command: str, **params) -> str:
    """Convert (command, **kwargs) to an IGV raw socket command string.

    IGV batch command format:
      new                              → "new"
      genome hg38                      → "genome hg38"
      goto chr9:100-200                → "goto chr9:100-200"
      setPreference SAM.COLOR_BY STRAND → "setPreference SAM.COLOR_BY STRAND"
      maxPanelHeight 350               → "maxPanelHeight 350"
      load /f.bam index=/f.bai name=X → "load /f.bam index=/f.bai name=X"
    """
    if not params:
        return command
    if command == "load":
        file_path = params.get("file", "")
        rest = " ".join(f"{k}={v}" for k, v in params.items() if k != "file")
        return f"load {file_path} {rest}".strip()
    if command == "setPreference":
        return f"setPreference {params['name']} {params['value']}"
    # Generic single-value commands: goto, genome, maxPanelHeight …
    val = next(iter(params.values()))
    return f"{command} {val}"


def igv_cmd(port: int, command: str, **params) -> str:
    """Send a single command to IGV via raw TCP socket.

    IGV 2.19.x (API v3.0): only a small subset of commands work as HTTP GET
    (/goto, /load). All others (new, genome, setPreference, maxPanelHeight …)
    are batch-only and must be sent as plain text over the raw socket.
    """
    cmd_str = _build_cmd_str(command, **params)
    try:
        with socket.create_connection(("localhost", port), timeout=10) as sock:
            sock.sendall((cmd_str + "\n").encode())
            sock.settimeout(5)
            resp = b""
            try:
                while True:
                    chunk = sock.recv(1024)
                    if not chunk:
                        break
                    resp += chunk
            except socket.timeout:
                pass
            text = resp.decode("utf-8", errors="replace").strip()
            if text.upper().startswith("ERROR"):
                print(f"[WARN] IGV: {command} → {text}")
            return text
    except OSError as e:
        print(f"[ERREUR] IGV command '{command}' failed: {e}")
        raise


def igv_batch(port: int, lines: list[str]) -> None:
    """Send batch commands to IGV via raw TCP socket (one connection per command).

    IGV 2.19.x (HTTP API v3.0): snapshot/snapshotDirectory/setPreference are not
    reachable via HTTP GET. IGV's "HTTP server" is actually a raw TCP socket server
    that reads the first line as a command. POST requests cause a BadStatusLine error.
    Sending raw text (no HTTP headers) works for all batch commands.
    """
    for cmd in lines:
        try:
            with socket.create_connection(("localhost", port), timeout=10) as sock:
                sock.sendall((cmd + "\n").encode())
                sock.settimeout(5)
                resp = b""
                try:
                    while True:
                        chunk = sock.recv(1024)
                        if not chunk:
                            break
                        resp += chunk
                except socket.timeout:
                    pass
                text = resp.decode("utf-8", errors="replace").strip()
                if text.upper().startswith("ERROR"):
                    print(f"[WARN] IGV: {cmd} → {text}")
        except OSError as e:
            print(f"[ERREUR] IGV socket '{cmd}' failed: {e}")
            raise


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
    """Generate 3 screenshots per variant (strand, softclip, squished).

    IGV 2.19.x (API v3.0): navigation/preference commands via GET,
    snapshot commands via POST /batch.
    """
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    chrom = normalize_chr(str(variant.get("CHR", "")))
    pos   = int(variant["POS"])
    dp    = int(variant.get("DP") or 100)
    name  = build_variant_name(variant)

    panel_height = dp * 2 + 150
    locus = f"{chrom}:{pos - _HALF_WINDOW}-{pos + _HALF_WINDOW}"

    strand_fname   = f"{name}_strand.png"
    softclip_fname = f"{name}_softclip.png"
    squished_fname = f"{name}_squished.png"

    # Navigate and set display via GET endpoints
    igv_cmd(port, "setPreference", name="SAM.DISPLAY_MODE",   value="SQUISHED")
    igv_cmd(port, "setPreference", name="SAM.SAMPLING_COUNT", value=str(dp))
    igv_cmd(port, "maxPanelHeight", value=str(panel_height))
    igv_cmd(port, "goto", locus=locus)
    time.sleep(delay)

    # All 3 captures in a single batch script (POST /batch)
    # strand ALWAYS before soft clips — reset order matters
    igv_batch(port, [
        f"snapshotDirectory {out_dir}",
        "setPreference SAM.COLOR_BY READ_STRAND",
        f"snapshot {strand_fname}",
        "setPreference SAM.SHOW_SOFT_CLIPPED true",
        "setPreference SAM.COLOR_BY NO_COLORING",
        f"snapshot {softclip_fname}",
        "setPreference SAM.SHOW_SOFT_CLIPPED false",
        f"snapshot {squished_fname}",
    ])

    return {
        "strand":   out_dir / strand_fname,
        "softclip": out_dir / softclip_fname,
        "squished": out_dir / squished_fname,
    }


def find_bams(sample_id: str, bam_dir: Path) -> dict | None:
    """Locate the 3 BAMs + indices for a sample. Returns None if any missing."""
    specs = {
        "raw":     bam_dir / f"{sample_id}.bam",
        "mutect2": bam_dir / f"{sample_id}.mutect2.bam",
        "chim":    bam_dir / f"{sample_id}.chim.bam",
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
    parser.add_argument("--config",  type=Path,
                        default=Path(__file__).parent / "config.yaml")
    args = parser.parse_args(argv)

    xlsx_path = args.xlsx.expanduser().resolve()
    bam_dir   = args.bam_dir.expanduser().resolve()
    out_dir   = args.out.expanduser().resolve()

    # Load exclusion terms from config
    exclusion_terms = []
    if args.config.exists():
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}
        exclusion_terms = cfg.get("igv_exclusion_terms", [])

    # 1. Charger les variants annotés
    variants = load_annotated_variants(xlsx_path, exclusion_terms=exclusion_terms)
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
        print(f"    → {result['squished'].name}")

    print(f"\nCaptures terminées → {sample_out}")


if __name__ == "__main__":
    main()
