"""igv_link.py — KODEX : logique pure pour piloter IGV (URL + locus + BAM set)."""

from pathlib import Path
from urllib.parse import quote

from igv_capture import normalize_chr
from openpyxl.styles import Alignment, Font
from openpyxl.worksheet.worksheet import Worksheet

# Clickable ▶ glyphs: oversized on purpose — they are the main hit target of the
# sheet, clicked hundreds of times during a review.
_GLYPH_FONT = Font(color="0563C1", underline="single", size=16, bold=True)
_GLYPH_ALIGN = Alignment(horizontal="center", vertical="center")


def build_locus(chrom: str, pos: int, half_window: int = 145) -> str:
    """Return an IGV locus string centered on pos, e.g. 'chr9:855-1145'."""
    chrom = normalize_chr(chrom)
    pos = int(pos)
    start = max(1, pos - half_window)
    end = pos + half_window
    return f"{chrom}:{start}-{end}"


def resolve_bam_dir(bam_dir: str | Path, anchor: Path) -> Path:
    """Resolve the configured bam_dir, expanding '~'.

    A *relative* bam_dir is resolved against `anchor` — the directory holding
    the XLSX being written. This keeps a run self-contained: '../BAM' follows
    the run wherever the directory is moved, whereas an absolute path silently
    breaks as soon as the run is relocated.
    """
    p = Path(bam_dir).expanduser()
    return p if p.is_absolute() else (Path(anchor) / p).resolve()


def resolve_reference_bam(
    bam_dir: Path,
    reference_lane: str,
    current_sample_id: str,
) -> Path | None:
    """Return the raw .bam used as the comparison track.

    Normally the sample whose id ends with _{lane}. When that sample *is* the
    current one — comparing a patient to itself tells you nothing — the first
    other raw BAM of the directory (alphabetical) is used instead. Same fallback
    when the lane is absent from the run.

    Returns None only when no other raw BAM exists at all. .chim.bam and
    .mutect2.bam are never considered.
    """
    bam_dir = Path(bam_dir)
    suffix = f"_{reference_lane}"
    others = []                                   # raw BAMs != current sample
    for bam in sorted(bam_dir.glob("*.bam")):
        stem = bam.name[:-len(".bam")]            # strip trailing .bam
        if stem.endswith(".chim") or stem.endswith(".mutect2"):
            continue
        if stem == current_sample_id:
            continue
        if stem.endswith(suffix):
            return bam
        others.append(bam)
    return others[0] if others else None


def reference_label(ref_bam: Path, reference_lane: str) -> str:
    """Track label for the comparison BAM.

    'REF_S2' when it really is the configured lane, 'REF_<sample_id>' when a
    fallback sample was picked — the reviewer must see *who* they compare to.
    """
    stem = Path(ref_bam).name[:-len(".bam")]
    if stem.endswith(f"_{reference_lane}"):
        return f"REF_{reference_lane}"
    return f"REF_{stem}"


def collect_bam_tracks(
    sample_id: str,
    bam_dir: Path,
    reference_lane: str,
) -> list[tuple[Path, str]] | None:
    """Return ordered (bam_path, label) tuples to load into IGV.

    Order: patient RAW, then CHIM and MUTECT2 if those files exist, then the
    reference raw BAM (REF_{lane}). Only the patient RAW BAM is mandatory —
    Run 2 ships raw+chim (no mutect2), run 1 ships all three. Returns
    None only if the patient RAW BAM is absent. The reference track is omitted
    if missing or if it is the current sample.
    """
    bam_dir = Path(bam_dir)
    raw = bam_dir / f"{sample_id}.bam"
    if not raw.exists():
        return None
    tracks: list[tuple[Path, str]] = [(raw, "RAW")]
    for suffix, label in ((".chim.bam", "CHIM"), (".mutect2.bam", "MUTECT2")):
        bam = bam_dir / f"{sample_id}{suffix}"
        if bam.exists():
            tracks.append((bam, label))
    ref = resolve_reference_bam(bam_dir, reference_lane, sample_id)
    if ref is not None:
        tracks.append((ref, reference_label(ref, reference_lane)))
    return tracks


def build_load_url(
    tracks: list[tuple[Path, str]],
    locus: str | None = None,
    port: int = 60151,
    genome: str = "hg38",
) -> str:
    """Build a single IGV HTTP GET that loads all tracks (one-time setup).

    Used by the _summary "Load BAMs" link. merge=false clears previously loaded
    tracks. File paths are URL-encoded (slashes kept); commas separate
    files/names. locus is optional — omit it to load without moving the view.
    """
    files = ",".join(quote(str(p.resolve()), safe="/") for p, _ in tracks)
    names = ",".join(label for _, label in tracks)
    locus_part = f"&locus={quote(locus, safe='-')}" if locus else ""
    return (
        f"http://localhost:{port}/load?"
        f"file={files}"
        f"&genome={genome}&merge=false"
        f"{locus_part}"
        f"&name={names}"
    )


def build_goto_url(locus: str, port: int = 60151) -> str:
    """Build a lightweight IGV navigation GET — no BAM (re)load, just goto.

    Used per-row in the data sheets: the _summary Load link loads the BAMs once,
    then each ▶ only moves the view. Avoids reloading tracks (and the F5/merge
    glitch) on every click.
    """
    locus_q = quote(locus, safe="-")          # keep '-', encode ':'
    return f"http://localhost:{port}/goto?locus={locus_q}"


def build_sort_base_url(port: int = 60151) -> str:
    """Build the IGV 'sort by base' GET — sorts at the current view, no goto.

    IGV's HTTP layer only exposes /load and /goto; 'sort' rides the undocumented
    /execute pass-through, which runs exactly one batch command per request (so
    goto+sort can't be chained in a single click). This is the per-sheet A1
    button: click ▶ to jump, then A1 to re-sort by base at that position.
    """
    return f"http://localhost:{port}/execute?command={quote('sort base')}"


def inject_summary_load_link(
    ws: "Worksheet",
    sample_id: str,
    bam_dir: Path,
    reference_lane: str,
    port: int = 60151,
    genome: str = "hg38",
) -> bool:
    """Append a one-click 'Load BAMs' hyperlink to the _summary sheet.

    Click it once per patient: it loads RAW (+CHIM/+MUTECT2 if present) plus the
    reference REF_{lane} into IGV. The per-row ▶ links then only navigate.
    Returns False (no-op) if the patient BAMs are missing.
    """
    tracks = collect_bam_tracks(sample_id, bam_dir, reference_lane)
    if tracks is None:
        return False

    url = build_load_url(tracks, locus=None, port=port, genome=genome)
    labels = ", ".join(label for _, label in tracks)
    r = ws.max_row + 2                          # blank line, then the link
    ws.cell(row=r, column=1, value="IGV")
    cell = ws.cell(row=r, column=2, value=f"▶ Charger les BAMs ({labels})")
    cell.hyperlink = url
    cell.font = Font(color="0563C1", underline="single")
    return True


def inject_igv_column(
    ws: "Worksheet",
    sample_id: str,
    bam_dir: Path,
    reference_lane: str,
    port: int = 60151,
    half_window: int = 145,
) -> bool:
    """Insert a leftmost IGV column. The A1 header is a '⇅ base' button that
    sorts the current view by base; each row below is a ▶ goto hyperlink.
    Returns False (no-op) if the patient's BAMs are missing — the sheet is left
    untouched.

    Each ▶ only navigates IGV (goto); the BAMs are loaded once via the _summary
    Load link. Sorting can't ride the same click (IGV runs one HTTP command per
    request), so A1 is a separate one-click 'sort base' button: jump with ▶, then
    sort with A1. The BAM presence check gates the column.
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
    hdr = ws.cell(row=1, column=1, value="⇅ base")   # A1 = click to sort by base
    hdr.hyperlink = build_sort_base_url(port=port)
    hdr.font = Font(color="0563C1", underline="single")
    for r in range(2, ws.max_row + 1):
        chrom = ws.cell(row=r, column=chr_col + 1).value   # +1: column shifted right
        pos = ws.cell(row=r, column=pos_col + 1).value
        if chrom in (None, "") or pos in (None, ""):
            continue
        locus = build_locus(str(chrom), int(pos), half_window=half_window)
        url = build_goto_url(locus, port=port)
        cell = ws.cell(row=r, column=1, value="▶")
        cell.hyperlink = url
        cell.font = _GLYPH_FONT
        cell.alignment = _GLYPH_ALIGN
    ws.column_dimensions["A"].width = 6
    return True


def inject_aggregate_igv_columns(
    ws: "Worksheet",
    bam_dir: Path,
    reference_lane: str,
    port: int = 60151,
    genome: str = "hg38",
    half_window: int = 145,
) -> int:
    """Insert two leftmost IGV columns into the aggregated multi-sample sheet.

    Unlike a per-sample data sheet, every row here may belong to a different
    patient, so loading BAMs can't be a one-time _summary link — it has to be
    per row:

      - A1 is the '⇅ base' sort-by-base button (one command per click).
      - Column A ('▶ BAM'), per row: loads that row's sample BAMs (RAW +
        CHIM/MUTECT2 if present + REF_{lane}) AND jumps to the variant locus.
        Click it when switching patients.
      - Column B ('▶'), per row: lightweight goto — re-navigate within the
        already-loaded patient without reloading tracks.

    Needs SAMPLE, CHR and POS header columns. Rows whose sample has no BAMs are
    left blank (tracks resolved once per sample, cached). Returns the number of
    rows that received links (0 if no BAMs found at all, or headers missing).
    """
    bam_dir = Path(bam_dir)
    headers = {ws.cell(row=1, column=c).value: c for c in range(1, ws.max_column + 1)}
    sample_col = headers.get("SAMPLE")
    chr_col, pos_col = headers.get("CHR"), headers.get("POS")
    if sample_col is None or chr_col is None or pos_col is None:
        return 0

    ws.insert_cols(1, 2)
    sample_col += 2          # original columns shifted right by 2
    chr_col += 2
    pos_col += 2

    a1 = ws.cell(row=1, column=1, value="⇅ base")   # A1 = sort by base
    a1.hyperlink = build_sort_base_url(port=port)
    a1.font = Font(color="0563C1", underline="single")
    ws.cell(row=1, column=2, value="goto")

    track_cache: dict[str, list | None] = {}
    n = 0
    for r in range(2, ws.max_row + 1):
        sample = ws.cell(row=r, column=sample_col).value
        chrom = ws.cell(row=r, column=chr_col).value
        pos = ws.cell(row=r, column=pos_col).value
        if sample in (None, "") or chrom in (None, "") or pos in (None, ""):
            continue
        sample = str(sample)
        if sample not in track_cache:
            track_cache[sample] = collect_bam_tracks(sample, bam_dir, reference_lane)
        tracks = track_cache[sample]
        if tracks is None:
            continue
        locus = build_locus(str(chrom), int(pos), half_window=half_window)
        load = ws.cell(row=r, column=1, value="▶ BAM")
        load.hyperlink = build_load_url(tracks, locus=locus, port=port, genome=genome)
        load.font = _GLYPH_FONT
        load.alignment = _GLYPH_ALIGN
        goto = ws.cell(row=r, column=2, value="▶")
        goto.hyperlink = build_goto_url(locus, port=port)
        goto.font = _GLYPH_FONT
        goto.alignment = _GLYPH_ALIGN
        n += 1
    ws.column_dimensions["A"].width = 10
    ws.column_dimensions["B"].width = 6
    return n
