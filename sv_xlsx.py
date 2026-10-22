"""sv_xlsx.py — KODEX : SV VCF (Manta / Delly) → Excel avec liens IGV."""

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import openpyxl
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter

import yaml

from igv_capture import normalize_chr
from igv_link import (
    build_goto_url,
    build_locus,
    build_load_url,
    reference_label,
    resolve_bam_dir,
    resolve_reference_bam,
)


_SVTYPE_FILL = {
    "DEL": PatternFill("solid", fgColor="FFD0D0"),
    "DUP": PatternFill("solid", fgColor="D0E8FF"),
    "INV": PatternFill("solid", fgColor="FFFFD0"),
    "BND": PatternFill("solid", fgColor="F0D0FF"),
    "INS": PatternFill("solid", fgColor="D0FFD8"),
}
_LOWQUAL_FILL = PatternFill("solid", fgColor="E0E0E0")
_LINK_FONT = Font(color="0563C1", underline="single")
_BOLD = Font(bold=True)

COLUMNS = [
    "▶ POS", "▶ END",
    "CHROM", "POS", "END", "SVTYPE", "SVLEN", "FILTER", "QUAL", "GT",
    "PR_REF", "PR_ALT", "SR_REF", "SR_ALT", "SUPPORT",
    "USER_CLASS", "USER_ANNOT", "USER_COM",
]


# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------

def _open_vcf(path: Path):
    """Open VCF, plain or bgzip (bgzip is gzip-backward-compatible)."""
    if path.suffix in (".gz", ".bgz"):
        return gzip.open(path, "rt")
    return open(path)


def _info(info_str: str) -> dict:
    out = {}
    for item in info_str.split(";"):
        k, _, v = item.partition("=")
        out[k] = v if v else True
    return out


def _fmt(keys: str, vals: str) -> dict:
    ks = keys.split(":")
    vs = vals.split(":")
    return dict(zip(ks, vs + [""] * max(0, len(ks) - len(vs))))


def _split2(val: str) -> tuple[int, int]:
    parts = val.split(",") if val and val != "." else []
    try:
        return int(parts[0] or 0), int(parts[1] if len(parts) > 1 else 0)
    except (ValueError, IndexError):
        return 0, 0


def _detect_caller(headers: list[str]) -> str:
    for h in headers:
        if h.startswith("##source="):
            src = h.lower()
            if "manta" in src:
                return "manta"
            if "delly" in src:
                return "delly"
    return "manta"  # assume Manta by default


def parse_vcf(path: Path) -> tuple[str, list[dict]]:
    """Parse a SV VCF. Returns (caller_name, [record_dicts])."""
    headers: list[str] = []
    records: list[dict] = []

    with _open_vcf(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##"):
                headers.append(line)
                continue
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, filt, info_str = cols[:8]
            fmt_keys = cols[8] if len(cols) > 8 else ""
            sample_str = cols[9] if len(cols) > 9 else ""

            nfo = _info(info_str)
            fmt = _fmt(fmt_keys, sample_str) if fmt_keys else {}

            svtype = nfo.get("SVTYPE", ".")
            raw_svlen = nfo.get("SVLEN", "")
            raw_end = nfo.get("END", pos)

            try:
                svlen = abs(int(str(raw_svlen).split(",")[0]))
            except (ValueError, TypeError):
                svlen = None
            try:
                end = int(raw_end)
            except (ValueError, TypeError):
                end = int(pos)

            records.append({
                "CHROM": normalize_chr(chrom),
                "POS": int(pos),
                "ID": vid,
                "REF": ref,
                "ALT": alt,
                "QUAL": qual,
                "FILTER": filt,
                "SVTYPE": svtype,
                "SVLEN": svlen,
                "END": end,
                "GT": fmt.get("GT", "./."),
                "MATEID": nfo.get("MATEID", ""),
                "IMPRECISE": "IMPRECISE" in nfo,
                "CONTIG": bool(nfo.get("CONTIG")),
                "HOMLEN": int(nfo["HOMLEN"]) if str(nfo.get("HOMLEN", "")).isdigit() else 0,
                "DEPTH": int(nfo["BND_DEPTH"]) if str(nfo.get("BND_DEPTH", "")).isdigit() else None,
                "_fmt": fmt,
            })

    caller = _detect_caller(headers)

    for r in records:
        f = r.pop("_fmt")
        if caller == "delly":
            r["PR_REF"] = int(f.get("DR") or 0)
            r["PR_ALT"] = int(f.get("DV") or 0)
            r["SR_REF"] = int(f.get("RR") or 0)
            r["SR_ALT"] = int(f.get("RV") or 0)
        else:  # manta / unknown
            r["PR_REF"], r["PR_ALT"] = _split2(f.get("PR", ""))
            r["SR_REF"], r["SR_ALT"] = _split2(f.get("SR", ""))
        r["FT"] = f.get("FT", "")
        r["SUPPORT"] = r["PR_ALT"] + r["SR_ALT"]

    return caller, records


# ---------------------------------------------------------------------------
# Excel writing
# ---------------------------------------------------------------------------

def _write_data_sheet(ws, records: list[dict], port: int) -> None:
    for ci, name in enumerate(COLUMNS, 1):
        c = ws.cell(row=1, column=ci, value=name)
        c.font = _BOLD

    ws.freeze_panes = "C2"
    ws.column_dimensions["A"].width = 8
    ws.column_dimensions["B"].width = 8
    ws.column_dimensions["C"].width = 7   # CHROM
    ws.column_dimensions["F"].width = 6   # SVTYPE
    ws.column_dimensions["H"].width = 10  # FILTER

    for ri, rec in enumerate(records, 2):
        fill = _SVTYPE_FILL.get(rec["SVTYPE"])
        if rec["FILTER"] not in ("PASS", ".", ""):
            fill = _LOWQUAL_FILL

        locus_pos = build_locus(rec["CHROM"], rec["POS"])
        ca = ws.cell(row=ri, column=1, value="▶")
        ca.hyperlink = build_goto_url(locus_pos, port=port)
        ca.font = _LINK_FONT

        locus_end = build_locus(rec["CHROM"], rec["END"])
        cb = ws.cell(row=ri, column=2, value="▶")
        cb.hyperlink = build_goto_url(locus_end, port=port)
        cb.font = _LINK_FONT

        row_data = [
            rec["CHROM"], rec["POS"], rec["END"], rec["SVTYPE"],
            rec["SVLEN"], rec["FILTER"], rec["QUAL"], rec["GT"],
            rec["PR_REF"], rec["PR_ALT"], rec["SR_REF"], rec["SR_ALT"], rec["SUPPORT"],
            "", "", "",
        ]
        for ci, val in enumerate(row_data, 3):
            cell = ws.cell(row=ri, column=ci, value=val)
            if fill:
                cell.fill = fill


JOINT_COLUMNS = [
    "▶ POS", "▶ END",
    "CHROM", "POS", "END", "SVTYPE", "SVLEN",
    "SCORE", "LECTURE", "COMMENTAIRE",
    "GERMLINE", "SOMATIC", "SEUIL_ILLUMINA", "FLAGS",
    "A_QUAL", "A_GT", "A_FT", "A_FILTER", "B_FILTER",
    "PR_REF", "PR_ALT", "SR_REF", "SR_ALT", "SUPPORT", "AF_PCT", "DEPTH", "EVID_READS",
    "ID",
    "USER_CLASS", "USER_ANNOT", "USER_COM",
]

_HET_FILL = PatternFill("solid", fgColor="C6E0B4")      # vert — génotypé 0/1 ou 1/1
_MOSAIC_FILL = PatternFill("solid", fgColor="FFE699")   # ambre — reads alt sans génotype
_NOISE_FILL = PatternFill("solid", fgColor="F2F2F2")    # gris — support insuffisant

MIN_MOSAIC_SUPPORT = 3          # reads alternatifs (PR+SR) sous lesquels on ne lit rien
def _read_verdict(gt: str, support: int) -> str:
    """Lecture croisée des deux passes pour un même SV.

    La passe A génotype (fraction allélique supposée 0,5) ; la passe B ne
    génotype pas et ne filtre donc rien. Un SV que A appelle 0/1 est un
    hétérozygote candidat ; un SV que A laisse en 0/0 alors que des reads
    alternatifs existent est exactement le cas que la passe tumor-only sert à
    rattraper — mosaïcisme possible.
    """
    if "1" in (gt or ""):
        return "hétérozygote probable"
    if support >= MIN_MOSAIC_SUPPORT:
        return "candidat mosaïque"
    return "support faible"


def evidence_counts(bam: Path) -> dict[str, int]:
    """Nombre de reads-preuve par SV, lu dans un evidence BAM de Manta.

    Chaque read porte un tag ZM listant les SV qu'il soutient, sous la forme
    ZM:Z:<ID>|<TYPE>[,<ID>|<TYPE>...] où TYPE vaut PR, SR ou SRM. Manta ne
    tagge que le breakend principal d'une jonction BND ; le report sur le
    MATEID est fait à l'écriture de la feuille.

    Renvoie un dict vide si samtools est absent ou le BAM illisible — la
    colonne reste alors vide plutôt que de bloquer la génération.
    """
    samtools = shutil.which("samtools")
    if samtools is None or not Path(bam).exists():
        return {}
    try:
        out = subprocess.run([samtools, "view", str(bam)],
                             capture_output=True, text=True, check=True).stdout
    except (subprocess.CalledProcessError, OSError):
        return {}
    counts: dict[str, int] = {}
    for line in out.splitlines():
        for field in line.split("\t"):
            if not field.startswith("ZM:Z:"):
                continue
            for item in field[5:].split(","):
                sv_id = item.split("|")[0]
                if sv_id:
                    counts[sv_id] = counts.get(sv_id, 0) + 1
    return counts


def find_evidence_bam(run_dir: Path) -> Path | None:
    """Localise l'evidence BAM d'un run — présent seulement si Manta a tourné
    avec --generateEvidenceBam."""
    hits = sorted((run_dir / "results" / "evidence").glob("evidence_*.bam"))
    return hits[0] if hits else None


def _alt_fraction(rec: dict) -> float | None:
    """Fraction allélique approximative du SV, en %.

    PR et SR comptent des choses différentes (paires vs reads splittés) : on
    retient le canal qui porte le plus de preuve alternative plutôt que de les
    additionner. Sert au tri — un mosaïque réel se distingue du bruit par sa
    fraction, pas par son nombre brut de reads à 1450x.
    """
    pairs = [(rec["PR_REF"], rec["PR_ALT"]), (rec["SR_REF"], rec["SR_ALT"])]
    ref, alt = max(pairs, key=lambda t: t[1])
    total = ref + alt
    return round(100 * alt / total, 2) if total else None


MIN_MOSAIC_AF = 1.0             # % — sous ce seuil, le bruit de fond domine à 1450x


def _call_status(filter_str: str) -> str:
    """PASS ou FAIL(raisons) — lecture directe du FILTER d'une passe."""
    if not filter_str or filter_str in ("PASS", "."):
        return "PASS"
    return "FAIL(" + filter_str.replace(";", ",") + ")"


def _flags(rec: dict) -> str:
    """Drapeaux techniques que ni le FILTER germinal ni le somatique ne résument.

    Ce sont eux qui séparent un vrai événement d'un artefact de bord de cible :
    un SV assemblé (CONTIG) avec des reads splittés est solide ; un SV imprécis,
    sans paire, en microhomologie longue et en région MQ0 ne l'est pas.
    """
    out = []
    if rec.get("CONTIG"):
        out.append("assemblé")
    if rec.get("IMPRECISE"):
        out.append("imprécis")
    if rec.get("HOMLEN", 0) >= 10:
        out.append(f"homologie {rec['HOMLEN']}pb")
    if rec["SR_ALT"] == 0 and rec["PR_ALT"] > 0:
        out.append("paires seules")
    if rec["PR_ALT"] == 0 and rec["SR_ALT"] > 0:
        out.append("splits seuls")
    return ", ".join(out)


def _score(rec: dict, af: float | None, support: int, filt_a: str) -> int:
    """Score 0-100 de confiance dans le SV, indépendant du FILTER germinal.

    Le FILTER de la passe A rejette tout ce qui n'est pas à 50 % de fraction
    allélique — inutilisable pour trier des mosaïques. Le score reconstruit donc
    la confiance à partir des preuves brutes : fraction, profondeur de preuve,
    qualité de l'assemblage et du mapping.
    """
    s = 0
    af = af or 0
    s += 35 if af >= 20 else 28 if af >= 10 else 20 if af >= 5 else 8 if af >= MIN_MOSAIC_AF else 0
    s += 25 if support >= 20 else 18 if support >= 10 else 10 if support >= 5 else 4 if support >= 3 else 0
    if rec.get("CONTIG"):
        s += 15
    if rec.get("IMPRECISE"):
        s -= 10
    if rec["PR_ALT"] > 0 and rec["SR_ALT"] > 0:
        s += 15                                   # paires ET splits : deux preuves
    if rec.get("HOMLEN", 0) >= 10:
        s -= 8
    if "MaxMQ0Frac" in filt_a:
        s -= 20                                   # région à mapping ambigu
    if "NoPairSupport" in filt_a:
        s -= 8
    if "1" in (rec.get("GT") or ""):
        s += 10
    return max(0, min(100, s))


def _comment(gt: str, af: float | None, support: int, score: int, rec: dict) -> str:
    """Phrase de lecture pour le relecteur — ce que la ligne dit, en clair."""
    if "1" in (gt or ""):
        txt = f"Génotypé {gt} par la passe germinale — hétérozygote candidat (AF {af} %)."
        ft = rec.get("FT") or ""
        if ft and ft != "PASS":
            txt += f" Génotype de qualité basse (FT={ft})."
        return txt
    if support < MIN_MOSAIC_SUPPORT or (af or 0) < MIN_MOSAIC_AF:
        return "Trop peu de preuve alternative pour être lu — bruit de fond attendu à cette profondeur."
    base = (f"Écarté par la passe germinale (fraction {af} % ≠ 0,5) mais {support} read(s) "
            f"alternatif(s) présents — mosaïcisme possible.")
    if score >= 60:
        return base + " Preuve technique solide, à vérifier sur IGV."
    if "MaxMQ0Frac" in (rec.get("FILTER") or ""):
        return base + " Région à mapping ambigu (MQ0) — prudence."
    if rec.get("IMPRECISE"):
        return base + " Jonction non assemblée — position approximative."
    return base


ILLUMINA_TUMOR_ONLY_MIN = 15    # userGuide "Unpaired tumor sample" : PR[1] >= 15 || SR[1] >= 15


def _illumina_threshold(rec: dict) -> str:
    """Seuil de précision recommandé par Illumina en analyse tumor-only.

    Le userGuide de Manta ne fournit aucun modèle de score en mode non apparié
    et propose à la place ce filtre de départ sur les compteurs de l'allèle ALT.
    Colonne séparée du SCORE : c'est le critère de l'éditeur, pas le nôtre.
    """
    if rec["PR_ALT"] >= ILLUMINA_TUMOR_ONLY_MIN or rec["SR_ALT"] >= ILLUMINA_TUMOR_ONLY_MIN:
        return "PASS"
    return "—"


def _evidence_for(rec: dict, evidence: dict[str, int] | None) -> int | None:
    """Reads-preuve du SV, en reportant le compte du breakend principal.

    Manta n'écrit le tag ZM que sur un seul des deux enregistrements d'une
    jonction BND ; sans ce report, la moitié des lignes afficherait 0 à tort.
    """
    if not evidence:
        return None
    return evidence.get(rec["ID"]) or evidence.get(rec.get("MATEID", "")) or 0


def _write_joint_sheet(ws, recs_a: list[dict], recs_b: list[dict], port: int,
                       evidence: dict[str, int] | None = None) -> None:
    """Feuille SV : une ligne par variant, les deux passes côte à côte.

    Jointure sur l'ID Manta — les deux passes partagent les mêmes candidats et
    donc les mêmes identifiants. Un SV vu par une seule passe est conservé, avec
    les colonnes de l'autre vides.
    """
    for ci, name in enumerate(JOINT_COLUMNS, 1):
        ws.cell(row=1, column=ci, value=name).font = _BOLD
    ws.freeze_panes = "C2"
    for col, w in (("A", 8), ("B", 8), ("C", 7), ("F", 6), ("H", 7),
                   ("I", 21), ("J", 62), ("K", 24), ("L", 16), ("M", 26)):
        ws.column_dimensions[col].width = w

    by_id_b = {r["ID"]: r for r in recs_b}
    seen_b: set[str] = set()
    merged: list[tuple[dict | None, dict | None]] = []
    for a in recs_a:
        b = by_id_b.get(a["ID"])
        if b is not None:
            seen_b.add(a["ID"])
        merged.append((a, b))
    merged.extend((None, b) for b in recs_b if b["ID"] not in seen_b)

    ri = 1
    for a, b in merged:
        ri += 1
        ref = a or b
        support = max(a["SUPPORT"] if a else 0, b["SUPPORT"] if b else 0)
        verdict = _read_verdict(a["GT"] if a else "", support)
        fill = {
            "hétérozygote probable": _HET_FILL,
            "candidat mosaïque": _MOSAIC_FILL,
        }.get(verdict, _NOISE_FILL)

        for col, pos in ((1, ref["POS"]), (2, ref["END"])):
            c = ws.cell(row=ri, column=col, value="▶")
            c.hyperlink = build_goto_url(build_locus(ref["CHROM"], pos), port=port)
            c.font = _LINK_FONT

        src = a or b
        af = _alt_fraction(src)
        filt_a = a["FILTER"] if a else ""
        score = _score(src, af, support, filt_a)
        row = [
            ref["CHROM"], ref["POS"], ref["END"], ref["SVTYPE"], ref["SVLEN"],
            score, verdict,
            _comment(a["GT"] if a else "", af, support, score, src),
            _call_status(filt_a) if a else "absent",
            _call_status(b["FILTER"]) if b else "absent",
            _illumina_threshold(src),
            _flags(src),
            a["QUAL"] if a else "", a["GT"] if a else "", a["FT"] if a else "",
            filt_a, b["FILTER"] if b else "",
            src["PR_REF"], src["PR_ALT"], src["SR_REF"], src["SR_ALT"], support,
            af, src.get("DEPTH"), _evidence_for(src, evidence),
            ref["ID"],
            "", "", "",
        ]
        for ci, val in enumerate(row, 3):
            cell = ws.cell(row=ri, column=ci, value=val)
            if fill:
                cell.fill = fill

    last_col = get_column_letter(len(JOINT_COLUMNS))
    ws.auto_filter.ref = f"C1:{last_col}{ri}"


def _write_summary(ws, sources: dict, tracks: list, port: int, genome: str) -> None:
    for ci, h in enumerate(["Source", "Total", "PASS", "LowQual/autres"], 1):
        ws.cell(row=1, column=ci, value=h).font = _BOLD

    for ri, (name, recs) in enumerate(sources.items(), 2):
        n_pass = sum(1 for r in recs if r["FILTER"] == "PASS")
        n_low = sum(1 for r in recs if r["FILTER"] not in ("PASS", ".", ""))
        ws.cell(row=ri, column=1, value=name)
        ws.cell(row=ri, column=2, value=len(recs))
        ws.cell(row=ri, column=3, value=n_pass)
        ws.cell(row=ri, column=4, value=n_low)

    if tracks:
        url = build_load_url(tracks, locus=None, port=port, genome=genome)
        r = len(sources) + 3
        labels = ", ".join(label for _, label in tracks)
        cell = ws.cell(row=r, column=1, value=f"▶ Charger les BAM ({labels})")
        cell.hyperlink = url
        cell.font = _LINK_FONT
        for i, (path, label) in enumerate(tracks):
            ws.cell(row=r + 1 + i, column=1, value=label)
            ws.cell(row=r + 1 + i, column=2, value=str(path))

    ws.column_dimensions["A"].width = 28
    ws.column_dimensions["B"].width = 8
    ws.column_dimensions["C"].width = 8
    ws.column_dimensions["D"].width = 14


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# Manta standard output filenames inside results/variants/
_MANTA_VCFS = [
    ("diploidSV",          "diploidSV.vcf.gz"),
    ("tumorSV",            "tumorSV.vcf.gz"),      # passe tumor-only (mosaïque)
    ("candidateSV",        "candidateSV.vcf.gz"),
    ("candidateSmallIndels", "candidateSmallIndels.vcf.gz"),
]


def _collect_tracks(args, anchor: Path) -> list[tuple[Path, str]]:
    """Pistes IGV du _summary : BAM patient + un 2e BAM de comparaison.

    Le 2e BAM est soit --ref-bam, soit celui de la piste de référence
    (igv.reference_lane) trouvée dans igv.bam_dir — mêmes règles que
    igv_goto : un bam_dir relatif est ancré sur le dossier du XLSX, et on
    retombe sur un autre BAM du run si la lane est absente.
    """
    tracks: list[tuple[Path, str]] = []
    if args.bam and Path(args.bam).exists():
        bam = Path(args.bam).resolve()
        tracks.append((bam, bam.stem))

    if args.ref_bam:
        ref = Path(args.ref_bam).expanduser()
        if ref.exists():
            tracks.append((ref.resolve(), f"REF_{ref.stem}"))
        return tracks

    if not args.bam:
        return tracks                       # pas de patient → pas de piste de comparaison

    cfg = {}
    if args.config and Path(args.config).exists():
        with open(args.config) as f:
            cfg = (yaml.safe_load(f) or {}).get("igv", {}) or {}
    bam_dir_cfg = args.bam_dir or cfg.get("bam_dir")
    if not bam_dir_cfg:
        return tracks
    bam_dir = resolve_bam_dir(bam_dir_cfg, anchor)
    lane = args.reference_lane or cfg.get("reference_lane", "S2")
    ref = resolve_reference_bam(bam_dir, lane, Path(args.bam).stem if args.bam else "")
    if ref is not None:
        tracks.append((ref.resolve(), reference_label(ref, lane)))
    return tracks


def _main_joint(args, ap):
    """Deux passes Manta → un seul classeur, feuille SV jointe.

    Passe A (--run-dir) fournit diploidSV + les candidats non scorés ;
    passe B (--run-dir-b) fournit tumorSV. Les candidateSV/candidateSmallIndels
    sont identiques entre les deux passes : un seul exemplaire est écrit.
    """
    var_a = args.run_dir / "results" / "variants"
    var_b = args.run_dir_b / "results" / "variants"
    vcf_a = var_a / "diploidSV.vcf.gz"
    vcf_b = var_b / "tumorSV.vcf.gz"
    for v in (vcf_a, vcf_b):
        if not v.exists():
            ap.error(f"VCF introuvable : {v}")

    out = args.output or (var_a / f"{args.bam.stem if args.bam else 'manta'}_sv.xlsx")
    out.parent.mkdir(parents=True, exist_ok=True)

    _, recs_a = parse_vcf(vcf_a)
    _, recs_b = parse_vcf(vcf_b)
    print(f"  diploidSV {len(recs_a)} / tumorSV {len(recs_b)}", end=" ", flush=True)

    extras: dict[str, list[dict]] = {}
    for sheet, fname in (("candidateSV", "candidateSV.vcf.gz"),
                         ("candidateSmallIndels", "candidateSmallIndels.vcf.gz")):
        path = var_a / fname
        if path.exists():
            extras[sheet] = parse_vcf(path)[1]

    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    sources = {"diploidSV (passe A)": recs_a, "tumorSV (passe B)": recs_b}
    sources.update(extras)
    _write_summary(wb.create_sheet("_summary"), sources,
                   _collect_tracks(args, out.parent), args.port, args.genome)

    ev_bam = args.evidence_bam or find_evidence_bam(args.run_dir)
    evidence = evidence_counts(ev_bam) if ev_bam else {}
    if evidence:
        print(f"[evidence] {ev_bam.name} : {len(evidence)} jonction(s) taggée(s)")
    _write_joint_sheet(wb.create_sheet("SV"), recs_a, recs_b, args.port, evidence)
    for sheet, recs in extras.items():
        _write_data_sheet(wb.create_sheet(sheet), recs, args.port)

    wb.save(out)
    n_het = sum(1 for r in recs_a if "1" in (r["GT"] or ""))
    ids_a = {r["ID"] for r in recs_a}
    n_mos = sum(1 for r in recs_a
                if "1" not in (r["GT"] or "") and r["SUPPORT"] >= MIN_MOSAIC_SUPPORT)
    n_only_b = sum(1 for r in recs_b if r["ID"] not in ids_a)
    print(f"→ {out.name} : {n_het} htz, {n_mos} mosaïques, {n_only_b} vus par B seule")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="KODEX — SV VCF → Excel interprétable avec liens IGV"
    )
    ap.add_argument("vcfs", nargs="*", metavar="VCF",
                    help="VCF(s) positionnels (mode manuel)")
    ap.add_argument("--run-dir", type=Path, default=None,
                    help="Dossier runWorkflow Manta (auto-découverte des VCFs)")
    ap.add_argument("--run-dir-b", type=Path, default=None,
                    help="Dossier runWorkflow de la passe tumor-only — active la feuille SV jointe")
    ap.add_argument("--bam", type=Path, default=None,
                    help="BAM patient à charger dans IGV")
    ap.add_argument("--evidence-bam", type=Path, default=None,
                    help="evidence BAM Manta (--generateEvidenceBam) — sinon cherché dans le run-dir")
    ap.add_argument("--ref-bam", type=Path, default=None,
                    help="2e BAM (comparaison) — sinon déduit du config.yaml (igv.bam_dir + igv.reference_lane)")
    ap.add_argument("--bam-dir", type=Path, default=None,
                    help="Dossier des BAMs (défaut: igv.bam_dir du config)")
    ap.add_argument("--reference-lane", default=None,
                    help="Suffixe de la piste de référence (défaut: igv.reference_lane du config)")
    ap.add_argument("--config", type=Path, default=Path(__file__).parent / "config.yaml",
                    help="config.yaml KODEX (bloc igv:)")
    ap.add_argument("--output", "-o", type=Path, default=None,
                    help="Fichier .xlsx de sortie")
    ap.add_argument("--port", type=int, default=60151)
    ap.add_argument("--genome", default="hg38")
    args = ap.parse_args(argv)

    if args.run_dir_b:
        return _main_joint(args, ap)

    # Collect (sheet_name, path) pairs
    vcf_pairs: list[tuple[str, Path]] = []
    if args.run_dir:
        var_dir = args.run_dir / "results" / "variants"
        for sheet, fname in _MANTA_VCFS:
            p = var_dir / fname
            if p.exists():
                vcf_pairs.append((sheet, p))
    for s in args.vcfs:
        p = Path(s)
        stem = p.name.replace(".vcf.gz", "").replace(".vcf", "")
        vcf_pairs.append((stem, p))

    if not vcf_pairs:
        ap.error("Aucun VCF trouvé. Utiliser --run-dir ou fournir des VCFs positionnels.")

    # Default output path
    if args.output:
        out = args.output
    elif args.run_dir:
        stem = args.bam.stem if args.bam else "manta"
        out = args.run_dir / "results" / "variants" / f"{stem}_sv.xlsx"
    else:
        out = vcf_pairs[0][1].parent / f"{vcf_pairs[0][0]}_sv.xlsx"

    # Parse
    sources: dict[str, list[dict]] = {}
    for sheet, path in vcf_pairs:
        print(f"  {path.name} ...", end=" ", flush=True)
        caller, recs = parse_vcf(path)
        print(f"{len(recs)} variant(s) [{caller}]")
        sources[sheet] = recs

    # Write workbook
    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    _write_summary(wb.create_sheet("_summary"), sources,
                   _collect_tracks(args, out.parent), args.port, args.genome)

    for sheet, recs in sources.items():
        ws = wb.create_sheet(sheet)
        _write_data_sheet(ws, recs, args.port)
        if not recs:
            ws.cell(row=2, column=3, value="(aucun variant appelé dans cette run)")

    wb.save(out)
    total = sum(len(v) for v in sources.values())
    print(f"\n→ {out}")
    print(f"  {total} variant(s) — {len(sources)} feuille(s)")


if __name__ == "__main__":
    main()
