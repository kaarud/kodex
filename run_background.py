"""run_background.py — KODEX : bruit de fond du run à la position d'un variant.

Répond à « est-ce que ce variant se voit chez les autres ? », que les colonnes
NiourK ne peuvent pas capturer : RNCNT_* et PJCNT_* comptent les échantillons où
un variant a été **appelé**, jamais le bruit qu'aucun caller n'a appelé. Un
artefact présent dans tout le run mais appelé chez un seul patient sort donc avec
RNCNT_TOTAL = 1, soit le score de rareté maximal.

Ce qui discrimine n'est pas la présence de reads alternatifs chez les autres —
sur une cohorte d'ADN dégradés la plupart en ont — mais le **contraste** entre la
fraction du patient et celle du run à cette position précise.

Spec : docs/superpowers/specs/2026-09-03-run-background-noise-design.md
"""

from statistics import median

_BASES = ("A", "C", "G", "T")
_EMPTY = {"depth": 0, "A": 0, "C": 0, "G": 0, "T": 0, "ins": 0, "del": 0}


# ---------------------------------------------------------------------------
# Logique pure — testable sans BAM ni samtools
# ---------------------------------------------------------------------------

def classify_variant(ref: str, alt: str) -> str:
    """Quel compteur du pileup lit ce variant : 'snv', 'ins', 'del' ou 'mnv'.

    'mnv' (substitution multi-nucléotidique) n'a pas de compteur unique dans un
    pileup par position — la mesure ne s'y applique pas.
    """
    ref, alt = str(ref).upper(), str(alt).upper()
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(alt) > len(ref) and alt.startswith(ref):
        return "ins"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "del"
    return "mnv"


def alt_fraction(counts: dict, ref: str, alt: str) -> float | None:
    """Fraction des reads portant l'allèle alternatif, ou None si non mesurable.

    None signifie « pas de mesure » (aucun read, ou variant sans compteur
    correspondant) — jamais « bruit nul ». La distinction compte : un site non
    couvert ne doit pas se lire comme un site propre.
    """
    depth = counts.get("depth", 0)
    if not depth:
        return None
    kind = classify_variant(ref, alt)
    if kind == "snv":
        return counts.get(str(alt).upper(), 0) / depth
    if kind in ("ins", "del"):
        return counts.get(kind, 0) / depth
    return None


def summarise_noise(
    patient: dict,
    others: dict[str, dict],
    ref: str,
    alt: str,
    genotype_floor: float = 0.25,
    min_other_samples: int = 5,
    statistic: str = "median",
) -> dict:
    """Résume le bruit du run à cette position et le contraste du patient.

    Les échantillons dont la fraction dépasse `genotype_floor` sont des porteurs,
    pas du bruit : les inclure écraserait le contraste vers 1 et condamnerait un
    variant partagé légitime. Ils sont exclus du calcul et listés dans `carriers`.

    Médiane plutôt que maximum : sur une cohorte hétérogène, un seul échantillon
    dégradé ne doit pas condamner un variant vrai chez tous les autres. Le maximum
    est rendu à part — l'écart médiane/max distingue « un échantillon pourri »
    d'un « site chroniquement mauvais ».

    Le contraste est borné par la profondeur : un bruit mesuré à 0 sur 500× ne
    vaut pas « nul » mais « sous 1/500 ».
    """
    patient_af = alt_fraction(patient, ref, alt)

    noise, carriers, depths = [], [], []
    for sample, counts in sorted(others.items()):
        af = alt_fraction(counts, ref, alt)
        if af is None:
            continue
        depths.append(counts.get("depth", 0))
        if af >= genotype_floor:
            carriers.append(sample)
            continue
        noise.append(af)

    out = {
        "patient_af": patient_af,
        "n_measured": len(noise),
        "n_with_alt": sum(1 for af in noise if af > 0),
        "carriers": carriers,
        "median": median(noise) if noise else None,
        "max": max(noise) if noise else None,
        "contrast": None,
    }

    if patient_af is None or len(noise) < min_other_samples:
        return out

    background = out["max"] if statistic == "max" else out["median"]
    floor = 1 / median(depths) if depths and median(depths) else 1.0
    out["contrast"] = patient_af / max(background, floor)
    return out


# ---------------------------------------------------------------------------
# Lecture des BAM
# ---------------------------------------------------------------------------

def primary_ratio(counts: dict, ref: str, alt: str) -> float | None:
    """AF calculee sans les alignements supplementaires, rapportee a l'AF brute.

    Un variant reel garde sa fraction quand on retire les reads `0x800` : ses
    porteurs sont des alignements primaires. Un artefact de structure la perd,
    ses porteurs etant des morceaux de reads recolles par l'aligneur.

    Le rapport echappe au biais de VAF : le *taux* de reads supplementaires est
    mecaniquement eleve a basse VAF (quelques reads parasites pesent lourd sur un
    petit denominateur), le rapport non. Valide sur le temoin positif du run
    (mosaique 1,2 %, rapport 1,02) contre des artefacts de VAF voisine (0,04-0,30).

    None si la mesure n'est pas disponible.
    """
    pri = counts.get("pri")
    if not pri:
        return None
    raw = alt_fraction(counts, ref, alt)
    fil = alt_fraction(pri, ref, alt)
    if raw is None or fil is None or raw == 0:
        return None
    return fil / raw


def dp_colour(contrast: float | None, ratio: float | None,
              cfg: dict | None = None) -> str:
    """Couleur de la cellule DP : verdict combine des deux detecteurs.

    Les deux modes d'echec sont independants et chacun aveugle a l'autre :
    un artefact partage par le run garde un rapport primaire de 1, un artefact
    de structure propre a un echantillon garde un contraste eleve. Le pire des
    deux verdicts l'emporte.
    """
    # contrast_good calé sur le témoin positif du run : la seule mosaïque connue
    # pour réelle ici mesure x7,3 ; un seuil à 10 l'aurait classée « à vérifier ».
    c = {"contrast_bad": 1.5, "contrast_weak": 3.0, "contrast_good": 5.0,
         "ratio_bad": 0.3, "ratio_weak": 0.8, **(cfg or {})}
    if contrast is None and ratio is None:
        return "grey"
    if (contrast is not None and contrast < c["contrast_bad"]) or \
       (ratio is not None and ratio < c["ratio_bad"]):
        return "red"
    if (contrast is not None and contrast < c["contrast_weak"]) or \
       (ratio is not None and ratio < c["ratio_weak"]):
        return "orange"
    if (contrast is None or contrast >= c["contrast_good"]) and \
       (ratio is None or ratio >= c["ratio_weak"]):
        return "green"
    return "yellow"


def cluster_positions(positions, gap: int = 2000):
    """Regroupe les positions en regions (chrom, start, end) triees.

    Lire un BAM position par position coute ~200 ms sur ce panel (9 000 reads a
    la position, 1,9 M sur 20 kb) ; un seul parcours de region coute 30 fois
    moins par read. Les variants TSC1/TSC2 sont groupes par exon, donc quelques
    regions suffisent. `gap` borne l'elargissement : au-dela, on parcourrait des
    milliers de colonnes inutiles.
    """
    from collections import defaultdict

    by_chrom = defaultdict(set)
    for chrom, pos in positions:
        by_chrom[str(chrom)].add(int(pos))

    out = []
    for chrom in sorted(by_chrom):
        ps = sorted(by_chrom[chrom])
        start = prev = ps[0]
        for pos in ps[1:]:
            if pos - prev > gap:
                out.append((chrom, start, prev))
                start = pos
            prev = pos
        out.append((chrom, start, prev))
    return out


def count_positions(bam, chrom: str, wanted: set, start: int, end: int,
                    min_mapq: int = 0, min_baseq: int = 0) -> dict:
    """Compteurs pour toutes les positions voulues d'une region, en un parcours."""
    got = {}
    ref = _match_contig(bam, chrom)
    if ref is None:
        return got
    for column in bam.pileup(ref, start - 1, end, truncate=True,
                             min_mapping_quality=min_mapq,
                             min_base_quality=min_baseq,
                             max_depth=1_000_000, stepper="all"):
        pos = column.reference_pos + 1
        if pos not in wanted:
            continue
        counts = dict(_EMPTY)
        pri = dict(_EMPTY)          # sans les alignements supplementaires
        for read in column.pileups:
            targets = (counts,) if read.alignment.is_supplementary else (counts, pri)
            for c in targets:
                c["depth"] += 1
                if read.indel < 0:
                    c["del"] += 1
                elif read.indel > 0:
                    c["ins"] += 1
            if read.is_del or read.query_position is None:
                continue
            base = (read.alignment.query_sequence or "")[read.query_position:
                                                         read.query_position + 1]
            if base.upper() in _BASES:
                for c in targets:
                    c[base.upper()] += 1
        counts["pri"] = pri
        counts["pos"] = pos          # d'ou vient reellement la mesure (cf. fenetre indel)
        got[pos] = counts
    return got


def count_indel_window(bam_path, chrom: str, pos: int, kind: str,
                       window: int = 10, **kw) -> dict:
    """Compteurs a la position ou l'indel est reellement aligne, dans +/- window.

    VCF et BAM ne left-alignent pas toujours un indel au meme endroit : sur
    TSC2 c.2220+1_2220+2del, le VCF ancre en chr16:2072363 et le BAM porte les
    16 deletions cinq bases en amont. Chercher a la seule position VCF renvoyait
    zero pour un variant reel.

    On retient la position de la fenetre qui porte le plus d'evenements du bon
    signe, la position VCF gagnant a egalite. La fenetre reste courte : au-dela,
    on capterait un indel voisin sans rapport.
    """
    best = count_alleles(bam_path, chrom, pos, **kw)
    if kind not in ("ins", "del"):
        return best
    for off in range(-window, window + 1):
        if off == 0:
            continue
        c = count_alleles(bam_path, chrom, pos + off, **kw)
        if c[kind] > best[kind]:
            best = c
    return best


def count_alleles(bam_path, chrom: str, pos: int,
                  min_mapq: int = 0, min_baseq: int = 0) -> dict:
    """Compteurs bruts à une position, façon popup IGV : A/C/G/T, ins, del, depth.

    Pas de filtre sur MAPQ ni sur la qualite de base : voir le bruit que les
    callers masquent est l'objectif, et comparer brut a filtre est en soi un
    signal (cas hairpin : VAF 1,8 % -> 0,09 % apres MAPQ).

    Duplicats, secondaires et supplementaires sont en revanche exclus
    (stepper="all"). Un duplicat PCR n'est pas une observation independante :
    les compter gonflait la profondeur d'un facteur ~15 sur ce panel, avec des
    reads massivement REF, ce qui divisait d'autant les fractions alternatives.

    Convention VCF pour les indels : ils sont ancres sur la base *precedente*.
    Une deletion se lit donc `read.indel < 0` a la position d'ancrage, pas
    `is_del` — qui teste si la base de cette colonne est elle-meme absente et
    renvoyait zero pour toutes les deletions.
    """
    import pysam

    counts = dict(_EMPTY)
    opened = not hasattr(bam_path, "pileup")
    bam = pysam.AlignmentFile(str(bam_path), "rb") if opened else bam_path
    try:
        chrom = _match_contig(bam, chrom)
        if chrom is None:
            return counts
        for column in bam.pileup(chrom, pos - 1, pos, truncate=True,
                                 min_mapping_quality=min_mapq,
                                 min_base_quality=min_baseq,
                                 max_depth=1_000_000, stepper="all"):
            for read in column.pileups:
                counts["depth"] += 1
                if read.indel < 0:
                    counts["del"] += 1
                elif read.indel > 0:
                    counts["ins"] += 1
                if read.is_del or read.query_position is None:
                    continue
                base = (read.alignment.query_sequence or "")[read.query_position:
                                                             read.query_position + 1]
                if base.upper() in _BASES:
                    counts[base.upper()] += 1
    finally:
        if opened:
            bam.close()
    return counts


def _match_contig(bam, chrom: str) -> str | None:
    """Concilie les conventions 'chr9' / '9' entre le XLSX et le BAM."""
    refs = set(bam.references)
    for candidate in (chrom, chrom.removeprefix("chr"), f"chr{chrom}"):
        if candidate in refs:
            return candidate
    return None


def sample_of(bam_path) -> str:
    """'100007_S6.bam' -> '100007_S6'. Les pistes dérivées ne sont pas des échantillons."""
    from pathlib import Path
    return Path(bam_path).name[:-len(".bam")]


def raw_bams(bam_dir) -> dict[str, "Path"]:
    """BAM d'alignement du run, par échantillon — .chim et .mutect2 exclus.

    Mutect2 réassemble localement et fabrique des reads qui n'existent pas dans
    l'alignement ; les compter comme du bruit de fond mélangerait deux questions.
    """
    from pathlib import Path

    out = {}
    for bam in sorted(Path(bam_dir).glob("*.bam")):
        stem = bam.name[:-len(".bam")]
        if stem.endswith(".chim") or stem.endswith(".mutect2"):
            continue
        out[stem] = bam
    return out


def count_matrix(variants, bams: dict, progress=None, indel_window: int = 10,
                 detail=None, **kw) -> dict:
    """{echantillon: {(chrom, pos): compteurs}} — un parcours par region et par BAM.

    `variants` : (chrom, pos, kind). Les positions voisines d'un indel sont
    ajoutees a la collecte (VCF et BAM ne left-alignent pas toujours pareil), et
    le meilleur compteur de la fenetre est retenu ensuite.

    `progress(i, n, echantillon, secs, cumul)` est appele deux fois par BAM :
    `secs=None` a l'ouverture — c'est le signal « je travaille sur celui-ci » —
    puis avec sa duree une fois fini. `detail` recoit la trace par region.
    """
    import time

    import pysam

    exact, needed = [], set()
    for chrom, pos, kind in variants:
        chrom, pos = str(chrom), int(pos)
        exact.append((chrom, pos, kind))
        needed.add((chrom, pos))
        if kind in ("ins", "del"):
            for off in range(-indel_window, indel_window + 1):
                needed.add((chrom, pos + off))

    regions = cluster_positions(needed)
    wanted_by_chrom = {}
    for chrom, pos in needed:
        wanted_by_chrom.setdefault(chrom, set()).add(pos)

    if detail:
        detail(f"{len(needed)} positions a lire, groupees en {len(regions)} regions")

    out = {}
    t_start = time.monotonic()
    for i, (sample, path) in enumerate(sorted(bams.items()), 1):
        if progress:
            progress(i, len(bams), sample, None, 0.0)
        t_bam = time.monotonic()
        raw = {}
        with pysam.AlignmentFile(str(path), "rb") as bam:
            for j, (chrom, start, end) in enumerate(regions, 1):
                if detail:
                    detail(f"      region {j}/{len(regions)} "
                           f"{chrom}:{start}-{end} ({end - start + 1} pb)")
                raw.setdefault(chrom, {}).update(
                    count_positions(bam, chrom, wanted_by_chrom[chrom],
                                    start, end, **kw))
        per_sample = {}
        for chrom, pos, kind in exact:
            counts = raw.get(chrom, {}).get(pos, dict(_EMPTY))
            if kind in ("ins", "del"):
                for off in range(-indel_window, indel_window + 1):
                    other = raw.get(chrom, {}).get(pos + off)
                    if other and other[kind] > counts[kind]:
                        counts = other
            per_sample[(chrom, pos)] = counts
        out[sample] = per_sample
        if progress:
            progress(i, len(bams), sample, time.monotonic() - t_bam,
                     time.monotonic() - t_start)
    return out


def format_comment(sample: str, chrom: str, pos: int, ref: str, alt: str,
                   patient: dict, others: dict, summary: dict) -> str:
    """Commentaire de cellule, dans l'esprit du popup IGV : le patient, puis le run."""
    def line(counts):
        return (f"DP={counts['depth']}  A:{counts['A']} C:{counts['C']} "
                f"G:{counts['G']} T:{counts['T']}  DEL:{counts['del']} INS:{counts['ins']}")

    pct = lambda v: "—" if v is None else f"{100 * v:.2f} %"
    # Un indel n'est pas toujours left-aligne pareil dans le VCF et dans le BAM :
    # count_matrix retient la position de la fenetre qui porte l'evenement. La
    # dire evite de croire qu'on a mesure la position affichee.
    used = patient.get("pos")
    where = f"{chrom}:{pos}" + (f" (lu en {chrom}:{used})"
                                if used is not None and used != pos else "")
    out = [f"{where} {ref}>{alt} — {sample} (patient)",
           f"  {line(patient)}",
           f"  -> alt {pct(summary['patient_af'])}", ""]

    if summary["n_measured"] == 0:
        out.append("Aucun autre echantillon mesurable a cette position.")
        return "\n".join(out)

    out.append(f"Run — {summary['n_measured']} autres echantillons")
    for s, counts in sorted(others.items()):
        af = alt_fraction(counts, ref, alt)
        if af is None or af >= 0.25:
            continue
        out.append(f"  {s:<18s} {line(counts)}  -> {pct(af)}")

    out += ["",
            f"  mediane {pct(summary['median'])} | max {pct(summary['max'])}"
            f" | avec alt : {summary['n_with_alt']}/{summary['n_measured']}"]
    if summary["carriers"]:
        out.append(f"  porteurs exclus du bruit : {', '.join(summary['carriers'])}")
    out.append("  contraste "
               + ("non mesure (trop peu d'echantillons)" if summary["contrast"] is None
                  else f"x{summary['contrast']:.1f}"))
    ratio = primary_ratio(patient, ref, alt)
    if ratio is not None:
        verdict = ("reads primaires — OK" if ratio >= 0.8 else
                   "ARTEFACT DE STRUCTURE probable" if ratio < 0.3 else
                   "a verifier")
        fil = alt_fraction(patient.get("pri", {}), ref, alt)
        out.append(f"  AF hors supplementaires : {pct(fil)}"
                   f"  (rapport {ratio:.2f} — {verdict})")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Annotation des classeurs
# ---------------------------------------------------------------------------

_VARIANT_COLS = ("CHR", "POS", "REF", "ALT")


def _fills():
    """Palette de la colonne DP — memes teintes que les colonnes de vcf_filter."""
    from openpyxl.styles import PatternFill
    solid = lambda rgb: PatternFill(start_color=rgb, end_color=rgb, fill_type="solid")
    return {"green":  solid("C6EFCE"),   # se detache du bruit, reads primaires
            "yellow": solid("FFEB9C"),   # au-dessus du bruit, site charge
            "orange": solid("FFC7A0"),   # marge faible sur un des deux axes
            "red":    solid("FFC7CE"),   # bruit du run ou artefact de structure
            "grey":   solid("D9D9D9")}   # non mesure — different de "mesure nulle"


_DP_FILLS = _fills()


def make_log(verbose: bool = False, stream=None):
    """Journal horodate et vide a chaque ligne.

    Sans `flush`, un pipe ou un `tee` retient la sortie par blocs de 4 ko : le
    script paraissait mort alors qu'il tournait. L'horodatage relatif dit d'un
    coup d'oeil si la derniere etape a pris dix secondes ou dix minutes — c'est
    la seule facon de distinguer « lent » de « bloque ».

    Les messages de niveau 1 ne sortent qu'en mode verbeux.
    """
    import sys
    import time

    t0 = time.monotonic()
    out = stream if stream is not None else sys.stderr

    def log(msg="", level=0):
        if level and not verbose:
            return
        d = int(time.monotonic() - t0)
        print(f"[{d // 60:02d}:{d % 60:02d}] {msg}", file=out, flush=True)

    return log


def _eta(done: int, total: int, elapsed: float) -> str:
    """Temps restant estime, extrapole sur la moyenne des BAM deja traites."""
    if done <= 0 or done >= total:
        return ""
    left = elapsed / done * (total - done)
    return f" — reste ~{int(left) // 60:d} min {int(left) % 60:02d} s"


def annotate_workbooks(xlsx_paths, bam_dir, out_suffix="_bg",
                       extra_genes=("RHEB", "PKD1"), log=None, verbose=False, **kw):
    """Ecrit le bruit de fond du run en commentaire sur la cellule DP.

    Une copie `<nom><out_suffix>.xlsx` est produite : les classeurs d'origine
    portent les annotations de relecture, on n'ecrit jamais dedans. Le suffixe
    tombe hors du glob `*_filtered.xlsx` de l'agregation.

    Les positions sont collectees sur l'ensemble des classeurs avant le moindre
    pileup : un site partage par plusieurs patients n'est lu qu'une fois par BAM.
    """
    from pathlib import Path
    from openpyxl import load_workbook
    from openpyxl.comments import Comment

    from igv_capture import extract_sample_id

    xlsx_paths = [Path(p) for p in xlsx_paths]
    excluded = {g.upper() for g in extra_genes}
    if log is None:
        log = make_log(verbose)

    # Le detail passe par le meme journal, avec un niveau : c'est `make_log` qui
    # decide de l'ecrire ou non. Un `log` fourni par un appelant peut n'accepter
    # qu'un argument — on retombe alors sur un silence, pas sur une erreur.
    def detail(msg):
        try:
            log(msg, 1)
        except TypeError:
            pass

    log(f"1/3 collecte — {len(xlsx_paths)} classeur(s)")
    bams = raw_bams(bam_dir)
    if not bams:
        raise SystemExit(f"Aucun BAM dans {bam_dir}")

    # 1. collecte des positions, tous classeurs confondus
    wanted, books = [], {}
    for path in xlsx_paths:
        detail(f"      lecture {path.name}")
        wb = load_workbook(path)
        books[path] = wb
        for ws in wb.worksheets:
            if ws.title.startswith("_"):
                continue
            idx = {c.value: c.column for c in ws[1]}
            if not all(k in idx for k in _VARIANT_COLS) or "DP" not in idx:
                continue
            for row in range(2, ws.max_row + 1):
                sym = str(ws.cell(row, idx["SYMBOL"]).value or "").upper() if "SYMBOL" in idx else ""
                if sym in excluded:
                    continue
                chrom = ws.cell(row, idx["CHR"]).value
                pos = ws.cell(row, idx["POS"]).value
                ref = ws.cell(row, idx["REF"]).value
                alt = ws.cell(row, idx["ALT"]).value
                if not chrom or pos in (None, ""):
                    continue
                wanted.append((str(chrom), int(pos),
                               classify_variant(str(ref), str(alt))))
    n_pos = len(set((c, p) for c, p, _ in wanted))
    log(f"2/3 pileup — {n_pos} positions uniques x {len(bams)} echantillons")

    def on_bam(i, n, sample, secs, cumul):
        if secs is None:
            log(f"  [{i:>2}/{n}] {sample} ...")
        else:
            log(f"  [{i:>2}/{n}] {sample} — {secs:.0f} s{_eta(i, n, cumul)}")

    # 2. un seul passage par BAM
    matrix = count_matrix(wanted, bams, progress=on_bam, detail=detail, **kw)

    # 3. ecriture
    log("3/3 ecriture des commentaires")
    written = 0
    for path, wb in books.items():
        file_sample = extract_sample_id(path)
        for ws in wb.worksheets:
            if ws.title.startswith("_"):
                continue
            idx = {c.value: c.column for c in ws[1]}
            if not all(k in idx for k in _VARIANT_COLS) or "DP" not in idx:
                continue
            for row in range(2, ws.max_row + 1):
                sym = str(ws.cell(row, idx["SYMBOL"]).value or "").upper() if "SYMBOL" in idx else ""
                if sym in excluded:
                    continue
                chrom = ws.cell(row, idx["CHR"]).value
                pos = ws.cell(row, idx["POS"]).value
                if not chrom or pos in (None, ""):
                    continue
                key = (str(chrom), int(pos))
                ref = str(ws.cell(row, idx["REF"]).value)
                alt = str(ws.cell(row, idx["ALT"]).value)
                # Un agregat porte plusieurs patients dans une meme feuille :
                # seule la colonne SAMPLE dit de qui est la ligne.
                sample = (str(ws.cell(row, idx["SAMPLE"]).value).strip()
                          if "SAMPLE" in idx else file_sample)
                patient = matrix.get(sample, {}).get(key)
                if patient is None:
                    continue
                others = {s: c[key] for s, c in matrix.items()
                          if s != sample and key in c}
                summary = summarise_noise(patient, others, ref, alt)
                text = format_comment(sample, str(chrom), int(pos), ref, alt,
                                      patient, others, summary)
                cell = ws.cell(row, idx["DP"])
                # Taille calculee : un commentaire fixe tronquait la fin du texte
                # des que le run comptait plus d'une dizaine d'echantillons —
                # or la mediane, le contraste et le rapport primaire sont
                # justement les dernieres lignes.
                lines = text.split("\n")
                cell.comment = Comment(
                    text, "KODEX run_background",
                    width=min(1000, 30 + 7 * max(len(l) for l in lines)),
                    height=min(900, 30 + 14 * len(lines)))
                colour = dp_colour(summary["contrast"],
                                   primary_ratio(patient, ref, alt))
                if colour in _DP_FILLS:
                    cell.fill = _DP_FILLS[colour]
                written += 1
        out = path.with_name(f"{path.stem}{out_suffix}.xlsx")
        wb.save(out)
        log(f"  -> {out.name}")
    log(f"{written} cellules DP commentees")
    return written


def main(argv=None):
    import argparse
    from pathlib import Path

    ap = argparse.ArgumentParser(
        description="KODEX — bruit de fond du run en commentaire sur DP.")
    ap.add_argument("xlsx", type=Path, nargs="+", help="classeurs *_filtered.xlsx")
    ap.add_argument("--bam-dir", type=Path, required=True)
    ap.add_argument("--suffix", default="_bg",
                    help="suffixe du classeur de sortie (defaut: _bg)")
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--min-baseq", type=int, default=0)
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="detaille chaque classeur lu et chaque region parcourue")
    a = ap.parse_args(argv)
    annotate_workbooks(a.xlsx, a.bam_dir, out_suffix=a.suffix,
                       log=make_log(a.verbose), verbose=a.verbose,
                       min_mapq=a.min_mapq, min_baseq=a.min_baseq)


if __name__ == "__main__":
    main()
