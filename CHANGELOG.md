# Changelog

Toutes les modifications notables de KODEX sont documentées ici.
Format : [Keep a Changelog](https://keepachangelog.com/fr/1.0.0/) — versioning [SemVer](https://semver.org/).

---

## [1.0.0] — 2026-03-27

### Ajouté
- `vcf_filter.py` — filtrage multi-critères des variants NiourK (gnomAD, RNCNT, PJCNT, IMPACT, ClinVar, LOVD, SpliceAI)
- Feuilles de sortie : `HIGH_IMPACT`, `CLINVAR_P_LP`, `FQ_MODERATE`, `LOVD`, `MOSAIC`, `RHEB_PKD1`, `_summary`, `_annotations`
- Score de qualité par variant (3 axes : technique, biologique, rareté) → dégradé couleur sur cellule Import DefGen
- Marquage PKD1/RHEB : fond bleu + commentaire pour les `downstream_gene_variant` dans ces gènes (GRCh38)
- Alerte discordance AF/run : HGVSc rouge si variant récurrent HTZ ici en mosaïque (ou inversement)
- Déduplication inter-feuilles : HGVSc bleu + commentaire si présent dans plusieurs filtres
- Formules `USER_COM` synchronisées entre feuilles (`=MID(TEXTJOIN(...))`)
- `aggregate_annotations.py` — consolidation multi-échantillons en un seul XLSX
- `igv_capture.py` — captures IGV automatisées (3 vues : strand, softclip, squished)
- `igv_bookmark.py` — injection des variants annotés comme ROI dans IGV + export BED
- `config.yaml` — configuration centralisée (filtres, scoring, exclusions IGV, discordance AF)
- `environment.yml`, `install.sh`, `run.sh` — packaging Ubuntu/micromamba
- Suite de tests pytest (56 tests)
- Documentation : `docs/pipeline_vcf_filter.md`, `docs/NOTICE_UTILISATION.md`
- Renommage du projet en **KODEX** (lignée K — CHU Angers)

---

## [Unreleased]

- Import des variants classifiés dans la base MySQL `tsc_db.variants`
- Classification ACMG/ClinGen intégrée
- Lien avec VCF_Interpreter pour la revue interactive
