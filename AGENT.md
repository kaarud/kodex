# AGENT.md — guide du dépôt KODEX pour un agent

Ce fichier sert à un agent (ou à un humain pressé) qui veut comprendre vite ce
que fait le dépôt et où se trouve chaque fonction. Le mode d'emploi détaillé est
dans `README.md` ; les choix de conception sont dans `docs/superpowers/specs/`.

## Ce que fait KODEX

Post-traitement des sorties du pipeline NiourK (CHU Angers, Nextflow) pour la
réanalyse NGS court-read TSC1/TSC2 (panel Twist TSC, Illumina NextSeq, GRCh38).
Il prend un classeur NiourK de plusieurs milliers de variants et produit un
classeur filtré, scoré et annoté d'une dizaine de variants par feuille, plus des
vues de couverture (CNV) et des outils de vérification sur BAM.

**KODEX signale des variants à revoir ; il ne pose pas de diagnostic.** Une vue
se juge sur la visibilité qu'elle donne au signal, pas comme un caller
autonome. Toute alerte se vérifie sur IGV et par une méthode orthogonale.

## Carte des points d'entrée

| Fichier | Rôle | Entrée → sortie |
|---|---|---|
| `vcf_filter.py` | Cœur : chargement, exclusions, feuilles par filtre, score 3 axes, mise en forme XLSX | `*.xlsx` NiourK → `*_filtered.xlsx` |
| `aggregate_annotations.py` | Agrège les `*_filtered.xlsx` d'un run, feuille de synthèse avec liens IGV | dossier → `<RUN>_annotations.xlsx` |
| `run_background.py` | Bruit de fond du run à la position de chaque variant (contraste patient / run) | `*_filtered.xlsx` + BAMs → `*_bg.xlsx` (commentaires de cellule) |
| `bg_filter.py` | Relit ces commentaires et garde la signature du témoin positif | `*_bg.xlsx` → `bg_candidats.xlsx` |
| `minibam.py` | Tranches BAM légères autour des variants annotés | annotations + BAMs → mini-BAMs |
| `sv_xlsx.py` | SV Manta / Delly → Excel avec liens IGV et score de preuve | VCF SV → `.xlsx` |
| `igv_capture.py`, `igv_bookmark.py`, `igv_goto.py`, `igv_link.py` | Pilotage d'IGV Desktop (port 60151) : captures, régions d'intérêt, liens go-to ; `igv_link.py` est la logique pure | — |
| `rpkm_exon_curve.R` | KODEX CNV : RPKM par exon, une ligne par patient | classeur KoverGraph → PDF/PNG/HTML |
| `perbase_coverage_curve.R` | KODEX CNV : profondeur base par base, normalisations échantillon / cohorte, jonction TSC2-PKD1 | mosdepth `*.per-base.bed.gz` → rapport |
| `perbase_coverage_curve_memory_safe.R` | Même courbe en pleine résolution pour gros runs, sans HTML | idem |
| `coverage_common.R` | Fonctions R partagées : seuils DeCovA 0,8 / 1,2, dossiers horodatés, rendu | — |
| `run.sh`, `install.sh` | CLI et installation | — |

## Configuration

Deux profils de schéma NiourK, décrits en YAML et non en dur dans le code :
`config.yaml` (run 2, défaut) et `config_legacy.yaml` (run 1, `--legacy`).
Alias de colonnes, colonne gnomAD de filtrage, exclusions, poids du score,
colonnes exportées : tout passe par là.

## Pièges connus

- `perbase_coverage_curve.R` : `globalize_layout()` ne doit pas tourner sur une
  disposition génomique, sinon PKD1 passe avant TSC2 (protégé et testé).
- PKD1 exons 1-33 sont dans la duplication segmentaire 16p (pseudogènes
  PKD1P1-P6) : un ratio élevé y est un artefact de mapping, pas une duplication.
- Dans chaque exon, les bases sont placées par position génomique croissante :
  le profil intra-exon est inversé pour les gènes brin moins (TSC1).
- `bg_filter.py` : `max_with_alt_frac = 0,35` est calé sur le témoin positif de
  run 2 ; en dessous de `min_dp = 100` reads, le contraste n'est pas mesuré.
- Le filtrage run 2 utilise `GNOMAD4_AF` avec des seuils calibrés sur l'ancienne
  popmax gnomAD 2.1 (voir README).

## Tests

```bash
micromamba run -n kodex python -m pytest tests/ -q          # pytest, fixtures dans tests/conftest.py
cd tests && for f in *.R; do micromamba run -n kodex Rscript "$f"; done   # testthat, chemins relatifs à tests/
```

## Données

Aucune donnée patient dans le dépôt : `*.xlsx`, `*.vcf*`, `*.bam`, `data/`,
`output/` et `panel/` (géométrie du panel, dérivée des données NiourK) sont
ignorés par `.gitignore`. Ne jamais les committer. Aucun chemin de données
n'est codé en dur : BAM, profondeurs et classeurs se passent en argument.

Nommage des runs, aligné sur le manuscrit : **run 1** = export Niourk v2.5
(janvier 2026, profil `config_legacy.yaml`), **run 2** = export Niourk2 v1.0.0
sous Singularity (mai 2026, profil `config.yaml`).

## Où lire la suite

- `docs/superpowers/specs/` — conception de chaque fonctionnalité
- `docs/superpowers/plans/` — plans d'implémentation
- `docs/LECTURE_CLASSEUR.md` — lecture du classeur de sortie (feuilles, SCORE, couleurs, dépannage)
- `docs/NOTICE_CNV.md` — rapport HTML de couverture
- `docs/igv_bam_interpretation.md` — lecture des BAM dans IGV
- `TODO.md`, `CHANGELOG.md`
