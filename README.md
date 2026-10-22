# KODEX

**KODEX** est un plugin de post-traitement des variants NiourK pour la réanalyse NGS TSC1/TSC2.
Il déchiffre les fichiers XLSX produits par le pipeline NiourK (CHU Angers) et génère des tableaux filtrés, annotés et prêts pour la revue clinique.

> Nom : dans la lignée des outils CHU Angers (Koala, NiourK, KaveShare), KODEX évoque la lecture et le déchiffrement des signes génomiques (*augur* + *codex*).

---

## Fonctionnalités

- **Double schéma NiourK** — gère les deux générations d'export : **run 2** (Niourk2 v1.0.0, conteneur Singularity, mai 2026) par défaut, **run 1** (Niourk v2.5, janvier 2026) via `--legacy`. Les différences de schéma (renommage de colonnes, source gnomAD, colonnes USER_*) sont décrites en YAML, pas en dur dans le code.
- **Filtrage multi-critères** — gnomAD, RNCNT, PJCNT, IMPACT, CLINVAR, LOVD, SpliceAI
- **Score de qualité par variant** — 3 axes : technique, biologique, rareté → score combiné et couleur RGB sur la cellule `SCORE`, feuilles triées par score (lecture : [`docs/LECTURE_CLASSEUR.md`](docs/LECTURE_CLASSEUR.md))
- **Détection de la mosaïque** — feuille MOSAIC dédiée (AF 5–35%)
- **Marquage PKD1/RHEB** — variants `downstream_gene_variant` dans les coordonnées GRCh38 de PKD1 et RHEB → fond bleu + commentaire
- **Discordance AF/run** — alerte rouge sur HGVSc si un variant récurrent HTZ apparaît en mosaïque (ou inversement)
- **Déduplication inter-feuilles** — HGVSc bleu + commentaire si le variant apparaît dans plusieurs filtres
- **Formules USER_COM** — `=MID(TEXTJOIN(...))` pour synchroniser les annotations entre feuilles
- **Feuille `_summary`** — tableau de bord avec comptages par filtre
- **Feuille `_annotations`** — synthèse des variants déjà annotés par l'utilisateur
- **Agrégation cohorte** — `aggregate_annotations.py` consolide tous les `*_filtered.xlsx` en un seul fichier
- **IGV capture** — `igv_capture.py` génère 3 screenshots par variant (strand, softclip, squished)
- **IGV bookmark** — `igv_bookmark.py` injecte les variants annotés comme régions d'intérêt dans IGV
- **Bruit de fond du run** — `run_background.py` mesure, à la position de chaque variant, la fraction alternative chez les autres échantillons du run (contraste patient / run) ; `bg_filter.py` en extrait les variants qui ont la signature du témoin positif
- **Mini-BAM** — `minibam.py` découpe les BAM autour des variants retenus (quelques Mo au lieu de plusieurs Go, ouvrables dans IGV)
- **SV** — `sv_xlsx.py` exporte les VCF Manta / Delly en Excel avec liens IGV et score de preuve
- **KODEX CNV** (R) — courbes de couverture pour repérer délétions et duplications d'exons :
  - `rpkm_exon_curve.R` : courbe RPKM par exon (une ligne par patient), seuils DeCovA 0,8 / 1,2
  - `perbase_coverage_curve.R` : profondeur base par base (mosdepth), normalisée par échantillon et par cohorte, rapport PDF/PNG/HTML incluant la jonction TSC2-PKD1
  - `perbase_coverage_curve_memory_safe.R` : même courbe en pleine résolution pour les gros runs (PNG/PDF seulement)

> KODEX signale des variants **à revoir** ; il ne pose pas de diagnostic. Toute alerte se vérifie sur IGV et par une méthode orthogonale.

---

## Installation

### Prérequis

- [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) ou Conda
- Python ≥ 3.11

### Installation automatique

```bash
git clone https://github.com/kaarud/kodex.git
cd kodex
bash install.sh
```

Le script crée l'environnement `kodex`, installe les dépendances et lance les tests.

### Installation manuelle

```bash
micromamba env create -f environment.yml
micromamba activate kodex
```

---

## Usage

### Filtrage d'un fichier NiourK

```bash
bash run.sh sample_filtered.xlsx
```

Ou directement (schéma **run 2** par défaut) :

```bash
micromamba run -n kodex python vcf_filter.py sample.xlsx
```

### Choix du schéma (run 2 vs run 1)

Par défaut KODEX attend le schéma **run 2** (`config.yaml`). Pour traiter un ancien
fichier **run 1**, ajoutez `--legacy` (raccourci pour `--config config_legacy.yaml`) :

```bash
# Nouveau schéma run 2 (défaut)
micromamba run -n kodex python vcf_filter.py <RUN>_Singularity_<id>_S1.xlsx

# Ancien schéma run 1
micromamba run -n kodex python vcf_filter.py --legacy <RUN>_<id>_S7.xlsx
```

Un `--config` explicite l'emporte toujours sur `--legacy`. Si le fichier ne correspond
pas au profil choisi, KODEX s'arrête proprement avec un message indiquant la colonne
attendue manquante et suggère `--legacy` — sans *traceback*.

### Mode beta (sans pré-exclusions)

```bash
bash run.sh sample.xlsx --beta
```

### Agrégation de la cohorte

```bash
bash run.sh --aggregate /chemin/vers/xlsx/
```

### KODEX CNV — courbes de couverture

Dans l'environnement `kodex` (R ≥ 4.3, voir `environment.yml`) :

```bash
# Courbe RPKM par exon, depuis le classeur KoverGraph du run
Rscript rpkm_exon_curve.R --xlsx /chemin/vers/<RUN>_Blanket.xlsx --genes TSC1,TSC2

# Profondeur base par base (fichiers mosdepth *.per-base.bed.gz) — rapport complet par défaut
Rscript perbase_coverage_curve.R --depth-dir /chemin/vers/depth --genes TSC1,TSC2

# Une seule vue : --raw, --normalised, --cohort ; lissage --smooth exon|<bp>
Rscript perbase_coverage_curve.R --depth-dir /chemin/vers/depth --cohort --smooth exon

# Gros run, pleine résolution sans HTML
Rscript perbase_coverage_curve_memory_safe.R --depth-dir /chemin/vers/depth --out output/perbase_coverage
```

La géométrie du panel est lue dans `panel/panel.tsv` (non versionné : dérivé des
données NiourK) ; `--rebuild-panel true --kover-dir <dossier KoverGraph>` la régénère.
Les sorties vont dans `output/` sous un dossier horodaté.

**Rapport HTML** : `coverage_report.html` a besoin du dossier `lib/` posé à côté
de lui ; il ne s'ouvre pas directement depuis Google Drive. Voir
[`docs/NOTICE_CNV.md`](docs/NOTICE_CNV.md).

### Bruit de fond du run

```bash
micromamba run -n kodex python run_background.py sample_filtered.xlsx --bam-dir /chemin/vers/bam/
micromamba run -n kodex python bg_filter.py sample_filtered_bg.xlsx -o candidats.xlsx
```

### Captures IGV

```bash
micromamba run -n kodex python igv_capture.py sample_filtered.xlsx \
    --bam-dir /chemin/vers/bam/ \
    --out /chemin/vers/igv_captures/
```

### Bookmark IGV (Regions of Interest)

```bash
micromamba run -n kodex python igv_bookmark.py sample_filtered.xlsx
```

### Lien IGV par position (go-to)

Le navigateur se fait en deux temps, pour ne charger les BAMs qu'une seule fois
par patient :

1. **Une fois** — la feuille `_summary` contient un lien **« ▶ Charger les BAMs »**
   qui charge dans IGV Desktop le BAM brut du patient (RAW), plus `CHIM` et
   `MUTECT2` s'ils existent, plus le BAM brut d'un échantillon de référence
   (« normal », par défaut lane `S2`).
2. **Par variant** — chaque feuille de données reçoit une colonne **IGV** dont
   chaque ▶ ne fait que naviguer (`goto`) sur la position du variant, sans
   recharger les pistes.
3. **Tri par base** — l'en-tête **A1** de chaque feuille de données (« ⇅ base »)
   est un bouton qui trie l'alignement courant par base. IGV n'exécute qu'une
   seule commande par requête HTTP : on ne peut donc pas enchaîner `goto` + `sort`
   en un clic. Le geste est : ▶ pour sauter sur le variant, puis A1 pour trier.

Un softclip présent à la fois sur le patient et la référence est probablement un
artefact. Seul le BAM RAW du patient est obligatoire — les runs qui n'ont pas de
`.mutect2.bam` (p. ex. run 2) chargent simplement ce qui existe.

IGV doit être ouvert (port 60151). Les BAMs sont lus en local depuis `igv.bam_dir`.
Le premier clic « Charger les BAMs » ouvre un onglet navigateur (réponse `OK`
d'IGV) ; il peut être fermé. Chaque ▶ ouvre aussi un onglet, mais la navigation
est instantanée.

Équivalent en ligne de commande :

```bash
# par sample + position
micromamba run -n kodex python igv_goto.py 100008_S1 chr16:2098100

# depuis une ligne d'un fichier filtré (index 0-based)
micromamba run -n kodex python igv_goto.py sample_filtered.xlsx --sheet HIGH_IMPACT --row 7
```

Configuration (`config.yaml`, bloc `igv:`) : `port`, `genome`, `bam_dir`,
`reference_lane`, `half_window`.

---

## Configuration

Deux profils de schéma sont fournis :

- `config.yaml` — profil **run 2** (défaut)
- `config_legacy.yaml` — profil **run 1** (via `--legacy`)

| Section | Description |
|---|---|
| `column_aliases` | Renommage colonnes brutes → vocabulaire interne (069 : `GENE→SYMBOL`, `CSQ→Consequence`, `REGION→REGIONS`) |
| `af_filter_column` | Colonne gnomAD pilotant le filtrage de fréquence (069 : `GNOMAD4_AF` ; 063 : `MAX_AF_GNOMADEX2.1`) |
| `synthesize_user_columns` | Crée des colonnes `USER_*` vides si absentes (vrai en 069) |
| `exclusions` | Seuils gnomAD, RNCNT, PJCNT, SB, AF mosaïque |
| `scoring` | Poids des axes technique / biologique / rareté |
| `igv_exclusion_terms` | Termes excluant un variant des captures IGV |
| `af_discordance` | Seuils pour l'alerte HTZ↔mosaïque |
| `output_columns` | Colonnes exportées (jeu adapté à chaque schéma) |

> **Note (run 2)** : le filtrage de fréquence utilise `GNOMAD4_AF` (AF globale gnomAD v4),
> alors que les seuils historiques (`0.005` / `0.001`) avaient été calibrés sur la popmax
> `MAX_AF_GNOMADEX2.1`. Les seuils sont conservés tels quels ; une recalibration éventuelle
> est une tâche distincte.

---

## Tests

```bash
bash run.sh --test
# ou directement :
micromamba run -n kodex python -m pytest tests/ -v

# tests R (testthat), depuis tests/
cd tests && for f in *.R; do micromamba run -n kodex Rscript "$f"; done
```

---

## Structure du projet

```
vcf_filter.py            — filtrage principal, scoring, mise en forme XLSX
aggregate_annotations.py — agrégation multi-échantillons
igv_capture.py           — captures IGV automatisées
igv_bookmark.py          — injection des variants comme ROI dans IGV
igv_goto.py              — CLI : charge les BAMs (patient + référence) et navigue sur une position dans IGV
igv_link.py              — logique pure : résolution des BAMs, fenêtre locus, construction de l'URL IGV
run_background.py        — bruit de fond du run à la position de chaque variant
bg_filter.py             — sélection des variants à signature de témoin positif
minibam.py               — tranches BAM légères autour des variants annotés
sv_xlsx.py               — SV Manta / Delly → Excel avec liens IGV
coverage_common.R        — fonctions R partagées (seuils DeCovA, sorties horodatées, rendu)
rpkm_exon_curve.R        — KODEX CNV : courbe RPKM par exon
perbase_coverage_curve.R — KODEX CNV : couverture base par base, rapport PDF/PNG/HTML
perbase_coverage_curve_memory_safe.R — variante pleine résolution économe en mémoire
config.yaml              — configuration profil run 2 (défaut)
config_legacy.yaml       — configuration profil run 1 (--legacy)
environment.yml          — dépendances micromamba/conda
install.sh               — installation automatisée
run.sh                   — point d'entrée CLI
tests/                   — tests pytest (196) et testthat
docs/                    — documentation détaillée
AGENT.md                 — guide rapide du dépôt pour un agent
```

---

## Contexte

KODEX est développé dans le cadre de la réanalyse NGS de la cohorte **PRELUDE-TSC** (CHU Angers).
Il traite les sorties du pipeline **NiourK** (KAVE - CHU Angers, Nextflow), aussi bien le run
historique **run 1** (schéma `--legacy`) que le run récent **run 2** (schéma par défaut) —
panel Twist TSC, Illumina NextSeq, GRCh38.

---

## Citation

Si vous utilisez KODEX dans vos travaux, merci de le citer (voir `CITATION.cff` ou le bouton *Cite this repository* sur GitHub).

---

## Licence

KODEX est distribué sous licence [MIT](LICENSE).
Les données patients et les fichiers de résultats NiourK ne sont **jamais** inclus dans ce dépôt (voir `.gitignore`).

---

## Crédits

Une grande partie du code de KODEX a été générée avec l'assistance de l'IA
[Claude](https://www.anthropic.com/claude) (Anthropic). Chaque contribution a
été relue et validée par l'auteur avant intégration.
