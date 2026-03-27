# KODEX

**KODEX** est un plugin de post-traitement des variants générés par la pipeline angevine NiourK pour la réanalyse NGS TSC1/TSC2.
Il déchiffre les fichiers XLSX produits par le pipeline NiourK (CHU Angers) et génère des tableaux filtrés, annotés et prêts pour la revue clinique. En somme il permet de passer d'un fichier contenants des millers de variants à plusieurs feuilles dans un classeur excel contenant une dizaine de variants.

---

## Fonctionnalités

- **Filtrage multi-critères** — gnomAD, RNCNT, PJCNT, IMPACT, CLINVAR, LOVD, SpliceAI
- **Score de qualité par variant** — 3 axes : technique, biologique, rareté → couleur dégradée sur la cellule Import DefGen
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

---

## Installation

### Prérequis

- [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) ou Conda
- Python ≥ 3.11

### Installation automatique

```bash
git clone https://github.com/<user>/kodex.git
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

Ou directement :

```bash
micromamba run -n kodex python vcf_filter.py sample.xlsx --config config.yaml
```

### Mode beta (sans pré-exclusions)

```bash
bash run.sh sample.xlsx --beta
```

### Agrégation de la cohorte

```bash
bash run.sh --aggregate ~/data/TSC_063/xlsx/
```

### Captures IGV

```bash
micromamba run -n kodex python igv_capture.py sample_filtered.xlsx \
    --bam-dir ~/data/TSC_063/bam/ \
    --out ~/data/TSC_063/igv_captures/
```

### Bookmark IGV (Regions of Interest)

```bash
micromamba run -n kodex python igv_bookmark.py sample_filtered.xlsx
```

---

## Configuration

Tous les paramètres sont dans `config.yaml` :

| Section | Description |
|---|---|
| `filters` | Seuils gnomAD, RNCNT, PJCNT, SB, AF mosaïque |
| `scoring` | Poids des axes technique / biologique / rareté |
| `igv_exclusion_terms` | Termes excluant un variant des captures IGV |
| `af_discordance` | Seuils pour l'alerte HTZ↔mosaïque |

---

## Tests

```bash
bash run.sh --test
# ou directement :
micromamba run -n kodex python -m pytest tests/ -v
```

---

## Structure du projet

```
vcf_filter.py            — filtrage principal, scoring, mise en forme XLSX
aggregate_annotations.py — agrégation multi-échantillons
igv_capture.py           — captures IGV automatisées
igv_bookmark.py          — injection des variants comme ROI dans IGV
config.yaml              — configuration centralisée
environment.yml          — dépendances micromamba/conda
install.sh               — installation automatisée
run.sh                   — point d'entrée CLI
tests/                   — suite de tests pytest (56 tests)
docs/                    — documentation détaillée
```

---

## Contexte



---

## Citation

Si vous utilisez KODEX dans vos travaux, merci de le citer (voir `CITATION.cff` ou le bouton *Cite this repository* sur GitHub).

---

## Licence

KODEX est distribué sous licence [GNU GPL v3](LICENSE).
Les données patients et les fichiers de résultats NiourK ne sont **jamais** inclus dans ce dépôt (voir `.gitignore`).
