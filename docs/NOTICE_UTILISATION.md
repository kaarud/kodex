---
title: "Notice d'utilisation — Pipeline de filtrage NGS PRELUDE-TSC"
subtitle: "Installation, configuration et exploitation sur serveur Ubuntu"
author: "PRELUDE-TSC / CHU Angers"
date: "Mars 2026"
lang: fr
---

# 1. Pré-requis

| Composant | Version minimum | Notes |
|-----------|----------------|-------|
| Ubuntu | 20.04 LTS | Testé sur 22.04 et 24.04 |
| micromamba | 1.5+ | Installé automatiquement par `install.sh` |
| Espace disque | ~500 Mo | Environnement conda + données |
| Accès réseau | Au premier lancement | Téléchargement des paquets conda-forge |

Le pipeline ne nécessite **aucun droit root** — tout s'installe dans l'espace utilisateur.

# 2. Installation

## 2.1. Récupération du code

```bash
# Depuis le dépôt Git (recommandé)
git clone <url_du_depot> 07_ngs_reanalysis
cd 07_ngs_reanalysis

# Ou par copie directe
scp -r user@machine:~/projets/PRELUDE-TSC/07_ngs_reanalysis .
cd 07_ngs_reanalysis
```

## 2.2. Installation automatique

```bash
bash install.sh
```

Le script :

1. Installe **micromamba** si absent (dans `~/micromamba` par défaut)
2. Crée l'environnement `ngs_reanalysis` (Python 3.11, pandas, openpyxl, PyYAML, requests, pytest)
3. Vérifie que toutes les dépendances sont opérationnelles
4. Lance la suite de tests (56 tests)

Pour installer micromamba dans un chemin spécifique :

```bash
bash install.sh /opt/micromamba
```

## 2.3. Vérification manuelle

```bash
micromamba activate ngs_reanalysis
python -c "import pandas, openpyxl, yaml, requests; print('OK')"
python -m pytest tests/ -q
```

## 2.4. Structure des fichiers

```
07_ngs_reanalysis/
├── vcf_filter.py              # Script principal de filtrage
├── aggregate_annotations.py   # Agrégation multi-échantillons
├── igv_capture.py             # Captures IGV automatisées
├── igv_bookmark.py            # Bookmarks IGV (régions d'intérêt)
├── config.yaml                # Configuration centralisée
├── environment.yml            # Définition de l'environnement conda
├── install.sh                 # Script d'installation
├── run.sh                     # Script de lancement simplifié
├── tests/                     # Suite de tests (56 tests)
│   ├── conftest.py
│   ├── test_cli.py
│   ├── test_exclusions.py
│   ├── test_load.py
│   ├── test_tiers.py
│   └── test_igv_capture.py
└── docs/                      # Documentation
```

---

# 3. Utilisation

## 3.1. Filtrage d'un échantillon (commande principale)

```bash
# Avec le wrapper
bash run.sh /chemin/vers/fichier_niourk.xlsx

# Ou directement
micromamba activate ngs_reanalysis
python vcf_filter.py /chemin/vers/fichier_niourk.xlsx
```

**Entrée** : fichier XLSX NiourK contenant la feuille `Variants` (sortie du pipeline NiourK v2.5).

**Sortie** : fichier `*_filtered.xlsx` dans le même répertoire, contenant :

- `_summary` : rapport de pré-filtrage, hash de configuration, date
- Jusqu'à **12 onglets** d'analyse (HIGH_IMPACT, CLINVAR_P_LP, MOSAIC, etc.)
- Onglet **RHEB_PKD1** (gènes voisins)
- Coloration conditionnelle sur toutes les colonnes clés
- Score de qualité (gradient rouge→vert) sur la cellule Import DefGen
- Formules dynamiques USER_COM pour les doublons inter-onglets

## 3.2. Mode beta (sans pré-exclusions)

```bash
bash run.sh /chemin/vers/fichier.xlsx --beta
```

Désactive les 6 filtres d'exclusion — tous les variants TSC1/TSC2 sont classés directement dans les onglets. Utile pour vérifier qu'aucun variant n'est exclu à tort.

## 3.3. Agrégation multi-échantillons

Après avoir filtré et annoté individuellement chaque échantillon :

```bash
# Agréger tous les *_filtered.xlsx d'un répertoire
bash run.sh --aggregate /chemin/vers/dossier_xlsx/

# Ou des fichiers spécifiques
micromamba activate ngs_reanalysis
python aggregate_annotations.py fichier1_filtered.xlsx fichier2_filtered.xlsx
```

**Sortie** : `TSC_063_annotations.xlsx` — tableau de synthèse de tous les variants annotés (USER_CLASS/ANNOT/COM non vides) avec la colonne SAMPLE et SOURCE_SHEET.

## 3.4. Captures IGV (poste de travail uniquement)

> **Pré-requis** : IGV Desktop ouvert avec le port TCP activé (60151 par défaut).

```bash
# Captures automatiques (3 vues × N variants)
python igv_capture.py fichier_filtered.xlsx --bam-dir /chemin/vers/bams/

# Bookmarks IGV (régions d'intérêt)
python igv_bookmark.py fichier_filtered.xlsx
```

Les scripts IGV ne sont pertinents que sur un poste avec écran — pas sur le serveur.

## 3.5. Tests

```bash
bash run.sh --test
# Ou : python -m pytest tests/ -v
```

---

# 4. Configuration (`config.yaml`)

Toute modification de seuils se fait dans ce fichier — **aucun code Python à modifier**.

## 4.1. Gènes cibles

```yaml
genes: [TSC1, TSC2]
```

## 4.2. Seuils d'exclusion

```yaml
exclusions:
  gnomad_max_af: 0.005       # > 0.5% gnomAD → exclu
  rncnt_max_total: 6         # > 6 récurrences dans le run → exclu
  pjcnt_max_total: 10        # > 10 dans la méta-analyse → exclu
  exclude_hom: true           # homozygotes récurrents → exclus
  exclude_clinvar_benign: true
  mosaic_artifact:
    enabled: true
    rncnt_low_max: 6
    af_max_pct: 35.0
```

## 4.3. Onglets d'analyse

Chaque onglet peut être activé/désactivé et ses seuils ajustés :

```yaml
sheets:
  high_impact:
    enabled: true
  clinvar_p_lp:
    enabled: true
  mosaic:
    enabled: true
    af_min_pct: 5.0
    af_max_pct: 35.0
    rncnt_max_total: 2
  deep_mosaic:
    enabled: true
    af_min_pct: 0.5
    af_max_pct: 5.0
    min_dp: 200
  extra_genes: [RHEB, PKD1]
```

## 4.4. Score de qualité — système RGB 3 axes

La couleur de la cellule Import DefGen encode **trois dimensions indépendantes** :

| Canal | Dimension | Critères |
|-------|-----------|----------|
| **R** (rouge) | Risque artefact = 1 − qualité technique | SB (strand bias), CALLNB, DP/AF |
| **G** (vert)  | Signal biologique | ClinVar, IMPACT, CADD, REVEL, SpliceAI |
| **B** (bleu)  | Rareté du variant | RNCNT, PJCNT, gnomAD, LOVD |

Poids par défaut (ajustables dans `config.yaml`) :

```yaml
scoring:
  technical:                     # → canal R (inverted)
    sb:       { weight: 2.0, good_min: 0.3, good_max: 0.7 }
    callnb:   { weight: 2.0, excellent: 4 }
    dp_af:    { weight: 0.5, min_alt_reads: 5 }
  biological:                    # → canal G
    clinvar:  { weight: 2.5, pathogenic: 1.0, absent: 0.4 }
    impact:   { weight: 2.0, HIGH: 1.0, MODERATE: 0.75 }
    cadd:     { weight: 1.0, good: 30 }
    revel:    { weight: 1.0, good: 0.75 }
    spliceai: { weight: 1.0, good: 0.8 }
  rarity:                        # → canal B
    rncnt:    { weight: 1.5, good_max: 1 }   # ≤ 1 = patient unique dans la base
    pjcnt:    { weight: 1.0, good_max: 0 }
    gnomad:   { weight: 1.5, ultra_rare: 0.0001, rare: 0.001 }
    lovd:     { weight: 0.5 }
    gnomad_pjcnt_benign_threshold: 10  # absent gnomAD + PJCNT ≥ 10 → artefact panel
```

## 4.5. Signalement des discordances AF/run

```yaml
af_discordance:
  htz_recurrent_in_mosaic:     # HTZ récurrent → ici mosaïque
    enabled: true
    min_rncnt_htz: 3
    mosaic_af_max: 35.0
  mosaic_recurrent_in_htz:     # Mosaïque récurrent → ici HTZ
    enabled: true
    min_rncnt_low: 3
    htz_af_min: 35.0
    htz_af_max: 65.0
```

## 4.6. Termes d'exclusion IGV

```yaml
igv_exclusion_terms:
  - "SB"
  - "ReMM low"
  - "Fréquence"
```

---

# 5. Interprétation du fichier de sortie

## 5.1. Code couleur des cellules

Les colonnes (SB, CALLNB, IMPACT, CLINVAR…) utilisent un code couleur classique :

| Couleur | Signification |
|---------|---------------|
| Vert | Bon / favorable (SB équilibré, gnomAD absent, 4 callers) |
| Jaune | Attention / modéré |
| Orange | Alerte (1 seul caller, deep mosaic) |
| Rouge | Défavorable (SB biaisé fort, ClinVar Benign) |

## 5.2. Import DefGen — Score RGB 3 axes

La cellule Import DefGen utilise un **code couleur RGB** encodant trois dimensions simultanément :

| Canal | Dimension | Élevé = |
|-------|-----------|---------|
| **R** rouge | Risque artefact (1 − qualité technique) | Séquençage douteux |
| **G** vert  | Signal biologique | Variant pertinent |
| **B** bleu  | Rareté | Variant novel/unique |

### Signatures visuelles

| Couleur | Hex | R | G | B | Interprétation | Exemple |
|---------|-----|---|---|---|----------------|---------|
|  **Cyan** | `#00FFFF` | 0 | max | max | **Candidat idéal** | ClinVar P, SB=0.5, 4 callers, gnomAD absent, RNCNT=1 |
|  **Vert pur** | `#00FF00` | 0 | max | 0 | Bio forte mais **non rare** — récurrent dans nos runs ou artefact panel | ClinVar P, PJCNT=12, gnomAD absent (artefact panel) |
|  **Jaune** | `#FFFF00` | max | max | 0 | **Faux-positif dangereux** — biologie forte + artefact technique | ClinVar P, SB=0.97, 1 caller, mais gnomAD fréquent |
|  **Rouge** | `#FF0000` | max | 0 | 0 | **Artefact pur** | SB extrême, 1 caller, MODIFIER, gnomAD > 0.1% |
|  **R_O_se** | `#FD6C9E` | max | 0 | 0 | **Artefact pur** | SB extrême, 1 caller, MODIFIER, gnomAD > 0.1% |
|  **Noir** | `#000000` | 0 | 0 | 0 | Sans intérêt biologique, séquençage propre | MODIFIER, commun, technique parfaite |
|  **Bleu-vert** | `#1C88A0` | faible | moyen | moyen | VUS/splice modéré, rare | SpliceAI=0.72, gnomAD ultra-rare, RNCNT=1 |


### Règle spéciale — Paradoxe gnomAD

> **Absent gnomAD + PJCNT ≥ 10** → contribution gnomAD annulée.
> L'absence dans la population générale n'est pas un argument de rareté si le variant est très fréquent dans nos runs — c'est un **artefact panel** (variant spécifique à notre technique, invisible en gnomAD car jamais séquencé avec notre panel). La cellule passe en vert pur (B=0) plutôt qu'en cyan.

Au survol de la cellule Import DefGen : commentaire détaillé avec scores Tech/Bio/Rareté et diagnostic ligne par ligne de chaque critère.

## 5.3. Police bleue sur HGVSc

Un variant en **police bleue** apparaît dans **plusieurs onglets** du classeur. Le commentaire au survol liste les onglets concernés.

## 5.4. Fond rouge sur HGVSc

Un variant avec **fond rouge sur HGVSc** présente une **discordance AF/run** :

- HTZ récurrent dans le run mais ici en mosaïque → artefact probable
- Mosaïque récurrent mais ici en HTZ → **variant potentiellement pathogène, à vérifier sur IGV**

## 5.5. Colonnes USER_CLASS / USER_ANNOT / USER_COM

| Colonne | Usage |
|---------|-------|
| USER_CLASS | Classe ACMG (1–5) ou terme d'exclusion (SB, ReMM low) |
| USER_ANNOT | Annotation libre du biologiste |
| USER_COM | Formule automatique : agrège les annotations des autres onglets |

---

# 6. Flux de travail recommandé

```
1. Récupérer le XLSX NiourK depuis le serveur CHU
   scp sl157:/data/results/TSC_063/sample.xlsx .

2. Filtrer
   bash run.sh sample.xlsx

3. Ouvrir dans Excel → interpréter les variants onglet par onglet
   → Annoter USER_CLASS (classe ACMG) et USER_ANNOT (commentaire)
   → Les formules USER_COM propagent automatiquement

4. (Optionnel) Vérifier sur IGV les variants signalés
   python igv_bookmark.py sample_filtered.xlsx
   python igv_capture.py sample_filtered.xlsx --bam-dir /chemin/bams/

5. Agréger après avoir traité tous les échantillons
   bash run.sh --aggregate /dossier/xlsx/
```

---

# 7. Dépannage

| Problème | Solution |
|----------|----------|
| `micromamba: command not found` | Relancer `bash install.sh` ou sourcer le profil : `eval "$(micromamba shell hook -s bash)"` |
| `ModuleNotFoundError: No module named 'pandas'` | Activer l'environnement : `micromamba activate ngs_reanalysis` |
| `KeyError: 'Variants'` | Le fichier XLSX n'est pas au format NiourK (feuille `Variants` manquante) |
| Tests en échec | Vérifier la version Python (3.11+) et relancer `micromamba update -n ngs_reanalysis -f environment.yml` |
| Formule `#NOM?` dans Excel | Normal si ouvert avec LibreOffice — nécessite Microsoft Excel 2010+ |
| IGV non joignable | Vérifier qu'IGV Desktop est ouvert et que `Enable remote control` est coché dans View > Preferences > Advanced |

---

# 8. Sécurité et isolation

- L'environnement `ngs_reanalysis` est **entièrement isolé** du système : aucune dépendance système, aucun paquet Python global utilisé
- Les scripts ne font **aucun accès réseau** (sauf `igv_capture.py` qui communique en localhost avec IGV Desktop)
- Les données ne sont **jamais modifiées en place** — le fichier d'entrée n'est pas touché, un nouveau `*_filtered.xlsx` est créé
- Le fichier `config.yaml` est versionné et son hash MD5 est enregistré dans chaque fichier de sortie (`_summary`) pour traçabilité
- Aucune donnée patient n'est stockée dans le code ou la configuration
