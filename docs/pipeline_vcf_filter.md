---
title: "Pipeline de filtrage et analyse des variants NGS — PRELUDE-TSC"
subtitle: "Documentation technique — vcf_filter.py et outils associés"
author: "PRELUDE-TSC / CHU Angers"
date: "Mars 2026"
lang: fr
---

# 1. Vue d'ensemble

Le pipeline `07_ngs_reanalysis` est un ensemble d'outils Python conçus pour la **réanalyse systématique des variants NGS court-read** (Illumina) issus du pipeline NiourK v2.5 (KAVE, CHU Angers). Il prend en entrée les fichiers XLSX annotés par NiourK et produit un classeur Excel multi-onglets structuré, prêt pour l'interprétation biologique.

**Contexte clinique** : cohorte PRELUDE-TSC — patients atteints de sclérose tubéreuse de Bourneville (STB), gènes cibles TSC1 et TSC2.

**Run courant** : TSC_063 — 22 échantillons, Illumina NextSeq 101 bp PE, panel Twist TSC CHU, alignement GRCh38, annotation VEP v106 + SpliceAI + CADD + REVEL + SPiP.

## Architecture des scripts

| Script | Rôle |
|---|---|
| `vcf_filter.py` | Filtrage, classification en onglets, scoring, formatage Excel |
| `igv_capture.py` | Capture automatisée de screenshots IGV (3 vues par variant) |
| `igv_bookmark.py` | Injection de régions d'intérêt (ROI) dans IGV Desktop |
| `aggregate_annotations.py` | Agrégation multi-échantillons des variants annotés |
| `config.yaml` | Configuration centralisée (seuils, poids, onglets) |

## Flux de travail

```
NiourK XLSX (par échantillon)
    │
    ▼
vcf_filter.py ──→ *_filtered.xlsx (classeur multi-onglets)
    │                    │
    │                    ├──→ igv_bookmark.py ──→ ROI dans IGV Desktop
    │                    │
    │                    └──→ igv_capture.py ──→ Screenshots PNG (3 vues/variant)
    │
    ▼
aggregate_annotations.py ──→ TSC_063_annotations.xlsx (synthèse multi-échantillons)
```

---

# 2. vcf_filter.py — Filtrage et classification

## 2.1. Chargement et normalisation

- Lecture de la feuille `Variants` du XLSX NiourK
- Suppression automatique du préfixe échantillon NiourK (`<sample_id>\n<col_name>` → `<col_name>`)
- Filtrage par gènes cibles (`genes: [TSC1, TSC2]` dans le YAML)
- Dédoublonnage par transcrit de référence (`REF_TRANS = VRAI`) : pour chaque position génomique, seul le transcrit canonique est conservé

## 2.2. Pré-exclusions séquentielles

Sept filtres d'exclusion sont appliqués **séquentiellement** (chaque filtre opère sur les survivants du précédent). Tous les seuils sont configurables dans `config.yaml` :

| # | Filtre | Seuil par défaut | Logique |
|---|--------|------------------|---------|
| 1 | **gnomAD** | > 0.5 % (0.005) | Variant trop fréquent en population générale |
| 2 | **RNCNT_TOTAL** | > 4 | Récurrent dans le run courant (22 échantillons) |
| 3 | **PJCNT_TOTAL** | > 10 | Récurrent dans la méta-analyse de 30 runs |
| 4 | **RNCNT_HOM** | > 0 | Homozygote récurrent → probable artefact |
| 5 | **ClinVar Benign** | Benign/Likely_benign | Classification ClinVar bénigne |
| 6 | **Artefact HTZ→mosaïque** | Ratio HTZ > 80 %, ≥ 3 HTZ, AF < 35 % | Faux positif mosaïque (hétérozygotes récurrents à basse fréquence) |
| 7 | **Artefact mosaïque récurrent** | RNCNT_LOW > 5, AF < 35 % | Hotspot de mosaïque artefactuel |

Un **pré-rapport** est généré sur la sortie standard et inclus dans la feuille `_summary` du classeur, détaillant le nombre de variants exclus à chaque étape.

## 2.3. Classification en onglets (fiche d'analyse STB)

Les variants survivants sont classés dans **12 onglets indépendants** (un variant peut apparaître dans plusieurs onglets) :

| Onglet | Critère | Objectif clinique |
|--------|---------|-------------------|
| **HIGH_IMPACT** | IMPACT = HIGH | Stop-gained, frameshift, splice acceptor/donor |
| **CLINVAR_P_LP** | ClinVar Pathogenic/Likely_pathogenic (hors Conflicting) | Variants déjà classés pathogènes |
| **CLINVAR_CONFLICTING** | ClinVar Conflicting | Interprétations discordantes entre labos |
| **CLINVAR_VUS** | ClinVar Uncertain_significance | Variants de signification incertaine |
| **CLINVAR_NOT_PROVIDED** | ClinVar not_provided | Soumis mais non classés |
| **FQ_MODERATE** | gnomAD < 0.1 % + IMPACT MODERATE | Missense rares |
| **FQ_LOW** | gnomAD < 0.1 % + IMPACT LOW | Synonymes rares (potentiel splice) |
| **FQ_MODIFIER** | gnomAD < 0.1 % + MODIFIER + (SpliceAI ≥ 0.5 OU SPiP Altered) | Introniques avec signal d'épissage |
| **4_CALLERS** | CALLNB ≥ 4 | Consensus des 4 variant callers (Mutect2, Strelka2, GATK HC, DeepVariant) |
| **LOVD** | Entrée LOVD non vide | Variant répertorié dans la base LOVD |
| **MOSAIC** | 5 % < AF ≤ 35 %, RNCNT ≤ 2 | Mosaïcisme somatique standard |
| **DEEP_MOSAIC** | 0.5 % < AF ≤ 5 %, RNCNT ≤ 2, DP ≥ 200, SB non extrême | Mosaïcisme profond (sous le seuil classique) |

Un onglet supplémentaire **RHEB_PKD1** regroupe tous les variants des gènes RHEB et PKD1 (gènes voisins de TSC1/TSC2 sur les mêmes loci) sans filtre de qualité.

---

# 3. Formatage Excel avancé

## 3.1. Coloration conditionnelle des cellules

Chaque cellule de données est colorée selon des règles spécifiques à sa colonne, avec une palette cohérente (vert → jaune → orange → rouge) :

| Colonne | Vert | Jaune | Orange | Rouge |
|---------|------|-------|--------|-------|
| **AF_PCT** | 40–60 % (hétérozygote) | 5–40 % (mosaïque) | 0.5–5 % (deep mosaic) | — |
| **SB** | 0.4–0.6 | 0.3–0.4 / 0.6–0.7 | — | < 0.3 / > 0.7 |
| **CALLNB** | 4 callers | 2 callers | 1 caller | 0 |
| **RNCNT_TOTAL** | 0 (unique) | 3–4 | — | — |
| **IMPACT** | — | MODERATE | LOW | HIGH |
| **REVEL** | < 0.5 (bénin) | 0.5–0.75 | — | ≥ 0.75 (pathogène) |
| **CADD_PHRED** | — | 20–30 | — | ≥ 30 |
| **SpliceAI** | — | 0.5–0.8 | — | ≥ 0.8 |
| **USER_CLASS** | 1–2 (bénin) | 3 (VUS) | 4 (probablement patho) | 5 (pathogène) |

## 3.2. Score de qualité variant (Import DefGen)

Chaque variant reçoit un **score de qualité global** (0–100 %) affiché sous forme de gradient de couleur sur la cellule `Import DefGen` :

- **Rouge** (0 %) → **Jaune** (50 %) → **Vert** (100 %)

Le score est la moyenne pondérée de deux axes indépendants :

### Axe technique (fiabilité du calling)

| Critère | Poids | Évaluation |
|---------|-------|------------|
| Strand Bias (SB) | 2.0 | 0.3–0.7 = bon, < 0.2 ou > 0.8 = mauvais |
| Nombre de callers (CALLNB) | 2.0 | 4 = excellent, 3 = bon, 2 = modéré, 1 = suspect |
| Récurrence run (RNCNT) | 1.5 | 0 = unique, 1–2 = modéré, > 4 = suspect |
| Récurrence inter-projets (PJCNT) | 1.0 | 0 = unique, > 5 = suspect |
| Cohérence DP/AF | 0.5 | ≥ 5 reads alt attendus = bon |

### Axe biologique (pertinence clinique)

| Critère | Poids | Évaluation |
|---------|-------|------------|
| ClinVar | 2.5 | Pathogenic = 100 %, VUS = 60 %, Absent = 40 %, Benign = 0 % |
| Impact fonctionnel (VEP) | 2.0 | HIGH = 100 %, MODERATE = 75 %, LOW = 25 %, MODIFIER = 0 % |
| gnomAD fréquence | 1.5 | Absent = 100 %, ultra-rare < 0.01 % = 67 %, rare < 0.1 % = 33 % |
| CADD | 1.0 | ≥ 30 = 100 %, 20–30 = 50 % |
| REVEL | 1.0 | ≥ 0.75 = 100 %, 0.5–0.75 = 50 % |
| SpliceAI | 1.0 | ≥ 0.8 = 100 %, 0.5–0.8 = 50 % |
| LOVD | 0.5 | Présent = 100 % |

### Commentaire détaillé au survol

Au survol de la cellule `Import DefGen`, un commentaire structuré s'affiche :

```
===================================
  SCORE : 78% — Variant solide, à valider
  Tech 4/5  |  Bio 3/4
===================================

(+) SB équilibré, 3 callers, ClinVar P/LP, HIGH impact (stop_gained)
(-) RNCNT=2 récurrent

───────────────────────────────────
TECHNIQUE (4/5) :
  ✓ SB=0.45 (équilibré)
  ✓ CALLNB=3 (3 callers concordants)
  ~ RNCNT=2 (peu fréquent)
  ✓ PJCNT=0 (jamais vu inter-projets)
  ✓ DP=350, AF=47.2% (~165 reads alt)

───────────────────────────────────
BIOLOGIQUE (3/4) :
  ✓ ClinVar: Pathogenic
  ✓ IMPACT: HIGH (stop_gained)
  ✗ SpliceAI=0.12
  ✓ Absent gnomAD (novel)
```

Les verdicts sont configurables :

| Score | Verdict par défaut |
|-------|--------------------|
| ≥ 70 % | Variant solide, à valider |
| ≥ 50 % | À évaluer — signaux mixtes |
| ≥ 30 % | Probablement artefact ou inintéressant |
| < 30 % | Artefact probable |

**Tous les poids et seuils du scoring sont modifiables dans `config.yaml`** sans toucher au code.

## 3.3. Détection de doublons inter-onglets

Un variant peut apparaître dans plusieurs onglets (ex. : un variant HIGH + ClinVar Pathogenic + LOVD apparaît dans 3 feuilles). Le pipeline détecte ces doublons et applique :

- **Police bleue** sur la cellule `HGVSc` pour signaler visuellement le doublon
- **Commentaire au survol** listant les autres onglets : *« Aussi dans : CLINVAR_P_LP, LOVD »*

## 3.4. Formules dynamiques USER_COM

Pour chaque variant présent dans plusieurs onglets, une **formule Excel dynamique** est injectée dans la colonne `USER_COM`. Cette formule agrège automatiquement les annotations `USER_CLASS` et `USER_ANNOT` saisies dans les autres onglets :

```
=MID(IF(OR('HIGH_IMPACT'!J2<>"", 'HIGH_IMPACT'!K2<>""),
  " | HIGH_IMPACT: "&'HIGH_IMPACT'!J2&" "&'HIGH_IMPACT'!K2, "")
&IF(OR('LOVD'!J2<>"", 'LOVD'!K2<>""),
  " | LOVD: "&'LOVD'!J2&" "&'LOVD'!K2, ""), 4, 9999)
```

**Avantage** : lorsque le biologiste annote un variant dans un onglet (ex. classe ACMG « 5 » dans `USER_CLASS`), cette annotation se propage automatiquement dans tous les autres onglets contenant le même variant, évitant les saisies redondantes.

Compatible Excel 2010+ (pas de TEXTJOIN).

## 3.5. Détection des variants downstream dans les gènes voisins

Les gènes TSC2 et PKD1 sont adjacents sur le chromosome 16. Un variant annoté `downstream_gene_variant` pour TSC2 peut en réalité se situer dans PKD1 (ou dans la région intergénique entre les deux).

Le pipeline identifie ces cas par coordonnées génomiques (GRCh38) et applique :

- **Fond bleu + police bleue** sur la cellule `Consequence`
- **Police bleue** sur `HGVSc` dans la feuille RHEB_PKD1
- **Commentaire** : *« Position dans PKD1 »*

---

# 4. igv_capture.py — Captures IGV automatisées

## 4.1. Principe

Pour chaque variant annoté par le biologiste (cellules `USER_CLASS`, `USER_ANNOT` ou `USER_COM` non vides), le script pilote IGV Desktop via son API TCP socket et génère **3 captures d'écran** :

| Vue | Paramètre IGV | Objectif |
|-----|---------------|----------|
| **Strand** | `COLOR_BY READ_STRAND` | Détecter un biais de brin (artefact) |
| **Soft-clip** | `SHOW_SOFT_CLIPPED true` | Repérer les réarrangements/insertions |
| **Squished** | `DISPLAY_MODE SQUISHED` | Vue d'ensemble de la couverture |

## 4.2. Fonctionnement

1. Extraction de l'ID échantillon depuis le nom de fichier (`TSC_063_345822_S2_filtered.xlsx` → `345822_S2`)
2. Localisation automatique des 3 BAMs (RAW, MUTECT2, CHIM) et de leurs index
3. Initialisation IGV : reset session, chargement du génome hg38, chargement des 3 BAMs
4. Pour chaque variant :
   - Calcul de la hauteur de panel dynamique (`DP × 2 + 150` pixels)
   - Navigation au locus (±145 pb autour de la position)
   - Génération des 3 captures dans `~/data/TSC_063/igv_captures/<sample_id>/`

## 4.3. Termes d'exclusion IGV

Les variants dont les annotations contiennent certains termes sont **exclus** des captures IGV (configurables dans `config.yaml`) :

- `SB` — biais de brin avéré
- `ReMM low` — score ReMM faible
- `Fréquence` / `freq` / `fréq` — variant trop fréquent

La correspondance est **insensible à la casse** et fonctionne en sous-chaîne.

---

# 5. igv_bookmark.py — Régions d'intérêt IGV

Injecte les variants annotés comme **régions d'intérêt (ROI)** dans IGV Desktop, permettant une navigation rapide via le `Region Navigator`.

Pour chaque variant annoté :

- Commande `region` envoyée à IGV (description incluant le nom du variant et les annotations utilisateur)
- Export d'un **fichier BED** (`*.roi.bed`) avec descriptions complètes (compatible `Regions > Import Regions` d'IGV)

Les annotations personnelles du biologiste (`USER_CLASS`, `USER_ANNOT`, `USER_COM`) apparaissent dans la colonne Description du Region Navigator.

---

# 6. aggregate_annotations.py — Synthèse multi-échantillons

## 6.1. Objectif

Après l'interprétation individuelle de chaque échantillon, ce script agrège les variants annotés de **tous les `*_filtered.xlsx`** en un tableau de synthèse unique.

## 6.2. Fonctionnalités

- Parcours de toutes les feuilles (sauf `_summary`) de chaque classeur filtré
- Sélection des lignes ayant au moins un champ `USER_CLASS`, `USER_ANNOT` ou `USER_COM` non vide
- Ignore les cellules contenant des formules (`=...`) — seules les annotations manuelles sont retenues
- Dédoublonnage par position génomique (`CHR + POS + REF + ALT`)
- Ajout de la colonne `SAMPLE` (ID échantillon) et `SOURCE_SHEET` (onglets d'origine)
- Écriture d'une feuille `_annotations` dans chaque `*_filtered.xlsx` (résumé local)
- Écriture d'un fichier de synthèse global `TSC_063_annotations.xlsx`

## 6.3. Colonnes de sortie prioritaires

Le tableau de synthèse ordonne les colonnes par importance clinique : identification du variant (SAMPLE, SYMBOL, HGVSc/p), annotations utilisateur, métriques techniques (AF, DP, SB, CALLNB), bases de données (ClinVar, LOVD), scores prédictifs (CADD, REVEL, SpliceAI, SPiP), fréquences populationnelles (gnomAD).

---

# 7. Configuration centralisée (config.yaml)

L'ensemble du pipeline est piloté par un **fichier YAML unique**. Toute modification de seuils, poids ou comportement se fait sans toucher au code Python.

## Sections du fichier de configuration

| Section | Contenu |
|---------|---------|
| `genes` | Gènes cibles du filtrage principal (`[TSC1, TSC2]`) |
| `exclusions` | Seuils des 7 pré-exclusions |
| `sheets` | Activation et paramétrage de chaque onglet |
| `scoring` | Poids et seuils du score de qualité (axes technique + biologique) |
| `igv_exclusion_terms` | Termes excluant un variant des captures/bookmarks IGV |
| `output_columns` | Ordre et sélection des colonnes dans le classeur de sortie |

---

# 8. Suite de tests

Le pipeline est couvert par **55 tests automatisés** (pytest) :

| Fichier | Tests | Couverture |
|---------|-------|------------|
| `test_load.py` | Chargement, préfixe NiourK, dédoublonnage REF_TRANS |
| `test_exclusions.py` | Chaque filtre d'exclusion individuellement |
| `test_tiers.py` | Classification en onglets |
| `test_cli.py` | Intégration end-to-end, scoring, gradient couleur, doublons, formules |
| `test_igv_capture.py` | Chargement variants annotés, exclusion terms, captures IGV, bookmarks |

Les tests utilisent des **fixtures synthétiques** (XLSX générés en mémoire) et des mocks pour IGV, garantissant une exécution rapide et reproductible sans données réelles.

---

# 9. Résumé des fonctionnalités clés

1. **Filtrage rigoureux** — 7 niveaux d'exclusion séquentielle, tous configurables
2. **12 onglets d'analyse** — classification multi-critères couvrant impact, ClinVar, fréquence, mosaïcisme, épissage, consensus callers
3. **Score de qualité bi-axial** — gradient visuel rouge→vert avec commentaire détaillé au survol
4. **Détection de doublons inter-onglets** — signalisation visuelle et propagation automatique des annotations
5. **Formules dynamiques** — agrégation Excel temps réel des annotations entre onglets
6. **Captures IGV automatisées** — 3 vues par variant (strand, soft-clip, squished), pilotage TCP
7. **Bookmarks IGV** — navigation rapide via Region Navigator avec annotations personnelles
8. **Détection gènes voisins** — signalisation des `downstream_gene_variant` potentiellement dans PKD1/RHEB
9. **Agrégation multi-échantillons** — synthèse globale des variants annotés de la cohorte
10. **Configuration YAML unique** — zéro modification de code pour adapter les seuils
