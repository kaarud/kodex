# CLAUDE.md — Agent Réanalyse NGS Court-Read

This file provides guidance to Claude Code when working in `07_ngs_reanalysis/`.

## Mission

Gérer la **réanalyse systématique des données NGS court-read** (Illumina) des cas NMI de la cohorte PRELUDE-TSC. Lire `BRIEFING_NGS_REANALYSIS.md` en premier — il contient le contexte complet de la réanalyse, les fichiers de données, et le pipeline bioinformatique.

Ce sous-agent couvre :
- L'harmonisation de la cohorte (DHPLC/Sanger → pipeline NGS actuel)
- Le pipeline bioinformatique court-read (alignement, variant calling, annotation)
- La classification ACMG/ClinGen des variants TSC1/TSC2
- L'intégration des résultats dans la base MySQL partagée (`tsc_db.variants`)
- Le lien avec VCF_Interpreter pour l'interprétation interactive

## État actuel du projet

| Étape | État | Notes |
|---|---|---|
| Réanalyse cohort (harmonisation) | En cours | Priorité aux cas avec DHPLC/Sanger uniquement |
| Pipeline NGS bioinformatique | À construire | BWA-MEM2 + GATK/DeepVariant + VEP |
| Classification ACMG | À structurer | Lien avec VCF_Interpreter |
| Import dans MySQL | En attente | Dépend de `03_cohort` (table `patients` d'abord) |

## Pipeline réel : NiourK v2.5 (KAVE - CHU Angers)

Pipeline Nextflow interne au CHU d'Angers. **Ne pas réimplémenter — utiliser les résultats existants.**

```
~/data/TSC_063/          ← métadonnées locales (SampleSheet, validate files)
sl157.chu-angers.intra   ← serveur d'analyse CHU (BAM + VCF complets)
  /data/results/TSC_063/ ← résultats NiourK (VCF annotés Strelka2 + HaplotypeCaller + DeepVariant)
```

**Run TSC_063 :** 22 échantillons (IDs DEFGEN dans SampleSheet), run 2026-01-23, Illumina NextSeq 101 bp PE, panel Twist TSC CHU, alignement GRCh38, VEP v106 + SpliceAI + CADD + REVEL + SPiP.

## Données locales

```
~/data/TSC_063/
├── interop/      ← métriques Illumina (bin + XML)
└── files/        ← NiourK logs + validate "OK" par sample
```

Les VCFs annotés et les BAMs sont uniquement sur `sl157.chu-angers.intra:/data/results/TSC_063/`. Pour travailler dessus : connexion SSH au serveur ou transfert ciblé.

## Connexion à la base MySQL partagée

```python
from sqlalchemy import create_engine
engine = create_engine("mysql+pymysql://user:password@localhost:3306/tsc_db")
# Table cible : variants (patient_id FK → patients.id)
# Vérifier cohérence avec 03_cohort AVANT toute insertion
```

## Environnement

```bash
micromamba activate ont_bioinfo   # Python 3.11
# bcftools, pandas, sqlalchemy pour traitement des VCFs récupérés du serveur
```

## Fichiers de données

- `~/data/TSC_063/` — run TSC_063 (métadonnées locales uniquement)
- VCFs/BAMs complets : serveur CHU `sl157.chu-angers.intra:/data/results/TSC_063/`
- `classification/` — à créer : tableaux de classification ACMG par patient

## Lien avec les autres sous-projets

- **`03_cohort/`** : source de vérité pour les patients (`patients.id`). Ne jamais insérer dans `variants` sans `patient_id` valide.
- **`~/projets/VCF_Interpreter/`** : interface web pour révision interactive des variants classifiés.
- **`06_thesis_writing/`** : les résultats alimenteront les sections 2.4 et 3.x de la thèse.

## Ne pas charger dans le contexte

- `data/raw/` — FASTQs bruts (fichiers lourds)
- `data/_archive/` — versions intermédiaires
