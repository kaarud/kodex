# Design — Onglets fiche d'analyse panel STB

**Date:** 2026-03-24
**Statut:** Validé

## Contexte

La fiche d'analyse du panel STB (document 1587.docx) définit un workflow d'interprétation
des variants TSC1/TSC2 en étapes séquentielles. Chaque étape doit produire un **onglet Excel
indépendant** — le chevauchement entre onglets est voulu (argument de pathogénicité).

## Architecture

```
XLSX NiourK (sheet "Variants")
  → load_variants() [filtre TSC1/TSC2]
  → apply_exclusions() [pré-filtre bruit]
  → classify_sheets() [11 onglets indépendants sur survivants]
  → write_output() [Excel multi-onglets]
```

## Pré-exclusions (inchangées)

Appliquées séquentiellement pour éliminer le bruit évident :

| Filtre | Condition | Config key |
|--------|-----------|------------|
| gnomAD | MAX_AF_GNOMADEX2.1 > 0.005 | `exclusions.gnomad_max_af` |
| RNCNT réseau | RNCNT_TOTAL > 10 | `exclusions.rncnt_max_total` |
| PJCNT projet | PJCNT_TOTAL > 4 | `exclusions.pjcnt_max_total` |
| Homozygotes | RNCNT_HOM > 0 | `exclusions.exclude_hom` |
| ClinVar Benign | Benign / Likely_benign | `exclusions.exclude_clinvar_benign` |
| Artefact HTZ→low-AF | ratio HTZ élevé + AF basse | `exclusions.htz_to_mosaic_artifact` |
| Artefact mosaïque | RNCNT_LOW élevé + AF basse | `exclusions.mosaic_artifact` |

## 11 onglets indépendants

Chaque onglet filtre **indépendamment** les variants post-exclusion.

### 1. `HIGH_IMPACT`
- `IMPACT == "HIGH"`

### 2. `CLINVAR_P_LP`
- CLINVAR contient `Pathogenic` ou `Likely_pathogenic`
- Exclut les entrées contenant `Conflicting` (capturées par l'onglet dédié)

### 3. `CLINVAR_CONFLICTING`
- CLINVAR contient `Conflicting`

### 4. `CLINVAR_VUS`
- CLINVAR contient `Uncertain_significance`

### 5. `CLINVAR_NOT_PROVIDED`
- CLINVAR contient `not_provided`

### 6. `FQ_MODERATE`
- `MAX_AF_GNOMADEX2.1 < 0.001` (< 0.1%)
- `IMPACT == "MODERATE"`

### 7. `FQ_LOW`
- `MAX_AF_GNOMADEX2.1 < 0.001`
- `IMPACT == "LOW"`

### 8. `FQ_MODIFIER`
- `MAX_AF_GNOMADEX2.1 < 0.001`
- `IMPACT == "MODIFIER"`
- (`SpliceAI_num >= 0.50` **OU** `SPiP_Interpretation == "Altered"`)

### 9. `4_CALLERS`
- `CALLNB == 4`

### 10. `LOVD`
- `LOVD` non vide et ≠ `"."`

### 11. `MOSAIC`
- `AF_PCT > 5`
- `RNCNT_TOTAL <= 2`

## Modifications fichiers

### `config.yaml`
- Section `tiers` remplacée par section `sheets` avec les seuils de chaque onglet

### `vcf_filter.py`
- `classify_tiers()` → `classify_sheets()` : 11 filtres indépendants
- `build_pre_report()` : affiche les 11 onglets au lieu des 5 tiers
- `write_output()` : écrit les 11 onglets + `_summary`

### Tests
- `test_tiers.py` → adapté aux 11 nouveaux onglets
- Vérifier le chevauchement (un variant peut apparaître dans plusieurs onglets)

## Hors scope (futur)

- **Reads chimériques TSC1/TSC2** : nécessite analyse BAM + IGV. Critère SB : 0.2 < SB < 0.8. À paramétrer ultérieurement.
- **RHEB / Blanket / MLPA** : exclus de cette itération.
- **Communication IGV** : script futur.
