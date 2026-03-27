# Design Spec — vcf_filter.py
**Date :** 2026-03-23
**Projet :** PRELUDE-TSC / 07_ngs_reanalysis
**Objectif :** Filtrer les XLSX annotés NiourK v2.5 (TSC_063) pour identifier les variants candidats pathogènes TSC1/TSC2 à réviser manuellement (cible : ≤ 25 variants par tier).

---

## 1. Vue d'ensemble

```
TSC_063_XXXXX_SX.xlsx  +  config.yaml
              ↓
        vcf_filter.py
              ↓
  [1] Pre-report terminal — comptage par filtre et par tier
              ↓  (ajuster config.yaml si besoin, relancer)
  [2] TSC_063_XXXXX_SX_filtered.xlsx
        ├── _summary       ← miroir du pre-report + métadonnées
        ├── CLINVAR_LOVD
        ├── HIGH
        ├── MODERATE
        ├── SPLICING
        └── MOSAIC
```

**Usage :**
```bash
python vcf_filter.py TSC_063_346769_S1.xlsx --config config.yaml
```

---

## 2. Gènes cibles

- `SYMBOL` ∈ {**TSC1**, **TSC2**} uniquement

> **Multi-transcrit :** NiourK génère une ligne par transcrit par variant. Après filtrage, une déduplication par (CHR, POS, REF, ALT) est appliquée pour compter les variants uniques dans le pre-report. Tous les rangs sont conservés dans les onglets de sortie.

---

## 3. Filtres d'exclusion communs (tous les tiers)

Appliqués séquentiellement en premier. Un variant exclu ici n'apparaît dans aucun onglet.

| Ordre | Filtre | Condition d'exclusion | Raison |
|-------|--------|-----------------------|--------|
| 1 | Population gnomAD | `MAX_AF_GNOMADEX2.1` > 0.005 | Hard exclusion : variant trop fréquent (> 0.5 %) |
| 2 | Cohorte globale | `RNCNT_TOTAL` > 10 | Artefact/variant récurrent tous projets |
| 3 | Cohorte projet en cours | `PJCNT_TOTAL` > 4 | Trop récurrent dans le run TSC_063 |
| 4 | Homozygote cohorte | `RNCNT_HOM` > 0 | Variant déjà vu HOM dans la cohorte → incompatible avec pathogénicité TSC dominante |
| 5 | ClinVar bénin | `CLINVAR` contient "Benign" ou "Likely_benign" | Déjà classifié bénin |
| 6 | Artefact HTZ → low-AF | voir §3.1 | Faux mosaïque |
| 7 | Artefact mosaïque récurrente | voir §3.1 | Hotspot technique |

> **Valeurs manquantes :**
> - `MAX_AF_GNOMADEX2.1` null (variant absent de gnomAD) → **ne pas exclure** (variant novel, conserver)
> - `RNCNT_HOM` null → traiter comme 0 (ne pas exclure)
> - `CLINVAR` null → ne pas exclure

### 3.1 Artefacts par comportement inhabituel

> **Note :** La règle HOM→HTZ n'est pas nécessaire. Le filtre 4 (`RNCNT_HOM > 0`) exclut déjà tout variant ayant été observé homozygote dans la cohorte, ce qui est biologiquement incompatible avec un variant pathogène TSC dominant. L'artefact HOM→HTZ est donc subsumed par ce filtre.

#### Artefact HTZ → low-AF (faux mosaïque)
Variant habituellement hétérozygote dans la cohorte (RNCNT_HTZ dominant) mais détecté à bas allèle dans l'échantillon courant → faux mosaïque (déséquilibre allélique, région difficile, artefact de couverture).

```
RNCNT_HTZ / (RNCNT_HTZ + RNCNT_LOW) > 0.80
ET RNCNT_HTZ ≥ 3
ET AF_PCT échantillon < 35 %        ← AF_PCT en pourcentage (ex: 12.5 = 12.5 %)
→ EXCLU
```

> Valeurs manquantes : si RNCNT_HTZ ou RNCNT_LOW null → ne pas appliquer cette règle.

#### Artefact mosaïque récurrente (hotspot technique)
Variant récurrent à bas allèle dans la cohorte → hotspot d'artefact de séquençage.

```
RNCNT_LOW > 5
ET AF_PCT échantillon < 35 %
→ EXCLU
```

> **Définition RNCNT_LOW (NiourK) :** nombre d'échantillons de la cohorte où ce variant a été détecté avec AF < 35 % (mosaïques et bas allèles). Seuil NiourK interne = 35 % (cohérent avec la fenêtre mosaïque utilisée ici).

---

## 4. Tiers de filtrage

Tous les tiers partagent les filtres d'exclusion communs (§3). Un variant **peut apparaître dans plusieurs onglets** (intentionnel — chaque onglet a un angle d'analyse différent).

### Tier 1 — CLINVAR_LOVD (priorité haute)
Variants déjà annotés dans des bases expertes. Révision prioritaire avant les autres tiers.

- `CLINVAR` contient "Pathogenic" ou "Likely_pathogenic"
- **OU** `LOVD` non vide (non null, non ".")

> Pas de filtre sur les scores in silico pour ce tier — l'annotation experte prime.

### Tier 2 — HIGH
Variants à fort impact prédit.

- `IMPACT` = HIGH *(stop_gained, frameshift_variant, splice_acceptor_variant, splice_donor_variant, start_lost)*
- **OU** `SpliceAI_num` ≥ 0.50
- **OU** `SPiP_Interpretation` = "Altered" **ET** `SPiP_InterConfident` = "Yes"

### Tier 3 — MODERATE
Variants faux-sens et indels in-frame avec support de scores in silico.

- `IMPACT` = MODERATE
- **ET** (au moins un critère de score) :
  - Si variant faux-sens (`Consequence` contient "missense") : `REVEL` ≥ 0.75
  - Si REVEL absent (null) ou variant non faux-sens : `CADD_PHRED` ≥ 20
- Population plus stringente pour ce tier : `MAX_AF_GNOMADEX2.1` < 0.001 (0.1 %)

> Valeurs manquantes : REVEL null → brancher sur CADD. CADD null ET REVEL null → variante sans score, **ne pas inclure dans MODERATE**.

### Tier 4 — SPLICING
Capture tous les variants avec potentiel impact sur l'épissage, toutes conséquences confondues (synonymes, introniques profonds, faux-sens).

- `SpliceAI_num` ≥ 0.20
- **OU** `SPiP_Interpretation` = "Altered" (avec ou sans Confident, pour sensibilité maximale)

> Les variants déjà dans HIGH (SpliceAI ≥ 0.50) apparaîtront aussi ici — intentionnel. Le pré-rapport affichera le compte après déduplication.

### Tier 5 — MOSAIC
Vraies mosaïques candidates (variants somatiques à bas allèle).

- `AF_PCT` échantillon entre **5 % et 35 %**
- `RNCNT_LOW` ≤ 5 (rare dans la cohorte → pas un hotspot technique)
- `RNCNT_HTZ` < 3 (pas un variant habituellement HTZ → exclut l'artefact HTZ→low-AF)

> Les variants à fort impact (stop, frameshift) peuvent aussi apparaître en MOSAIC si leur AF est < 35 % — cette duplication est intentionnelle : l'état mosaïque change l'interprétation clinique.

---

## 5. Colonnes de sortie

Définies dans `output_columns` du `config.yaml`. Modifiables sans toucher au code.

| Groupe | Colonnes |
|--------|----------|
| **Identité** | SYMBOL, HGVSc, HGVSp, CHR, POS, REF, ALT, REF_TRANS, REGIONS, Feature |
| **Qualité / artefact** | AF_PCT, DP, SB, CALLNB, CALLAF_mutect2\|strelka\|gatkHC\|deepvariant, CALLAD_mutect2, CALLAD_strelka, CALLAD_gatkHC, CALLAD_deepvariant, CALLFILTER_mutect2, CALLFILTER_strelka, CALLFILTER_gatkHC, CALLFILTER_deepvariant |
| **Fonctionnel** | IMPACT, Consequence |
| **Scores** | REVEL, CADD_PHRED, SpliceAI_num, SPiP_Interpretation, SPiP_InterConfident, SIFT, PolyPhen |
| **Population** | GNOMAD3_AF, GNOMADEX2.1_AF, MAX_AF_GNOMADEX2.1 |
| **Cohorte** | RNCNT_TOTAL, PJCNT_TOTAL, RNCNT_LOW, PJCNT_LOW, RNCNT_HTZ, PJCNT_HTZ, RNCNT_HOM, PJCNT_HOM |
| **Clinique** | CLINVAR, LOVD, Existing_variation, PUBMED |
| **Travail** | USER_CLASS, USER_ANNOT, USER_COM |

> **CALLFILTER_* :** colonnes par caller (format texte, ex: "PASS", "LowGQX,NoPassedVariantGTs"). Affichées pour détecter les variants rejetés par un ou plusieurs callers — indicateur artefact.

---

## 6. config.yaml

```yaml
# vcf_filter config — PRELUDE-TSC
# AF_PCT est en POURCENTAGE dans les XLSX NiourK (ex: 12.5 = 12.5 %)
# MAX_AF_GNOMADEX2.1 est en FRACTION (ex: 0.005 = 0.5 %)

genes: [TSC1, TSC2]

exclusions:
  gnomad_max_af: 0.005           # fraction — hard exclusion > 0.5 %
  rncnt_max_total: 10            # cohorte globale
  pjcnt_max_total: 4             # projet en cours (run TSC_063)
  exclude_hom: true              # RNCNT_HOM > 0 → exclu (incompatible TSC dominant)
  exclude_clinvar_benign: true

  htz_to_mosaic_artifact:        # HTZ habituel vu en low-AF = faux mosaïque
    enabled: true
    min_htz_ratio: 0.80          # fraction des occurrences cohorte qui sont HTZ
    min_htz_count: 3
    af_max_pct: 35.0             # AF_PCT < 35 % dans l'échantillon (en %)

  mosaic_artifact:               # mosaïque récurrente = hotspot technique
    enabled: true
    rncnt_low_max: 5
    af_max_pct: 35.0             # AF_PCT < 35 %

tiers:
  clinvar_lovd:
    enabled: true

  high:
    enabled: true
    spliceai_threshold: 0.50     # fraction

  moderate:
    enabled: true
    gnomad_max_af: 0.001         # fraction — plus stringent pour ce tier (0.1 %)
    revel_min: 0.75
    cadd_min: 20                 # utilisé si REVEL null ou variant non faux-sens

  splicing:
    enabled: true
    spliceai_min: 0.20           # fraction

  mosaic:
    enabled: true
    af_min_pct: 5.0              # AF_PCT en %
    af_max_pct: 35.0
    rncnt_low_max: 5
    rncnt_htz_max: 3

# Seuils de référence ACMG/ClinGen SVI (documentation — cohérents avec les tiers ci-dessus)
# REVEL  : PP3 ≥ 0.75 | PP3_Fort ≥ 0.932 | BP4 ≤ 0.15   (faux-sens uniquement)
# CADD   : PP3 ≥ 20   | BP4 ≤ 15                          (tous types, surtout non-codant)
# SpliceAI : PP3 ≥ 0.50 (alerte 0.20) | BP4 ≤ 0.02       (tous types)
# SPiP   : PP3 = "Altered" + Confident="Yes" | BP4 = "Normal"
# SIFT/PolyPhen : indicatif uniquement (obsolètes pour ACMG)

output_columns:
  - SYMBOL
  - HGVSc
  - HGVSp
  - CHR
  - POS
  - REF
  - ALT
  - REF_TRANS
  - REGIONS
  - Feature
  - AF_PCT
  - DP
  - SB
  - CALLNB
  - IMPACT
  - Consequence
  - REVEL
  - CADD_PHRED
  - SpliceAI_num
  - SPiP_Interpretation
  - SPiP_InterConfident
  - SIFT
  - PolyPhen
  - GNOMAD3_AF
  - GNOMADEX2.1_AF
  - MAX_AF_GNOMADEX2.1
  - RNCNT_TOTAL
  - PJCNT_TOTAL
  - RNCNT_LOW
  - PJCNT_LOW
  - RNCNT_HTZ
  - PJCNT_HTZ
  - RNCNT_HOM
  - PJCNT_HOM
  - CLINVAR
  - LOVD
  - Existing_variation
  - PUBMED
  - USER_CLASS
  - USER_ANNOT
  - USER_COM
```

---

## 7. Pre-report (terminal)

Les comptages sont **cumulatifs** : chaque ligne représente les variants retirés parmi ceux qui n'avaient pas encore été exclus.

```
=== Pre-report TSC_063_346769_S1 — 2026-03-23 ===
config : config.yaml

Variants TSC1/TSC2 bruts (feuille Variants) : 342
  dont variants uniques (CHR+POS+REF+ALT)   : 198

--- Exclusions communes (séquentielles) ---
Après gnomAD > 0.5 %              :  -70  → 128 restants
Après RNCNT_TOTAL > 10            :  -18  → 110 restants
Après PJCNT_TOTAL > 4             :   -5  → 105 restants
Après RNCNT_HOM > 0               :   -8  →  97 restants
Après ClinVar Benign              :   -3  →  94 restants
Après artefact HTZ→low-AF         :   -4  →  90 restants
Après artefact mosaïque récurrente:   -3  →  87 restants

--- Tiers (sur 87 variants uniques restants) ---
CLINVAR_LOVD  :  5 variants
HIGH          : 12 variants
MODERATE      : 18 variants
SPLICING      :  9 variants  (dont 4 déjà dans HIGH)
MOSAIC        :  3 variants

→ Génération de l'Excel ? (o/n)
```

---

## 8. Onglet `_summary`

Miroir du pre-report, intégré dans le fichier Excel de sortie pour traçabilité. Contient :
- Nom du fichier source
- Date de génération
- Version du config.yaml utilisé (hash MD5 ou chemin absolu)
- Tableau des comptages (identique au pre-report)

---

## 9. Dépendances Python

```
openpyxl     # lecture/écriture XLSX
pandas       # filtrage et manipulation
pyyaml       # config.yaml
hashlib      # hash du config pour traçabilité (stdlib)
```

Environnement : `micromamba activate ont_bioinfo` (Python 3.11)
