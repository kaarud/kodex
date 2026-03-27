# Interprétation des BAM NGS court-read dans IGV
## Guide de revue des variants TSC1/TSC2 — PRELUDE-TSC

---

## Sources

| Ressource | Référence | Accès |
|---|---|---|
| IGV Desktop User Guide | Robinson et al. — documentation officielle | https://igv.org/doc/desktop/ |
| IGV : article fondateur | Robinson J.T. et al. (2011) *Nat Methods* 8:382–384 | doi:10.1038/nmeth.1903 |
| IGV 2.0 update | Thorvaldsdóttir H. et al. (2013) *Brief Bioinform* 14:178–192 | doi:10.1093/bib/bbs017 |
| GATK — variant evaluation in IGV | Broad Institute Best Practices | https://gatk.broadinstitute.org (chercher "IGV" dans la doc) |
| SAMtools/SAM format spec | Li H. et al. (2009) *Bioinformatics* 25:2078–2079 | doi:10.1093/bioinformatics/btp352 |
| Artefacts NGS revue | Meacham F. et al. (2011) *BMC Bioinformatics* 12:451 | doi:10.1186/1471-2105-12-451 |

---

## 1. Vue d'ensemble des 3 BAMs du pipeline NiourK

Chaque variant est visualisé sur 3 pistes dans l'ordre :

| Piste | Contenu | Utilité |
|---|---|---|
| **RAW** | BAM brut post-alignement, avant tout filtre caller | Signal de base : le variant est-il présent dans les reads bruts ? |
| **MUTECT2** | BAM après filtres Mutect2 | Confirme si Mutect2 retient les reads porteurs |
| **CHIM** | BAM chimeric reads (reads avec alignement split/discordant) | Indique si le variant est lié à un réarrangement structural |

---

## 2. Captures générées par igv_capture.py

Pour chaque variant, 3 captures sont produites :

| Fichier | Mode IGV | Ce qu'on cherche |
|---|---|---|
| `_strand.png` | **COLOR_BY READ_STRAND** (rouge = forward, bleu = reverse) | Répartition équilibrée des brins → variant réel |
| `_softclip.png` | **SHOW_SOFT_CLIPPED true** + pas de coloration | Localisation des bases non alignées |
| `_squished.png` | Vue squished sans coloration ni soft-clips | Vue globale de la couverture et du VAF |

---

## 3. Critères de lecture — variant réel vs artefact

### 3.1 Biais de brin (strand bias) — capture `_strand.png`

| Observation | Interprétation |
|---|---|
| Reads porteurs **sur les deux brins** (rouge + bleu) | ✅ Favorable : signal réel |
| Reads porteurs **sur un seul brin uniquement** | ⚠️ Artefact probable (oxydation, ligation asymétrique) |
| Ratio brin majeur / total > 90 % | ⚠️ Alerte forte — vérifier SB dans le XLSX |

La métrique **SB (Strand Bias)** dans le XLSX NiourK quantifie cela. Valeur proche de 0.5 = équilibré = favorable.

### 3.2 Soft-clips — capture `_softclip.png`

Les soft-clips (bases en gris clair hors alignement) apparaissent en extension des reads.

| Observation | Interprétation |
|---|---|
| Soft-clips **dispersés**, pas concentrés sur le variant | Normal |
| Soft-clips **concentrés exactement sur le variant** | ⚠️ Misalignment local — le variant peut être un artefact de bord de read |
| Soft-clips présents sur tous les reads à la même position | ⚠️ Indel complexe ou délétion mal alignée — vérifier avec BAM CHIM |
| Absence totale de soft-clips sur reads porteurs | ✅ Favorable |

### 3.3 Mapping quality — sur toutes les captures

IGV colore les reads par mapping quality (MQ) en mode par défaut :
- **Reads opaques / bien colorés** : MQ ≥ 30 → alignement unique fiable
- **Reads pâles / grisés** : MQ < 20 → reads multimapping, région répétée
- **Reads blancs** : MQ = 0 → non unique, alignement ambigu

> Si les reads porteurs du variant sont pâles : probable artefact de région répétée.

### 3.4 Concordance multi-callers

Vérifier dans le XLSX : colonnes `CALLFILTER_mutect2 / strelka / gatkHC / deepvariant`

| CALLNB | Interprétation |
|---|---|
| 4/4 callers | ✅ Signal robuste |
| 3/4 callers | ✅ Acceptable |
| 2/4 callers | ⚠️ Borderline — regarder IGV attentivement |
| 1/4 caller | ⚠️ Signal faible — artefact probable sauf si mosiaque |
| `PASS` dans CALLFILTER | Caller a passé ses propres filtres internes |

### 3.5 VAF et profondeur

| Observation | Interprétation |
|---|---|
| AF_PCT 40–60 % | Hétérozygote germinal classique |
| AF_PCT 90–100 % | Homozygote ou hémizygote (vérifier RNCNT_HOM) |
| AF_PCT 5–35 % + reads équilibrés | Probable mosaïque — voir onglet MOSAIC |
| AF_PCT < 5 % | Mosaïque profonde ou artefact — grande prudence |
| DP < 30× | Couverture insuffisante — interprétation difficile |

---

## 4. Signaux spécifiques aux variants TSC1/TSC2

### Variants de type splice (HGVSc contenant `+` ou `-`)

- Vérifier que les reads porteurs **chevauchent effectivement la jonction exon-intron**
- Soft-clips à la jonction = normal pour un variant splice donor/acceptor réel
- Vérifier SpliceAI_num et SPiP_Interpretation dans le XLSX

### Variants de type indel (del / ins / dup)

- Les indels génèrent souvent des soft-clips sur les reads adjacents → **normal**
- Vérifier que l'indel est visible dans le BAM RAW **et** MUTECT2
- Si visible dans CHIM seulement → probable réarrangement structural (exclure délétion exonique)

### Variants mosaïques (AF 5–35 %)

- Demandent une couverture ≥ 100× pour être fiables
- Le bilan des 3 callers est particulièrement important : un vrai mosaïque est souvent appelé par 2–3 callers
- RNCNT_LOW dans le XLSX = nombre de runs où ce variant a été vu à bas AF

---

## 5. Workflow de revue IGV recommandé

```
1. Ouvrir IGV Desktop, charger hg38
2. Charger les 3 BAMs (RAW + MUTECT2 + CHIM)
3. Lancer igv_bookmark TSC_063_*_filtered.xlsx
4. Ouvrir Regions > Region Navigator
5. Pour chaque variant :
   a. Vérifier brin (capture _strand) → biais ?
   b. Vérifier soft-clips (capture _softclip) → concentration sur variant ?
   c. Vérifier MQ des reads porteurs → pâles ?
   d. Comparer RAW vs MUTECT2 vs CHIM
   e. Croiser avec CALLNB, SB, DP dans le XLSX
6. Renseigner USER_CLASS + USER_ANNOT + USER_COM dans le XLSX filtré
```

---

## 6. Codes USER_CLASS recommandés (ACMG simplifié)

| Classe | Signification |
|---|---|
| 5 | Pathogène — variant confirmé TSC |
| 4 | Probablement pathogène |
| 3 | Variant de signification incertaine (VUS) |
| 2 | Probablement bénin |
| 1 | Bénin / artefact confirmé |
| M | Mosaïque — à préciser en USER_COM |

---

*Document interne PRELUDE-TSC — mise à jour 2026-03-24*
