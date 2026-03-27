# Classification des variants TSC1/TSC2 — Cadre actuel et évolution
## Focus : fréquence comme critère « supporting » et concept VSI+

*Document interne PRELUDE-TSC — 2026-03-24*

---

## Sources citées

| Référence | DOI / Accès |
|---|---|
| Richards et al. 2015 — Cadre ACMG/AMP | doi:10.1038/gim.2015.30 |
| Tavtigian et al. 2018 — Cadre bayésien | doi:10.1038/gim.2017.210 |
| Tavtigian et al. 2020 — Système de points | doi:10.1002/humu.24088 |
| Ghosh et al. 2018 — Révision critère BA1 | doi:10.1002/humu.23642 |
| Brnich et al. 2019 — Critère PS3/BS3 (fonctionnel) | doi:10.1186/s13073-019-0690-2 |
| Leon-Quintero et al. 2025 — Mosaïcisme somatique | doi:10.1111/cge.14636 |
| ANPGM BP-NGSDiag_001 v2 — Recommandations françaises | https://anpgm.fr/media/documents/BP-NGSDiag_001_Interpretation_Variants_v2.pdf |
| ClinGen SVI Working Group | https://clinicalgenome.org/working-groups/sequence-variant-interpretation/ |
| ClinGen VCEPs (liste) | https://clinicalgenome.org/affiliation/vcep/ |

---

## 1. Cadre de référence : ACMG/AMP 2015

Le cadre Richards et al. (2015) définit 5 classes de variants et un système de critères combinables :

```
Classe 5 — Pathogène (P)
Classe 4 — Probablement pathogène (LP)
Classe 3 — Variant de signification incertaine (VUS / VSI)
Classe 2 — Probablement bénin (LB)
Classe 1 — Bénin (B)
```

Les critères sont organisés en **forces** : *Very Strong (PVS1)*, *Strong (PS1-4)*, *Moderate (PM1-6)*, *Supporting (PP1-5)* du côté pathogène ; miroir du côté bénin (BA1, BS1-4, BP1-7).

La classification finale résulte de la **combinaison** de ces critères selon des règles de score définies.

---

## 2. Évolution vers un cadre bayésien (Tavtigian 2018–2020)

### 2.1 Le problème des règles combinatoires rigides

Le cadre 2015 pose des règles combinatoires binaires (ex. : « 1 Strong + 3 Supporting = LP »). Deux limites majeures :

1. **Perte d'information** : un variant à 3 critères Supporting + 1 Moderate n'est pas équivalent à un variant à 4 Supporting, mais les deux peuvent tomber dans la même case.
2. **Manque de gradation** : impossible d'exprimer un variant « presque LP » vs « au fond du VUS ».

### 2.2 Le cadre bayésien (Tavtigian et al. 2018, doi:10.1038/gim.2017.210)

Tavtigian et al. modélisent les critères ACMG comme des **ratios de vraisemblance (LR)** additifs dans un cadre bayésien. Chaque critère apporte un poids :

| Force | LR pathogène | Exemple |
|---|---|---|
| Very Strong (PVS1) | ~350 | Variant perte de fonction sur gène haploinsuffisant |
| Strong (PS) | ~18.7 | De novo confirmé (PS2) |
| Moderate (PM) | ~4.3 | Absent de gnomAD (PM2) |
| Supporting (PP) | ~2.1 | Score splicing élevé (PP3) |

La probabilité finale de pathogénicité est calculée à partir d'une **probabilité a priori** (prévalence de la maladie) combinée à tous les LR.

### 2.3 Le système de points (Tavtigian et al. 2020, doi:10.1002/humu.24088)

Pour simplifier l'usage clinique, le groupe propose une **échelle de points** :

| Force | Points pathogène | Points bénin |
|---|---|---|
| Very Strong | +8 | −8 |
| Strong | +4 | −4 |
| Moderate | +2 | −2 |
| Supporting | +1 | −1 |

**Seuils de classification par points :**

| Score total | Classe |
|---|---|
| ≥ 10 | Pathogène (5) |
| 6 à 9 | Probablement pathogène (4) |
| 0 à 5 | VUS (3) |
| −1 à −6 | Probablement bénin (2) |
| ≤ −7 | Bénin (1) |

---

## 3. La révision du critère de fréquence : BA1 et PM2

### 3.1 Le problème du BA1 « standalone »

Dans le cadre 2015, **BA1 (allele frequency > 5 %)** permet à lui seul de classer un variant en Bénin (classe 1), sans autre critère. Ghosh et al. (2018, doi:10.1002/humu.23642) identifient des cas où cette règle est trop forte :

- Variants fréquents dans certaines populations mais pathogènes dans d'autres contextes
- Faux positifs de gnomAD (erreurs d'annotation, pseudogènes)
- Gènes à pénétrance incomplète ou à expression variable

**Résultat : BA1 n'est plus standalone par défaut.** Le groupe SVI recommande que chaque VCEP définisse ses propres seuils de fréquence, et que BA1 soit appliqué avec le contexte gène-spécifique. Pour certains gènes, BA1 est ramené à une force « Supporting » (BS1_supporting).

### 3.2 Downgrade de PM2 : de Moderate à Supporting

Parallèlement, plusieurs VCEPs ont rétrogradé **PM2** (absent ou extrêmement rare dans gnomAD) de *Moderate* à *PM2_supporting*.

**Justification :** L'absence dans gnomAD n'est pas en soi informative sur la pathogénicité d'un variant — elle reflète surtout la rareté ou l'absence de données. Avec le système de points :

- PM2 Moderate → **+2 points**
- PM2 Supporting → **+1 point seulement**

### 3.3 Conséquence directe sur la classification

Un variant avec uniquement des critères de rareté (PM2) et de prédiction computationnelle (PP3 : CADD élevé, SpliceAI) :

**Avant :**
```
PM2 (Moderate, +2) + PP3 (Supporting, +1) = 3 points → VUS
Mais avec PM2_Moderate + d'autres critères → pouvait atteindre LP
```

**Après révision (PM2_supporting) :**
```
PM2_sup (+1) + PP3 (+1) = 2 points → reste VUS
Même avec 3 PP criteria (+3) = 5 points → reste VUS (seuil LP = 6)
```

Des variants qui atteignaient auparavant **6–7 points** (LP) tombent maintenant à **4–5 points** (VUS) uniquement par le downgrade de PM2.

---

## 4. Le concept de VSI+

### 4.1 Définition

**VSI+** (ou VUS+, ou classe 3+) est une **notation informelle**, utilisée dans certains laboratoires français et documentée dans le cadre ANPGM, pour désigner :

> Un variant de classe 3 (VUS) dont le score de points est élevé dans la zone VUS (typiquement 4–5 points), suggérant une probabilité de pathogénicité supérieure à la médiane des VUS, sans atteindre le seuil probablement pathogène (6 points).

⚠️ **Cette notation n'est pas standardisée par ClinGen, l'ACMG ou l'AMP.** Il n'existe pas de publication de référence avec DOI la définissant. Elle est utilisée pour la communication clinicien-biologiste mais ne doit pas figurer dans les bases de données sans annotation explicite de son caractère non-standard.

### 4.2 Pourquoi ce concept est apparu

Le downgrade de PM2 (§3.3) a créé une zone « grise haute » dans le VUS :

```
Score 0-2 : VUS « froid » — peu d'arguments
Score 3-5 : VUS « chaud » / VSI+ — arguments insuffisants mais réels
Score 6+  : Probablement pathogène (LP)
```

Le VSI+ traduit l'inconfort clinique face à un variant à 4–5 points : pas assez pour recommander une action clinique basée sur le variant seul, mais suffisamment suspect pour déclencher une surveillance ou des analyses complémentaires (étude parentale, études fonctionnelles).

### 4.3 Critères qui « font » un VSI+ en pratique (TSC1/TSC2)

| Critère | Points | Commentaire |
|---|---|---|
| PM2_supporting (absent gnomAD) | +1 | Presque toujours présent pour variant rare |
| PP3 (CADD ≥ 25, REVEL ≥ 0.5, SpliceAI ≥ 0.2) | +1 | Prédicteurs computationnels convergents |
| PP2 (gène à contrainte missense) | +1 | TSC1 et TSC2 ont scores de contrainte élevés |
| PM1 (domaine fonctionnel critique) | +2 | GAP domain TSC1/TSC2, domaine HEAT |
| BP4 négatif | 0 | Si absent : pas de pénalisation |

Un variant avec PM2_sup + PP3 + PP2 = **3 points** = VUS froid.
Avec PM1 en plus = **5 points** = VUS chaud / VSI+.

---

## 5. Classification des variants TSC1/TSC2 (STB) — État actuel

### 5.1 Absence de TSC-VCEP officiel

Au moment de la rédaction (mars 2026), **ClinGen n'a pas publié de VCEP TSC1/TSC2**. La classification repose donc sur le cadre ACMG/AMP général avec des adaptations de la littérature TSC.

### 5.2 Critères à fort impact pour TSC1/TSC2

#### PVS1 — Perte de fonction (Very Strong, +8 pts)

S'applique aux **variants non-sens, frameshift, délétion/duplication exonique, variants canoniques de site d'épissage (±1,+2)** de TSC1 et TSC2.

Conditions pour PVS1 plein :
- TSC1 et TSC2 sont des gènes suppresseurs de tumeur avec **haploinsuffisance documentée** comme mécanisme pathogène
- Le variant doit se trouver **hors de la dernière exon et hors des 50 derniers nucléotides du dernier intron** (NMD attendu)
- Vérifier l'absence d'isoforme fonctionnel alternative

> ⚠️ PVS1 doit être nuancé (PVS1_strong, PVS1_moderate) si le variant se situe en 3' du gène où le NMD est incertain ou si des données fonctionnelles sont absentes.

#### PS2/PM6 — Variant de novo (Strong / Moderate)

- **PS2** (+4 pts) : de novo confirmé par analyse parentale (duo/trio)
- **PM6** (+2 pts) : de novo présumé (non confirmé par parentale)

Dans TSC, ~60–70 % des variants pathogènes sont de novo. La confirmation parentale est donc un **critère de haute valeur** à systématiser.

#### PS3/BS3 — Études fonctionnelles (Brnich et al. 2019, doi:10.1186/s13073-019-0690-2)

- PS3 (+4 pts) : le variant abolit la fonction TSC1/TSC2 (activation mTORC1, test de complémentation, étude de protéine tronquée)
- BS3 (−4 pts) : le variant ne modifie pas la fonction

Les études fonctionnelles TSC disponibles couvrent principalement les variants missense et de la région GAP de TSC2.

#### PM2_supporting — Rareté dans gnomAD (+1 pt)

Applicable à quasi tous les variants pathogènes TSC1/TSC2 (absents ou < 0.001 % dans gnomAD). Après downgrade, contribue **+1 point uniquement**.

#### PP3 — Prédicteurs computationnels (+1 pt)

Convergence de plusieurs outils requis (pas un seul outil seul) :
- CADD PHRED ≥ 25
- REVEL ≥ 0.5 (pour missense)
- SpliceAI ≥ 0.5 ou SPiP_Interpretation = "Altered" (pour variants splice)

#### PP1/BS4 — Ségrégation familiale

- PP1_strong (+4 pts) si ≥ 7 méioses informatives
- PP1_moderate (+2 pts) si 5–6 méioses
- PP1_supporting (+1 pt) si 3–4 méioses
- BS4 (−4 pts) si non-ségrégation dans une famille

### 5.3 Cas spécifiques aux variants profonds introniques

TSC1 et TSC2 présentent des variants pathogènes en dehors des sites canoniques d'épissage (variants dits « deep intronic »), activant des pseudo-exons.

Classification :
- Utiliser SpliceAI_num ≥ 0.5 comme PP3 (+1 pt)
- SPiP_Interpretation = "Altered" comme PP3 (+1 pt) — **ne pas cumuler SpliceAI + SPiP comme deux PP3 séparés**
- La confirmation par RT-PCR sur ARN lymphocytaire ou fibroblastes = PS3 (+4 pts) si positive
- En l'absence de confirmation ARN : rester VUS (VSI+ si PM2_sup + PP3 + PP2 = 3–4 pts)

### 5.4 Cas spécifiques aux variants mosaïques

Référence : Leon-Quintero et al. 2025 (doi:10.1111/cge.14636).

Les variants mosaïques détectés par NGS (AF 5–35 %) suivent des règles modifiées :
- Un VAF ≥ 5 % en NGS panel constitutionnel peut être considéré comme **variant somatique mosaïque de novo** (équivalent PS2 adapté)
- Le niveau de preuve est réduit si le VAF < 10 % (confirmation par ddPCR ou amplicon deep-sequencing recommandée)
- La pathogénicité d'un variant mosaïque **TSC** peut être établie indépendamment de sa présence en germinal si le VAF est > 15 % avec confirmation multi-techniques

---

## 6. Résumé pratique pour la revue PRELUDE-TSC

| Type de variant | Critères clés | Classe probable |
|---|---|---|
| Stop / frameshift exonique | PVS1 + PM2_sup | P (5) dès confirmation PVS1 |
| Splice canonique ±1,±2 | PVS1 + PM2_sup | P (5) |
| Splice region ±3–8 | PP3 (SpliceAI/SPiP) + PM2_sup + PP2 | VSI (3) à LP (4) selon score |
| Deep intronic + SpliceAI ≥ 0.5 | PP3 + PM2_sup | VSI+ (3+) — confirmation ARN requise pour LP |
| Missense + domaine fonctionnel + CADD ≥ 25 | PM1 + PM2_sup + PP3 + PP2 | VSI+ (3+) à LP (4) |
| De novo confirmé + PVS1 | PS2 + PVS1 | P (5) |
| Mosaïque AF 15–35 % | Critères adaptés Leon-Quintero 2025 | LP à P selon VAF et confirmation |
| Synonyme sans impact splice | BP7 + BP4 | LB (2) |

---

*Aucun VCEP TSC n'étant disponible (ClinGen, mars 2026), ce guide s'appuie sur le cadre ACMG/AMP général avec les révisions SVI publiées. Les critères devront être mis à jour lors de la publication d'un TSC-VCEP officiel.*
