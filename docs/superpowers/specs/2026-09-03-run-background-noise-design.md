# Bruit de fond du run — mesure du contraste inter-échantillons

**Date :** 2026-09-03
**Statut :** spec, rien d'implémenté
**Portée :** TSC1 / TSC2 uniquement (voir *Hors périmètre*)

## Le problème

Sur la cohorte run 1, la vérification manuelle en IGV élimine la grande majorité
des variants candidats : ils sont retrouvés chez un témoin du run. Le motif dominant
est un désaccord de variant callers sur un bruit **partagé** — un caller appelle chez
un échantillon et pas chez l'autre, alors que les reads alternatifs sont présents
chez les deux.

Aucune colonne du XLSX NiourK ne capture ce mécanisme. `RNCNT_*` et `PJCNT_*` comptent
les échantillons où un variant a été **appelé** ; ils sont aveugles au bruit qu'aucun
caller n'a appelé. Conséquence mesurable : un artefact présent dans tout le run mais
appelé chez un seul patient obtient `RNCNT_TOTAL = 1`, c'est-à-dire le **score de
rareté maximal**. Le scoring actuel récompense ces variants pour la propriété qui
devrait les condamner.

### Les quatre cas de référence

Captures dans `<repo>/`.

| Cas | Position | AF patient | AF témoin S2 | Contraste |
|---|---|---|---|---|
| Variant htz propre | chr9:132 908 526 | T 751/1516 = 50 % | 0/396 | ∞ |
| Artefact CHIM uniquement | chr16:2 076 489 | G 49/1839 = 3 % | 1/740 = 0,1 % | ×30 |
| Artefact del couverture | chr9:132 896 302 | DEL 16/2073 = 0,8 % | DEL 1/547 = 0,2 % | ×4 |
| Artefact called S3 pas S2 | chr9:132 910 597 | DEL 13/552 = 2,4 % | DEL 11/556 = 2,0 % | **×1,2** |

Le contraste est donné ici comme rapport brut des AF observées ; la formule ci-dessous
applique un plancher de profondeur, donc `htz_propre` vaut en pratique ×250 environ et
non l'infini. L'ordre est inchangé.

Le cas à condamner se distingue par un rapport proche de 1. La présence ou l'absence
de reads alt chez les autres ne suffit pas — trois des quatre cas en ont. C'est le
**contraste** qui porte l'information.

## Ce que la mesure fait, et ce qu'elle ne fait pas

**Fait :** quantifier à quel point l'allèle alternatif du patient dépasse le bruit
que le run produit à cette position précise.

**Ne fait pas :** trancher les artefacts de reads chimériques. Le cas
`Artefact_variant_sur_bam_chim_uniquement` obtient ×30 et sera classé candidat par
cette métrique — à raison, puisque le contraste inter-échantillons est réellement
élevé. Ce qui l'accuse (reads porteurs saturés de mismatches, présents surtout dans
le BAM CHIM) relève d'une mesure distincte : qualité des reads porteurs. Elle est
hors périmètre ici et fera l'objet de sa propre spec.

Une colonne qui mesure une chose et la mesure bien vaut mieux qu'une couleur qui
prétend tout résumer.

## La métrique

Pour chaque variant retenu, à sa position exacte, sur chaque BAM du run :

```
compteurs(échantillon) = { A, C, G, T, DEL, INS, DP }      # pileup brut
AF_alt(échantillon)    = compteur[allèle_alt_du_patient] / DP
```

L'allèle comparé est celui appelé chez le patient — un `A>T` se compare aux `T` des
autres, une délétion aux `DEL` des autres. Comparer des événements différents au même
locus relève de la mesure de *désordre*, qui n'est pas celle-ci.

```
bruit    = médiane( AF_alt(e) pour e dans autres_échantillons
                    si AF_alt(e) < genotype_floor )
contraste = AF_alt(patient) / max(bruit, 1 / DP_médian_des_autres)
```

Le plancher au dénominateur évite la division par zéro et borne le contraste à ce que
la profondeur permet réellement de distinguer : avec 500× de couverture, un bruit
mesuré à 0 ne signifie pas « nul » mais « sous 1/500 ». Un site parfaitement propre
donne donc un contraste fini et proportionnel à la profondeur, pas l'infini.

### Médiane, pas maximum

Décision motivée par la cohorte : les échantillons de run 1 sont hétérogènes
(décennies de congélation, concentrations harmonisées au mieux, dernière goutte
utilisée quand c'était la dernière goutte). Le bruit de fond global est élevé et
**observé**, pas supposé.

Prendre le maximum reviendrait à laisser un seul échantillon dégradé condamner un
variant vrai chez tous les autres — perte de sensibilité inacceptable ici. La médiane
est robuste à quelques échantillons pourris.

Le maximum reste **affiché** dans le commentaire : l'écart médiane/max est lui-même
informatif. Médiane 0,1 % avec max 3 % signale un échantillon dégradé ; médiane 2 %
signale un site chroniquement mauvais.

Le choix est exposé en configuration (`statistic: median | max`) pour permettre de
mesurer l'effet du contraire sans toucher au code.

### `genotype_floor`

Un échantillon dont l'AF dépasse ~25 % à cette position ne produit pas du bruit : il
porte le variant. L'inclure dans la médiane écraserait le contraste à ~1 et peindrait
en rouge un variant partagé légitime. Ces échantillons sont exclus du calcul du bruit
et **signalés dans le commentaire** — un porteur parmi 21 TSC non apparentés est un
fait à regarder, pas à masquer.

### Résolution fine où ça compte

Les artefacts observés vivent entre 0,2 % et 3 %. La feuille `DEEP_MOSAIC` cible
0,5 % – 5 %. Bruit de fond et mosaïques profondes occupent le même intervalle d'AF :
l'échelle doit être fine entre ×1 et ×5, et grossière au-delà. Un variant à ×30 et un
à ×100 méritent le même traitement — ils seront ouverts dans les deux cas.

## Intégration au scoring

C'est la cible. Une entrée de plus dans `scoring.technical`, qui alimente le canal
rouge (risque artefact = 1 − qualité technique) :

```yaml
  technical:
    sb:        { weight: 2.0, ... }      # existant
    callnb:    { weight: 2.0, ... }      # existant
    dp_af:     { weight: 0.5, ... }      # existant

    run_background:
      weight: 3.0              # dominant — discriminant observé dans ~99 % des vérifications
      good: 10.0               # contraste ≥ ×10  → score plein
      moderate: 3.0            # ×3 à ×10         → demi-score
      artifact: 1.5            # < ×1.5           → zéro, risque maximal
      statistic: median        # median (sensible) | max (spécifique)
      genotype_floor: 25.0     # AF ≥ 25 % chez un autre = porteur, exclu du bruit
      min_other_samples: 5     # sous ce nombre de BAM lisibles, mesure non rendue
```

Le vocabulaire `weight` / `good` / `moderate` est celui de `cadd`, `revel`, `spliceai`.
Le mécanisme d'accumulation existant s'applique sans modification : `tech_max += w`
seulement si la mesure est présente, donc **un BAM absent ne pénalise pas** — il ne
compte pas.

Le poids de 3,0 (contre 4,5 pour les trois autres réunis) donne au bruit de fond ~40 %
du canal technique. C'est délibérément fort, et c'est un nombre dans un YAML : il sera
recalibré quand les chiffres de validation existeront.

## Dispositif provisoire d'évaluation

Tant que le pouvoir discriminant n'est pas établi, la mesure est rendue visible pour
être jugée — commentaire sur la cellule `DP`, reproduisant le popup IGV :

```
chr9:132 910 597 — S3 (patient)
  G:551  A:1  T:0  C:0   DEL:13  INS:0   DP=552   → DEL 2,4 %

Run run 1 — 21 autres échantillons
  S2   DP=556   DEL:11   → 2,0 %
  S4   DP=498   DEL:9    → 1,8 %
  ...
  médiane 1,9 %  |  max 2,4 %  |  avec reads alt : 19/21
  contraste ×1,2
```

**Pas de coloration de la cellule.** Les chiffres se lisent directement — la médiane,
le max et le contraste sont dans le commentaire, et un seuil de couleur n'ajouterait
rien qu'un risque de contresens (`DP` rouge se lirait « couverture insuffisante »).
La seule couleur reste celle du RGB global, où la mesure entre par le canal technique.

Une distinction à préserver dans le texte du commentaire : « aucun read alt mesuré »
et « rien mesuré » ne sont pas la même chose. Un BAM absent doit s'écrire comme tel,
jamais comme un zéro — sinon il se lit comme une bonne nouvelle.

**Ce dispositif est temporaire.** Si la mesure discrimine, elle vit dans le score YAML
et le commentaire `DP` disparaît. Il n'y aura pas de colonne dédiée permanente.

## Architecture

Module `run_background.py`, suivant le découpage de `minibam.py` — logique pure
séparée de l'exécution :

- **Pur, testable sans BAM ni samtools** : agrégation des compteurs, exclusion
  `genotype_floor`, médiane/max, calcul du contraste, mapping contraste → points de
  score.
- **Impur** : appel `samtools mpileup` sur les BAM du run restreints aux positions
  retenues. Les tranches produites par `minibam.py` conviennent directement.

Le pileup ne filtre ni sur MAPQ ni sur qualité de base par défaut — voir le bruit que
les callers masquent est l'objectif. Les seuils MAPQ/BQ sont des paramètres du module ;
la comparaison brut / filtré est en soi un signal (cas d'un artefact en épingle à cheveux :
VAF 1,8 % → 0,09 % après filtrage MAPQ).

## Validation — comment on saura que ça marche

La décision d'intégrer au score en dépend, donc le critère est posé d'avance.

Jeu d'évaluation : les variants déjà tranchés à la main, `USER_CLASS` / `USER_ANNOT`
faisant office de vérité terrain. Mesures :

1. Distribution du contraste chez les variants annotés artefact vs retenus.
2. Sensibilité à `good = 10` — combien de vrais variants tomberaient sous le seuil.
3. Les quatre cas de référence doivent s'ordonner comme observé en IGV :
   `htz_propre` (∞) > `CHIM` (×30) > `del_couverture` (×4) > `called_S3` (×1,2),
   ce dernier seul sous `artifact = 1.5`.

Le poids `3.0` est provisoire jusqu'à ce que (1) existe.

La validation est techniquement réalisable dès maintenant en local. Ce qui la retarde
est la contrainte de déploiement ci-dessous, pas une indisponibilité des données.

## Contrainte de déploiement — la vraie difficulté

En local, tout est co-localisé : les BAM vivent dans le même dossier parent que les
XLSX, et `resolve_bam_dir()` les trouve par un chemin relatif (`bam_dir: ../BAM`,
ancré sur le XLSX). C'est ce qui rend le lien IGV et `minibam.py` triviaux aujourd'hui.

En production au service de génétique, la topologie est éclatée :

- **Serveur X** — les scripts tournent, les XLSX sont générés. Pas de BAM.
- **Serveur Y** — les biologistes travaillent ; chacun a son dossier, où les BAM
  partagent le même dossier parent que les XLSX.

La relation BAM ↔ XLSX existe donc bien sur Y, mais **pas au moment où le XLSX est
produit**, sur X. Un pileup exécuté pendant `vcf_filter.py` n'aurait aucun BAM à lire.

Conséquence sur l'architecture : la mesure ne peut pas être une étape de la génération.
Elle doit être un **traitement séparé, exécuté là où sont les BAM**, enrichissant un
XLSX déjà produit — le modèle que suivent déjà `aggregate_annotations.py` et
`minibam.py`.

### Le point d'insertion existe déjà

Si `run_background` entre dans `scoring.technical`, le score ne peut plus être complet
sur X. Mais il n'a pas à l'être : `aggregate_annotations.py` importe déjà
`_apply_variant_scores` et `_apply_cell_colors` depuis `vcf_filter`, donc **le score est
recalculé au moment de l'agrégation**, pas figé à la génération. Et `minibam.py` lit
déjà la feuille `ANNOTATIONS`.

La chaîne se referme donc sans architecture nouvelle :

| Étape | Où | Produit |
|---|---|---|
| `vcf_filter.py` + `run_background` | local | XLSX par patient, mesure incluse au score |
| `aggregate_annotations.py` | local | feuille agrégée, score recalculé, colonne héritée |

En local — le cas actuel — les BAM sont co-localisés, donc `vcf_filter.py` peut appeler
la mesure directement et le pipeline reste en une passe. C'est seulement sur la topologie
X/Y éclatée que la mesure doit se scinder en post-traitement :

| Étape | Où | Produit |
|---|---|---|
| `vcf_filter.py` | X | XLSX par patient, score sans `run_background` |
| `aggregate_annotations.py` | X | feuille agrégée |
| `run_background.py` | Y | pileup → enrichit, score et couleurs recalculés |

Le module doit donc être appelable des deux façons : en ligne depuis `vcf_filter.py`
quand les BAM sont là, et en post-traitement sur un XLSX déjà produit quand ils ne le
sont pas. Le mécanisme `tech_max += w` rend la dégradation propre — sans mesure, le
score est simplement calculé sans elle.

Reste le portage du lien IGV et du chargeur de BAM vers les serveurs des
bioinformaticiens : leur ancrage local (`bam_dir` relatif au XLSX) devra être revu. Ce
chantier est commun à ces fonctionnalités et conditionne le calendrier — mais il ne
bloque pas la validation, qui se fait en local.

## Cible : les XLSX par patient, dès la génération

La mesure vit **dès le premier tri**, dans les `*_filtered.xlsx` produits par
`vcf_filter.py`. C'est là que la question « est-ce que ce variant se voit chez les
autres ? » se pose en premier — au moment de décider quoi regarder, pas seulement au
moment de conclure.

`aggregate_annotations.py` en **hérite** : la colonne se propage par `row.to_dict()`
lors de la collecte, et `_apply_variant_scores` recalcule le score sur la feuille
agrégée. Deux ajustements suffisent pour que ce soit visible à la relecture :
la colonne dans `_PRIORITY_COLS`, et le report du commentaire sur `DP`.

Autrement dit, aucune des deux étapes ne détient la mesure en propre — elle est produite
en amont et suivie partout, comme n'importe quelle autre colonne d'annotation.

### Nommer ce à quoi on compare

Principe repris de `reference_label()` : le relecteur doit voir *qui* est le comparateur,
jamais un « REF » anonyme masquant un fallback. La mesure suit la même règle — le
commentaire liste les échantillons interrogés par leur identifiant, jamais un « n=21 »
opaque.

Divergence assumée avec IGV : le lien charge **un** comparateur, parce qu'empiler
21 pistes serait illisible. Le pileup n'a pas cette contrainte et interroge **tous** les
BAM disponibles. C'est nécessaire — sur un comparateur unique, la médiane du bruit
dégénère en la comparaison par paire déjà faite à l'œil, et perd tout intérêt.

## Hors périmètre

- **RHEB et PKD1.** Conservés dans le classeur pour de bonnes raisons — un cas TSC
  publié avec variant RHEB, et les syndromes de délétion contiguë TSC2-PKD1 — mais
  exclus de cette analyse. Ils sont déjà chargés séparément (`vcf_filter.py:1611-1620`,
  feuille `RHEB_PKD1`) : il suffit de ne pas la traiter.
- **Qualité des reads porteurs** (NM des reads alt, proportion de `supplementary`,
  diversité des positions de début). Signal compagnon nécessaire — c'est lui qui
  rattrape le cas CHIM — mais mesure distincte, spec distincte.
- **Panel of Normals inter-runs.** À terme, la version cumulative de cette mesure a
  sa place en amont (NiourK), comme ressource qui se bonifie run après run. Prototyper
  en aval d'abord : demander une modification d'un pipeline clinique accrédité pour un
  signal non validé serait le mauvais ordre.
- **Détection de contexte homopolymère / répétition inversée.**
