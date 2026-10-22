# Lire le classeur `*_filtered.xlsx`

Guide de lecture de la sortie de `vcf_filter.py`. L'installation et les
commandes sont dans le [README](../README.md) ; les seuils se règlent dans
`config.yaml` (profil courant) ou `config_legacy.yaml` (`--legacy`), sans
toucher au code.

KODEX signale des variants **à revoir** ; il ne pose pas de diagnostic.

## Feuilles

| Feuille | Contenu |
|---|---|
| `_summary` | Comptages par filtre, hash de la configuration, date ; lien « ▶ Charger les BAMs » pour IGV |
| `HIGH_IMPACT` | Variants d'impact VEP HIGH |
| `CLINVAR_P_LP`, `CLINVAR_CONFLICTING`, `CLINVAR_VUS`, `CLINVAR_NOT_PROVIDED` | Par statut ClinVar |
| `FQ_MODERATE`, `FQ_LOW`, `FQ_MODIFIER` | Variants rares par impact VEP |
| `4_CALLERS` | Appelés par les quatre callers |
| `LOVD` | Présents dans LOVD |
| `MOSAIC` | AF entre 5 et 35 % |
| `DEEP_MOSAIC` | AF entre 0,5 et 5 %, profondeur ≥ 200×, biais de brin extrême exclu |
| gènes voisins (RHEB, PKD1) | Feuille fourre-tout, sans filtre qualité |

Les feuilles ne sont pas exclusives : un même variant peut apparaître dans
plusieurs d'entre elles. Chaque feuille est triée par `SCORE`, meilleur en tête.

Avant classement, des **pré-exclusions** retirent les variants trop fréquents
(gnomAD > 0,5 %, récurrence dans le run, PJCNT), homozygotes récurrents, ClinVar
bénins et artefacts de mosaïque récurrents. `--beta` les désactive pour vérifier
qu'aucun variant n'est exclu à tort.

## Colonne SCORE — couleur RGB sur trois axes

La cellule `SCORE` porte un score combiné (0–1) et une couleur qui encode trois
axes indépendants :

| Canal | Axe | Critères principaux |
|---|---|---|
| **R** rouge | Risque d'artefact (1 − qualité technique) | biais de brin, nombre de callers, DP / AF |
| **G** vert | Signal biologique | ClinVar, IMPACT, CADD, REVEL, SpliceAI, SPiP |
| **B** bleu | Rareté | récurrence dans le run, PJCNT, gnomAD, LOVD |

| Couleur | Lecture |
|---|---|
| Cyan (G + B) | Candidat idéal : biologie forte, rare, séquençage propre |
| Vert pur (G) | Biologie forte mais non rare : récurrent ou artefact de panel |
| Jaune (R + G) | Faux positif dangereux : biologie forte mais signal technique douteux |
| Rouge (R) | Artefact sans signal biologique |
| Noir | Sans intérêt, séquençage propre |

Le commentaire de la cellule détaille chaque critère et donne un verdict :
« Variant solide, à valider » (≥ 0,7), « À évaluer — signaux mixtes » (≥ 0,5),
« Probablement artefact ou inintéressant » (≥ 0,3), « Artefact probable ».

Quelques règles du score :

- **Paradoxe gnomAD** : absent de gnomAD mais PJCNT ≥ 10 → la contribution
  gnomAD est annulée. Un variant fréquent dans nos runs mais invisible en
  population est un artefact de panel, pas un variant rare.
- CADD ne pénalise pas les variants introniques profonds.
- Les prédicteurs d'épissage sont ignorés pour les frameshift et stop_gained.
- Une annotation absente (NaN) compte comme absente, pas comme un critère échoué.

## Signaux visuels sur HGVSc

| Signal | Signification |
|---|---|
| Police bleue | Le variant apparaît dans plusieurs feuilles (liste au survol) |
| Fond rouge | Discordance AF / run : HTZ récurrent ici en mosaïque (artefact probable), ou mosaïque récurrente ici en HTZ (**à vérifier sur IGV**) |
| Fond bleu (feuille gènes voisins) | `downstream_gene_variant` qui tombe dans PKD1 ou RHEB |
| Lien | HGVSc ouvre la fiche VarSome |

Les autres colonnes (SB, CALLNB, IMPACT, CLINVAR…) suivent un code classique :
vert favorable, jaune modéré, orange alerte, rouge défavorable.

## Colonnes à remplir

| Colonne | Usage |
|---|---|
| `USER_CLASS` | Classe ACMG (1–5) ou terme d'exclusion (`SB`, `ReMM low`, `Fréquence`) |
| `USER_ANNOT` | Annotation libre |
| `USER_COM` | Formule automatique qui reprend les annotations du même variant dans les autres feuilles |

Les termes d'exclusion (`igv_exclusion_terms`) retirent le variant des captures
IGV. KODEX n'écrase jamais un `*_filtered.xlsx` existant : une nouvelle sortie
prend le suffixe `_2`, `_3`…

## Colonne IGV

Avec IGV Desktop ouvert (port 60151), chaque ▶ navigue sur la position du
variant ; le lien de `_summary` charge une fois les BAM du patient et de
l'échantillon de référence. L'en-tête A1 (« ⇅ base ») trie l'alignement par base.
Les BAM sont cherchés dans `igv.bam_dir`, relatif au dossier des XLSX.

## Dépannage

| Problème | Solution |
|---|---|
| `micromamba: command not found` | Relancer `bash install.sh` ou `eval "$(micromamba shell hook -s bash)"` |
| `ModuleNotFoundError` | `micromamba activate kodex` |
| Colonne attendue manquante | Le fichier n'est pas du profil choisi : essayer `--legacy` |
| `#NOM?` dans les formules | Ouvert avec LibreOffice : il faut Microsoft Excel 2010+ |
| IGV ne répond pas | IGV Desktop ouvert, *Enable port* coché dans *View > Preferences > Advanced* |
| Rapport HTML de couverture vide | Voir [NOTICE_CNV.md](NOTICE_CNV.md) : le dossier `lib/` doit être à côté |

Le fichier d'entrée n'est jamais modifié, aucun accès réseau n'est fait (hors
IGV en localhost), et le hash de la configuration est écrit dans `_summary`
pour la traçabilité.
