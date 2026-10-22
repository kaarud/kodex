# KODEX — IGV "go-to position" : hyperlien Excel + CLI

**Date:** 2026-06-16
**Statut:** Spec validée — prêt pour plan d'implémentation

## Objectif

Permettre, depuis un `*_filtered.xlsx`, de sauter directement sur la position
d'un variant dans IGV Desktop en chargeant les BAMs du patient courant plus un
BAM de référence (« normal ») pour juger l'artefactness relative : un même
softclip présent sur le patient *et* sur la référence est probablement un
artefact.

Deux déclencheurs :

1. **Hyperlien Excel** (usage principal) — une colonne `IGV` par feuille de
   données ; un clic charge les BAMs et navigue.
2. **CLI `igv_goto.py`** (cœur / repli scriptable) — même comportement, invoqué
   à la main.

## Contexte

KODEX dispose déjà d'une intégration IGV (`igv_capture.py`,
`igv_bookmark.py`) qui parle au port IGV (défaut `60151`) via socket TCP brut et
HTTP GET. Fonctions réutilisées telles quelles :

- `check_igv(port)` — IGV joignable ?
- `igv_cmd(port, command, **params)` / `igv_batch(port, lines)` — envoi commandes
- `find_bams(sample_id, bam_dir)` — localise `{id}.bam`, `{id}.mutect2.bam`,
  `{id}.chim.bam` + index `.bai`
- `extract_sample_id(xlsx_path)` — `run_1_100008_S2_filtered.xlsx` → `100008_S2`
- `normalize_chr(chrom)` — préfixe `chr`
- `init_igv(port, bams, genome)` — `new` + `genome` + load des BAMs
- `_HALF_WINDOW = 145` — fenêtre ±145 pb (réutilisée pour centrer le `goto`)

Décisions de cadrage validées :

- **Emplacement des BAMs au clic** : copies **locales sur le Mac** (chemins de
  fichiers simples, comme le défaut `<run_dir>/bam/` de `igv_capture.py`).
- **Piste de référence** : **un seul échantillon fixe, BAM brut uniquement**
  (`.bam`, sans `.chim`/`.mutect2`).
- **Identification de la référence** : **par numéro de lane dans le run courant**
  (config `reference_lane: S2` → l'échantillon finissant par `_S2` dans le
  `bam_dir`). Survit aux nouveaux runs sans éditer l'id DEFGEN.
- **Approche déclencheur** : **C (CLI) comme cœur + A (hyperliens) par-dessus**.

## Architecture

Deux composants, tous deux appuyés sur `igv_capture.py` :

### 1. `igv_goto.py` (nouveau — cœur CLI)

Source de vérité du « que fait un clic ».

Invocations :

```bash
python igv_goto.py 100008_S2 chr16:2098100        # sample + position explicites
python igv_goto.py sample_filtered.xlsx --row 7    # résolution depuis une ligne XLSX
```

Comportement :

1. Résout les BAMs du patient courant via `find_bams(sample_id, bam_dir)`
   (RAW + CHIM + MUTECT2).
2. Résout la référence : scanne `bam_dir` pour l'échantillon finissant par la
   lane configurée (`_S2`) ; charge **uniquement son `.bam` brut** (label `REF_S2`).
   **Si le patient courant EST la référence, on saute la piste de référence**
   (pas de doublon).
3. Construit une fenêtre autour de POS (`half_window`, défaut 145) pour centrer
   l'affichage au lieu d'une base unique.
4. Envoie `new` (reset session) → load des 3 BAMs patient (labels RAW/CHIM/MUTECT2)
   + piste REF → `goto chr…:start-end`.
5. Erreurs explicites : IGV injoignable (`check_igv` faux), BAMs patient manquants
   (abandon), référence introuvable (charge le patient seul + avertit).

Arguments CLI : `sample_or_xlsx` (positionnel), `--row` (si xlsx), `--port`,
`--bam-dir`, `--genome`, `--config`, `--half-window`. Les défauts proviennent du
bloc `igv:` de la config.

### 2. Injection des hyperliens dans `vcf_filter.py`

Dans `write_output` (vcf_filter.py:1354), après écriture/formatage de chaque
feuille de données : ajout d'une colonne `IGV`.

- Le patient est constant par fichier (dérivé du nom via `extract_sample_id`),
  donc toutes les lignes d'une feuille partagent le même jeu de 4 BAMs ; seule la
  `locus` varie ligne à ligne.
- Texte de cellule : `▶` (ou `IGV`), avec `cell.hyperlink` pointant sur une URL
  GET unique IGV :

  ```
  http://localhost:60151/load?file=<cur.bam>,<cur.chim.bam>,<cur.mutect2.bam>,<ref.bam>
      &genome=hg38&merge=false&locus=chr16:2098000-2098290&name=RAW,CHIM,MUTECT2,REF_S2
  ```

- `merge=false` vide les pistes précédentes à chaque clic.
- Chemins = chemins absolus locaux, **URL-encodés**.
- Aucune dépendance runtime à `igv_goto.py` : l'URL est auto-suffisante.
- Si les BAMs du sample sont introuvables → colonne **silencieusement omise**
  (le reste du filtrage continue) + un avertissement console d'une ligne.

`write_output` reçoit les paramètres IGV nécessaires (port, genome, bam_dir,
reference_lane, half_window, sample_id) issus du bloc config + du nom de sortie.

## Construction de l'URL / des commandes (logique partagée)

Le cœur logique — résolution des BAMs patient + référence (avec dédup si
patient==référence), normalisation `chr`, fenêtre POS — est factorisé pour être
utilisé à la fois par `igv_goto.py` (commandes socket) et par l'injecteur de
hyperliens (URL GET). Cela garantit qu'un clic et la CLI font exactement la même
chose.

## Config (`config.yaml`)

```yaml
igv:
  enabled: true
  port: 60151
  genome: hg38
  bam_dir: <run_dir>/bam      # emplacement des BAMs locaux
  reference_lane: S2               # piste normale toujours chargée, matchée par suffixe _S2
  half_window: 145
```

## Gestion des erreurs

| Situation | CLI `igv_goto.py` | Hyperlien Excel |
|---|---|---|
| IGV non lancé | message clair + exit ≠ 0 | le clic échoue dans le navigateur (sans effet) |
| BAMs patient manquants | abandon + message | colonne omise + warning console à la génération |
| Référence introuvable | charge le patient seul + warning | URL sans la piste REF + warning console |
| Patient == référence | piste REF omise (pas de doublon) | piste REF omise de l'URL |

## Tests (pytest, sous `tests/`)

- **Constructeur d'URL** : liste de fichiers correcte, encodage, fenêtre `locus`
  (POS ± half_window), dédup de la référence quand patient == référence.
- **Résolution `reference_lane`** : depuis un `bam_dir` factice (fixtures), trouve
  le bon `_S2`, gère l'absence.
- **Injection de la colonne hyperlien** : sur un petit workbook openpyxl, vérifie
  présence colonne `IGV`, `cell.hyperlink` correct, omission propre si BAMs absents.
- Appels socket IGV **mockés** — pas d'IGV live dans les tests.

## Compromis connu

Cliquer un hyperlien `http://` dans Excel ouvre un onglet navigateur affichant la
réponse brute d'IGV à chaque clic (sans danger, mais s'accumule). Accepté pour
cette version. Évolution future possible (hors scope) : **Approche B** — schéma
d'URL custom `kodex-igv://…` capté par Hammerspoon (déjà utilisé dans chuchotte)
appelant `igv_goto.py`, supprimant l'onglet navigateur.

## Hors scope

- Approche B (Hammerspoon / schéma d'URL custom).
- Montage de partages serveur / streaming http des BAMs (les BAMs sont locaux).
- Pistes de référence multiples ou jeux de BAMs configurables au-delà de la lane.
