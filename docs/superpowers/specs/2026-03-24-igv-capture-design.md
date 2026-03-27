# Design — Script de capture IGV automatisée

**Date:** 2026-03-24
**Statut:** Validé

## Contexte

Après filtrage par `vcf_filter.py`, le biologiste annote les variants d'intérêt dans les colonnes `USER_CLASS`, `USER_ANNOT`, `USER_COM`. Le script `igv_capture.py` lit ce fichier filtré et génère automatiquement 2 captures IGV par variant annoté (strand bias + soft clips), avec les 3 BAMs empilés.

## Architecture

```
igv_capture.py <filtered_xlsx> [options]

  1. Lire TOUS les onglets du XLSX filtré (sauf _summary)
     → collecter variants avec USER_CLASS/ANNOT/COM non vide
     → dédupliquer par CHR+POS+REF+ALT
  2. Extraire sample_id depuis le nom de fichier
     (TSC_063_345822_S2_filtered.xlsx → 345822_S2)
  3. Localiser les 3 BAMs + leurs index (.bai) dans bam_dir :
     {id}.bam, {id}.mutect2.bam, {id}.chim.bam
     → vérifier existence BAM ET index pour chacun
     → si un BAM ou index manque : warning + skip de TOUS les variants du sample
  4. Vérifier que IGV répond sur le port (défaut 60151)
  5. Réinitialiser IGV : GET /new + GET /genome?name=hg38
  6. Charger les 3 BAMs (une seule fois)
  7. Pour chaque variant :
       a. Normaliser CHR (ajouter préfixe "chr" si absent)
       b. goto chr:pos-145 → pos+145  (fenêtre 290 pb centrée)
       c. sleep(delay)  ← délai APRÈS goto, avant snapshot
       d. Capture 1 : strand bias → PNG
       e. Capture 2 : soft clips → PNG (reset coloring après)
  8. Sauvegarder dans out_dir/{sample_id}/ (absolu, créé si absent)
```

## Commandes IGV (API HTTP port 60151)

### Initialisation (une fois par exécution)
```
GET /new
GET /genome?name=hg38   # genome IGV hébergé — si genome local : option --genome-path
GET /load?file=<abs_path>.bam&index=<abs_path>.bam.bai&name=RAW
GET /load?file=<abs_path>.mutect2.bam&index=<abs_path>.mutect2.bam.bai&name=MUTECT2
GET /load?file=<abs_path>.chim.bam&index=<abs_path>.chim.bam.bai&name=CHIM
```
**Convention index** : chercher `{file}.bai` en premier, fallback `{file}.bam.bai`.

### Par variant

La hauteur du panel est calculée dynamiquement depuis la colonne `DP` du XLSX :
```
panel_height = DP * 2 + 150   # mode squished ~2px/read + 150px coverage + headers
```
Exemples : DP=50 → 250px · DP=200 → 550px · DP=1000 → 2150px

```
# Adapter la hauteur au DP du variant (tous les reads visibles)
GET /setPreference?name=SAM.DISPLAY_MODE&value=SQUISHED
GET /setPreference?name=SAM.SAMPLING_COUNT&value=<DP>
GET /maxPanelHeight?value=<DP*2+150>

GET /goto?locus=chr9:132897468-132897758

sleep(delay)   # ← délai critique : IGV doit finir de charger les reads

# Capture 1 — strand bias
GET /setPreference?name=SAM.COLOR_BY&value=READ_STRAND
GET /snapshot?filename=<abs_out>/<SYMBOL>_<HGVSc_safe>_strand.png

# Capture 2 — soft clips
GET /setPreference?name=SAM.SHOW_SOFT_CLIPPED&value=true
GET /setPreference?name=SAM.COLOR_BY&value=NO_COLORING
GET /snapshot?filename=<abs_out>/<SYMBOL>_<HGVSc_safe>_softclip.png

# Reset pour le prochain variant
# IMPORTANT : strand capture TOUJOURS avant soft clips — ce reset est valide uniquement dans cet ordre
GET /setPreference?name=SAM.SHOW_SOFT_CLIPPED&value=false
```

## Nommage des fichiers de sortie

```
{out_dir}/{sample_id}/
  {SYMBOL}_{HGVSc_safe}_strand.png
  {SYMBOL}_{HGVSc_safe}_softclip.png
```

**`HGVSc_safe`** : caractères remplacés pour noms de fichiers valides :
`>→gt`, `/→_`, `*→star`, `:→_`, ` →_`

**Fallback si HGVSc vide** : utiliser `{CHR}_{POS}_{REF}_{ALT}` à la place.

**`out_dir` doit être résolu en chemin absolu** (`Path.expanduser().resolve()`) avant tout appel `/snapshot` — IGV écrit relatif à son propre répertoire de travail.

## CLI

```bash
igv_capture <filtered_xlsx> \
  [--port 60151] \
  [--bam-dir ~/data/TSC_063/bam/] \
  [--out ~/data/TSC_063/igv_captures/] \
  [--delay 2.0]
```

Defaults :
- `--port` : 60151
- `--bam-dir` : `~/data/TSC_063/bam/`
- `--out` : `~/data/TSC_063/igv_captures/`
- `--delay` : 2.0s (inséré après `/goto`, avant le premier snapshot — délai profond coverage)

## Gestion d'erreurs

| Cas | Comportement |
|-----|-------------|
| IGV non joignable (port 60151) | Message clair + exit(1) |
| BAM manquant (n'importe lequel des 3) | Warning niveau sample + skip de TOUS les variants du sample |
| CHR/POS invalide ou vide | Warning + skip du variant |
| HGVSc vide | Fallback nom fichier `CHR_POS_REF_ALT` |
| Aucun variant annoté | Message informatif + exit(0) |
| out_dir/{sample_id}/ existe déjà | Overwrite silencieux des PNGs existants |

## Prérequis utilisateur

- IGV Desktop ouvert manuellement avec GRCh38 chargé
- Le script envoie `/new` + `/genome?name=hg38` au démarrage → la session est réinitialisée
- Si le script est relancé à mi-run : les BAMs sont rechargés automatiquement

## Hors scope

- Lancement automatique d'IGV
- Reads chimériques (futur — paramètre SB 0.2–0.8 + IGV)
- Sélection interactive des variants via interface graphique
