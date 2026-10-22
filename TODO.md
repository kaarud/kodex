# TODO

## RPKM exon coverage curve / CNV

- [ ] **Bioconductor CNV calling (Approach C)** — add statistical exon-level CNV
  calling on top of the visual `rpkm_exon_curve.R`.
  - Use **raw read counts**, NOT the RPKM sheet (count models need library size +
    count magnitude; RPKM has those divided out).
  - Count source: `ExomeDepth::getBamCounts()`, `GATK CollectReadCounts`, or
    `featureCounts` over the panel BED, run on the BAMs (mosdepth `regions.bed.gz`
    is mean depth, close but not integer counts).
  - Candidate callers for a ~22-sample targeted TSC panel: **ExomeDepth**
    (beta-binomial + HMM, Bayes-factor calls) and/or **panelcn.MOPS**
    (cohort mixture-of-Poissons, panel-tuned); intersect both to raise specificity.
  - Overlay the resulting calls (with Bayes factor / quality) onto the existing
    coverage curve as an annotation layer.
- [ ] Optional: expose `--del` / `--dup` threshold flags in `rpkm_exon_curve.R`
  instead of the hardcoded 0.8 / 1.2 band.

## Per-base coverage curve (`perbase_coverage_curve.R`)

Spec: `docs/superpowers/specs/2026-06-17-perbase-coverage-curve-design.md`

- [x] **Per-sample-median normalization** — `--normalised` (default) divides each
  sample's per-base depth by that sample's own median over the plotted bases, so
  heights are comparable and real deletions drop below the cohort regardless of
  library size. `--raw` keeps absolute depth. (`normalise_by_sample_median`.)
- [x] **Cohort-reference normalisation (`--cohort`)** — second normalisation step
  after the per-sample median: divide each position by the cohort median there
  (DeCovA's *ratio2center*). Cancels shared per-base capture/GC bias so normal
  samples flatten to ~1.0 and CNVs deviate. (`normalise_by_cohort_reference`.)
- [ ] **Exhaustive DeCovA-like transformation (full `CNV_tool.pm` parity)** — the
  `--cohort` step above is only the first of DeCovA's pipeline. To reproduce its
  visual + clustering behaviour, port the rest of `lib/DeCovA/CNV_tool.pm`
  (ref: `github.com/HCL-HUBL/DeCovA`):
  - **Inter-sample library-size factor** — use summed depth across all regions
    (`Ref_depth = sample_total / mean_of_totals`) as the per-sample factor,
    instead of (or alongside) the per-sample median; closer to DeCovA's center.
  - **Region-level aggregation option** — DeCovA's CNV/clustering works on *mean
    depth per exon/BED interval*, not per base; add a region-mean mode to cut
    per-base noise before calling.
  - **Iterative reference refinement** — call CNVs, drop flagged regions/samples
    (`Appel`), recompute the normalisation factors, repeat until stable, so true
    CNVs don't contaminate the baseline.
  - **Region/sample exclusion + range** — port `--ex_region`, `--ex_sample`,
    `--ex_cov`, `--ex_DP`, and `--range` (keep only samples within
    median ± range·quartile for the center/spread calc).
  - **Center/spread CNV calls** — `level`/`spread` (avg|med, std|qtile) with
    `level_del` (0.8) / `level_dup` (1.2) and optional spread thresholds; emit
    per-region del/dup calls, not just a curve.
  - **Sex-awareness** — chrX/chrY handling (`--sex_file`, male chrX doubling,
    `ratioByGender`). N/A for the TSC1/TSC2 autosomal panel but needed for parity.
- [ ] **Strand-aware within-exon orientation** — v1 places bases within each exon
  by ascending genomic position, so minus-strand genes (e.g. TSC1) render mirrored
  vs mRNA 5'->3' (TSC2, plus-strand, is not). Exon-level order is already
  transcript-correct (by exon_num). To localize intra-exon features biologically,
  derive strand (from the gene or KoverGraph) and map `exonic_x = x_start +
  (gend - pos)` for minus-strand exons in `to_exonic_coords`. Harmless for the
  current cross-sample outlier QC view (shared x-mapping).


### Introns on the per-base curve — not started

Today `to_exonic_coords` keeps only positions that fall inside a panel exon and
concatenates the exons, so introns are dropped entirely. Adding them would show
the off-target / near-target coverage between exons, useful to spot a breakpoint
that sits outside the capture. Needs: intronic positions kept in
`assign_exons` (currently `nomatch = NULL` drops them), an x-axis that mixes
exonic and intronic stretches without letting a 30 kb intron swamp a 100 bp
exon (log-compressed intron width, or a fixed narrow width per intron), and a
decision on what the cohort reference means where depth is near zero.

### TSC2-PKD1 junction window — done (2026-09-17)

`genomic_layout()` lays exons on one axis by coordinate instead of by exon_num,
so PKD1 (minus strand, exon numbers running down as coordinates run up) is no
longer drawn backwards next to TSC2. `junction_marker()` draws the boundary and
names the 1116 bp gap it hides. Figures 4-6 of the default report cover
TSC2+PKD1 under the three renderings.

Note for whoever touches this next: `globalize_layout()` must NOT run on a
genomic layout — it re-offsets genes alphabetically and silently puts PKD1
before TSC2 (the gap then prints as -87074 bp). Guarded and tested.

Interpretation caveat carried on the figures themselves: PKD1 exons 1-33 sit in
the 16p segmental duplication with pseudogenes PKD1P1-P6 (>97% identity), so
short-read depth there is confounded; exons 34-46 are single-copy and border the
junction. The run bears this out — PKD1 E01 reaches a ratio of 5.5 in several
samples, which is a mapping artefact, not a duplication.


## Bugs

- [] fails to generate igv links when i change the reference bam, normally S2


---
User input
Quand je fais aggregate le script et que ça efface les anciennes remarques que j'ai fait sur le fichier ça me rend dingue à changer
il serait intéressant de savoir lesquels dans le runcount sont comptabilisés dans les commentaires
Il serait intéressant si l'information est disponible de mettre les données qualités de chaque patient, les concetration sur la page d'acceuil etc, à fusionner

Le plus important je pense ça serait que quand je fais un vcf_filter et qu'il voit qu'il y a déjà un tableur de créé de mettre un warning et de faire un backup de l'ancien ou de mettre un v2 v3 aux nouveaux
