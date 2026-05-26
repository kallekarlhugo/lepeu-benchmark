# LepEU — SNP Calling Benchmark in Lepidoptera

Pipeline and post-processing for benchmarking variant calling across 5 butterfly genomes,
comparing 2 aligners × 4 caller variants per species.

| Species | Aligners | Callers per species/aligner |
| --- | --- | --- |
| *Pieris napi* | bwa-mem2, NGM | HaplotypeCaller, bcftools (cohort, individual, cohort+HWE, cohort−HWE) |
| *Polyommatus icarus* | " | " |
| *Plebejus argus* | " | " |
| *Cyaniris semiargus* | " | " |
| *Lysandra bellargus* | " | " |

DeepVariant runs are defined in `runs.tsv` but are currently parked (known issue with grenepipe's per-(sample, contig) scatter on GPU; see [Known issues](#known-issues)).

---

## Pipeline overview

```
FASTQs ──► trim ──► map ──► sort ──► merge ──► dedup ──► BAMs ──► variant calling ──► VCF
                                                                                       │
                                                                                       ▼
                                                                                Post-processing
                                                                                       │
                                                              ┌────────────────────────┼────────────────────────┐
                                                              ▼                        ▼                        ▼
                                                          per-run filters       per-run stats              cross-run
                                                          + per-chrom           + pixy + vcftools FST      isec + multiqc
                                                          + annotation
```

Stage 1 (`launch-runs.sh`) does trim → map → dedup → call inside grenepipe.
Stage 2 (`postprocess.sh`) does everything downstream of the joint VCF.

---

## Step-by-step: what happens to each FASTQ

| Step | Tool | Input | Output |
| --- | --- | --- | --- |
| 1. **Trim** | fastp | `<sample>_R{1,2}.fastq.gz` | `runs/<run>/trimmed/<sample>_R{1,2}.fastq.gz` |
| 2. **Map** | bwa-mem2 OR NGM | trimmed FASTQs | `runs/<run>/mapping/sorted/<sample>.sorted.bam` |
| 3. **Merge** | samtools merge | sorted BAMs (if multiple lanes/units) | `runs/<run>/mapping/merged/<sample>.merged.bam` |
| 4. **Dedup** | Picard MarkDuplicates | merged BAM | `runs/<run>/mapping/dedup/<sample>.dedup.bam` |
| 5. **Final BAM** | symlinks | dedup BAM | `runs/<run>/mapping/final/<sample>.bam` (symlink → dedup) |
| 6. **Variant call** | bcftools / GATK HaplotypeCaller / DeepVariant | final BAMs + reference | per-(sample, contig) g.VCFs in `runs/<run>/calling/called/` |
| 7. **Joint genotype** | bcftools / GATK CombineGVCFs+GenotypeGVCFs / GLnexus | per-(sample, contig) g.VCFs | `runs/<run>/calling/genotyped-all.vcf.gz` |
| 8. **Normalize** | bcftools norm -m -any -f $REF --check-ref w | joint VCF | `postprocess/<run>/filtered/normalized.vcf.gz` |
| 9. **QUAL filter** | bcftools view -e 'QUAL<20' | normalized VCF | `postprocess/<run>/filtered/qual20.vcf.gz` |
| 10. **Biallelic SNPs** | bcftools view --types snps -m2 -M2 | qual20 VCF | `postprocess/<run>/filtered/qual20.biallelic.snps.vcf.gz` |
| 11. **Per-chromosome** | bcftools view -r `<chrom>` | qual20 and biallelic VCFs | `postprocess/<run>/per_chrom/<level>/<chrom>.vcf.gz` |
| 12. **VCF stats** | bcftools stats | filtered VCFs | `postprocess/<run>/stats/<level>.txt` |
| 13. **SFS** | bcftools query \| awk | biallelic VCF | `postprocess/<run>/sfs/{ac_histogram_*,summary}.tsv` |
| 14. **Pixy pi/dxy/FST** | pixy (Hudson FST) | normalized / qual20 / biallelic VCFs | `postprocess/<run>/pixy/<level>/{windows_10kb,windows_50kb,per_chrom,whole_genome}_*` |
| 15. **vcftools FST** | vcftools --weir-fst-pop | biallelic VCF + populations | `postprocess/<run>/vcftools_fst/{*_whole_genome,*_windows_10kb,*_windows_50kb,*_chr_*,fst_summary}.tsv` |
| 16. **Annotation** | bcftools csq -g GFF --force --phase a | biallelic VCF + species GFF | `postprocess/<run>/annotated/qual20.biallelic.snps.annotated.vcf.gz` (with BCSQ + IMPACT INFO) |
| 17. **Cross-run intersect** | bcftools isec --collapse none | all biallelic VCFs of a species | `postprocess/isec/<species>/{sites.txt, run_index.tsv, per_run_counts.tsv, per_caller_counts.tsv, per_aligner_counts.tsv, concordance_histogram.tsv, pattern_counts.tsv}` |
| 18. **Aggregate report** | MultiQC | all bcftools-stats outputs | `postprocess/multiqc/lepeu-postprocess.html` |

---

## Data layout

```
/lunarc/nobackup/projects/lepeu-lisbon/
├── data/                                          raw FASTQ reads
│
├── kalle_temp/LepEU/                              project workspace
│   ├── refs/                                      reference genomes (.fna + .fai + .dict + bwa indices)
│   │   └── lepeu_annotations/                     per-species AnnEvo GFFs
│   ├── <species>_samples.tsv                      one row per sample: sample, unit, population, fq1, fq2
│   ├── runs.tsv                                   the run matrix
│   ├── runs/<run_id>/                             Stage 1 outputs (per pipeline run)
│   ├── postprocess/<run_id>/                      Stage 2 outputs (per pipeline run)
│   ├── postprocess/isec/<species>/                cross-run intersection
│   ├── postprocess/multiqc/                       aggregate HTML report
│   ├── archive/                                   stashed intermediate state from prior attempts
│   ├── grenepipe/                                 modded grenepipe workflow (moved here from $HOME for shared access)
│   ├── launch-runs.sh                             Stage 1 launcher
│   ├── postprocess.sh                             Stage 2 launcher
│   └── build_status.sh                            per-run progress tracker → status.tsv
│
└── shared/conda/
    ├── grenepipe/                                 snakemake + grenepipe deps (read-only)
    └── lepeu-pixy/                                pixy + bcftools + vcftools + multiqc + tabix
```

### Reference accession → species

| Species | Reference accession |
| --- | --- |
| *Pieris napi* | GCA_905475465.2 |
| *Polyommatus icarus* | GCA_937595015.1 |
| *Plebejus argus* | GCA_905404155.3 |
| *Cyaniris semiargus* | GCA_905187585.1 |
| *Lysandra bellargus* | GCA_905333045.1 |

GFFs auto-match by `GCA_<accession>` pattern in `refs/lepeu_annotations/`.

### `<species>_samples.tsv` schema

```
sample        unit  population  fq1                       fq2
Pnapi_A1      1     A           /lunarc/.../Pnapi_A1_R1.fastq.gz  /lunarc/.../Pnapi_A1_R2.fastq.gz
...
```

### `runs.tsv` schema

```
run_id  reference  samples_table  mapping_tool  calling_tool  calling_mode  bams_from
```

- **anchors** have `bams_from = -` and run trim → map → dedup → call from scratch
- **dependents** have `bams_from = <anchor_run_id>` and reuse the anchor's BAMs via a generated `mappings.tsv`

`calling_mode` is `default` for non-bcftools rows, otherwise `combined`, `individual`, `combined_hwe`, or `combined_nohwe`.

---

## How to read the outputs

| What you want | File |
| --- | --- |
| Joint VCF for a run | `runs/<run>/calling/genotyped-all.vcf.gz` |
| QUAL≥20 filtered VCF (keeps invariants) | `postprocess/<run>/filtered/qual20.vcf.gz` |
| Biallelic SNPs only (post-norm, indels dropped) | `postprocess/<run>/filtered/qual20.biallelic.snps.vcf.gz` |
| Functionally annotated VCF (BCSQ + IMPACT) | `postprocess/<run>/annotated/qual20.biallelic.snps.annotated.vcf.gz` |
| Whole-genome pi/dxy | `postprocess/<run>/pixy/qual20/whole_genome_{pi,dxy}.tsv` |
| Whole-genome FST (canonical, W&C ratio-of-averages) | `postprocess/<run>/vcftools_fst/fst_summary.tsv` row `region=whole_genome`, column `weighted_fst` |
| Whole-genome FST (pixy Hudson, approximate) | `postprocess/<run>/pixy/qual20/whole_genome_fst.tsv` |
| Per-chromosome pi/dxy/FST | `postprocess/<run>/pixy/qual20/per_chrom_{pi,dxy,fst}.txt` |
| 10 kb windowed pi/dxy/FST | `postprocess/<run>/pixy/qual20/windows_10kb_{pi,dxy,fst}.txt` and `postprocess/<run>/vcftools_fst/<p1>_vs_<p2>_windows_10kb.windowed.weir.fst` |
| 50 kb windowed pi/dxy/FST | `postprocess/<run>/pixy/qual20/windows_50kb_*` and `postprocess/<run>/vcftools_fst/<p1>_vs_<p2>_windows_50kb.windowed.weir.fst` |
| Allele-count histogram / SFS | `postprocess/<run>/sfs/ac_histogram_{folded,unfolded}.tsv` |
| Singleton/doubleton counts | `postprocess/<run>/sfs/summary.tsv` |
| Consequence-bucket counts (intron, missense, …) | `postprocess/<run>/annotated/csq_consequence_counts.tsv` |
| HIGH/MODERATE/LOW/MODIFIER counts | `postprocess/<run>/annotated/csq_impact_counts.tsv` |
| Cross-caller concordance (per species) | `postprocess/isec/<species>/concordance_histogram.tsv` |
| Variants seen by each caller | `postprocess/isec/<species>/per_caller_counts.tsv` |
| Variants seen by each aligner | `postprocess/isec/<species>/per_aligner_counts.tsv` |
| Aggregate quality report | `postprocess/multiqc/lepeu-postprocess.html` |

### FST: which estimator?

- **pixy** uses **Hudson** (1992), set explicitly via `--fst_type hudson`. Per-window and per-chrom.
- **vcftools** uses **Weir & Cockerham** (1984), the only estimator `--weir-fst-pop` supports. Per-window, per-chrom, and whole-genome.

For whole-genome, **use vcftools' `weighted_fst`** — it's the ratio-of-averages aggregation (Bhatia et al. 2013), theoretically correct. The pixy `whole_genome_fst.tsv` is a weighted-mean approximation kept for sanity checking; the per-window and per-chrom pixy values are exact and aggregate cleanly.

---

## Running it

### Setup (one-time)

```bash
# Pixy env (writable)
mamba create -y -p /lunarc/.../shared/conda/lepeu-pixy \
    -c bioconda -c conda-forge \
    pixy bcftools vcftools multiqc tabix
```

### Stage 1 — variant calling

```bash
# Anchors first (one job per species × aligner)
./launch-runs.sh --anchors-only

# Wait for them to finish (each ~24–96 h on lu48), then dependents
./launch-runs.sh --dependents-only
```

Each row of `runs.tsv` becomes one sbatch job. The launcher copies grenepipe's config, overrides per-run fields, generates `populations.tsv` (for `combined_hwe`) and `mappings.tsv` (for dependents), and runs snakemake with `--executor local`.

### Stage 2 — post-processing

```bash
./postprocess.sh --list                  # status overview
./postprocess.sh --submit                # one sbatch per finished Stage-1 run
./postprocess.sh --isec [species ...]    # cross-run intersection, local
./postprocess.sh --multiqc               # aggregate HTML report, local
```

Each Stage 2 sbatch runs the 17-step pipeline above in one process. Idempotent — outputs that already exist are skipped on re-runs.

### Progress tracking

```bash
./build_status.sh                        # writes status.tsv with Y/N per step
```

Open `status.tsv` in a spreadsheet or `column -t -s$'\t' status.tsv | less -S`.

---

## Modifications to grenepipe

Local fork at `kalle_temp/LepEU/grenepipe/` differs from upstream <https://github.com/lczech/grenepipe>:

- Added NGM aligner support (`workflow/rules/mapping-ngm.smk`)
- Added DeepVariant + GLnexus calling pipeline (`workflow/rules/calling-deepvariant.smk`)
- `keep-invariant-sites: "true"` flows across all callers
- `multiqc.yaml`: added `snakemake-wrapper-utils` dependency
- `profiles/lunarc-slurm/`: `singularity-args` / `apptainer-args` set with `--nv` and Lunarc bind mounts

---

## Known issues

- **DeepVariant** is currently parked. Grenepipe's per-(sample, contig) scatter combined with CUDA init overhead is slow; will be revisited.
- **Pixy 2.0.0 numpy bug** crashes on 10 kb windows where one population has zero called samples in a window. Causes `pixy/unfiltered/` and `pixy/qual20.biallelic.snps/` 10 kb outputs to be missing for some runs. Per-chrom and qual20 10 kb / 50 kb outputs are unaffected.
- **Pixy on biallelic-SNPs-only VCFs gives biased pi/dxy** (invariant sites missing from denominator). Use `pixy/qual20/` for canonical pi/dxy.
- **bcftools csq** uses the AnnEvo GFFs, which don't cover all unplaced scaffolds. `--force` lets csq skip those gracefully (variants on unannotated contigs get no BCSQ entry; IMPACT=MODIFIER).
