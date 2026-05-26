#!/usr/bin/env bash
#
# postprocess.sh — post-processing for the LepEU benchmarking matrix.
#
# For each finished run (calling/genotyped-all.vcf.gz.done present):
#   0. norm + split multi-allelics on the raw VCF  -> filtered/normalized.vcf.gz
#      (essential for fair cross-caller comparison: HC and DV don't auto-normalize,
#      bcftools-called VCFs have a no-split norm. All downstream steps operate on
#      the normalized output.)
#   1. Filter on QUAL>=20                          -> filtered/qual20.vcf.gz
#      (invariant sites kept; QUAL='.' doesn't match <20)
#   2. Keep only biallelic SNPs (indels dropped)   -> filtered/qual20.biallelic.snps.vcf.gz
#      (norm already split multi-allelics in Step 0, so this is a simple type filter)
#   3. Extract per-chromosome subsets of BOTH filtered VCFs                -> per_chrom/...
#   4. bcftools stats on whole-genome filtered VCFs                        -> stats/*.txt
#   5. AC histograms (folded + unfolded SFS) + singleton/doubleton counts  -> sfs/*.tsv
#   6. pixy pi/fst/dxy on three filter levels: unfiltered / qual20 / qual20.biallelic.snps
#        - sliding windows at each size in $PIXY_WINDOW_SIZES (default: 10kb + 50kb)
#        - per-chromosome (full-length BED entries)
#        - whole-genome pi/dxy (exact, aggregated from per-chrom by summing
#          count_diffs/count_comparisons)
#        - whole-genome FST (APPROXIMATE, weighted mean by no_snps) —
#          for sanity-checking against vcftools' proper ratio-of-averages
#          (Step 7). Use vcftools' weighted_fst as the canonical value.
#        - FST estimator: Hudson (--fst_type hudson), not W&C
#        - WARNING: pixy on biallelic-SNPs-only VCFs gives BIASED pi/dxy because
#          invariant sites are missing from the denominator. Use the qual20
#          numbers as canonical; biallelic results are useful only for relative
#          comparison across runs.
#   7. vcftools Weir & Cockerham FST (proper ratio-of-averages "weighted" FST):
#        - whole-genome, sliding windows ($PIXY_WINDOW_SIZES — 10kb + 50kb default), per-chromosome
#        - reported in vcftools_fst/fst_summary.tsv. Two FST estimators in play now:
#          pixy = Hudson, vcftools = W&C — both useful, different assumptions.
#   8. bcftools csq annotation on the biallelic-SNPs VCF using the species-matched
#      GFF (parsed by GCA accession from the ref filename)
#        -> annotated/qual20.biallelic.snps.annotated.vcf.gz with TWO INFO fields:
#             BCSQ   — consequence terms from bcftools csq (missense, stop_gained, ...)
#             IMPACT — derived SnpEff/VEP-style bucket (HIGH/MODERATE/LOW/MODIFIER),
#                      the most severe consequence per record. Missing-consequence
#                      records (no transcript overlap) get IMPACT=MODIFIER.
#        -> annotated/csq_consequence_counts.tsv (counts per consequence)
#        -> annotated/csq_impact_counts.tsv     (counts per HIGH/MODERATE/LOW/MODIFIER)
#        -> annotated/impact.tsv.gz             (CHROM/POS/REF/ALT/IMPACT lookup,
#                                                tabix-indexed; intermediate kept
#                                                for downstream use)
#        Filter on impact downstream with e.g.:
#             bcftools view -i 'INFO/IMPACT="HIGH"' annotated.vcf.gz
#
# Aggregate everything with multiqc into one report.
#
# Usage:
#   ./postprocess.sh --list                   # show finished/pending status
#   ./postprocess.sh --submit                 # submit one SLURM job per finished run
#   ./postprocess.sh --submit run1 run2 ...   # only those runs
#   ./postprocess.sh --run <run_id>           # invoked inside sbatch; processes one run
#   ./postprocess.sh --multiqc                # aggregate all stats files into one report
#   ./postprocess.sh --isec                   # cross-run bcftools isec per species
#                                             # (concordance + per-caller / per-aligner unique counts)
#   ./postprocess.sh --isec Picarus Pnapi     # only those species
#
# Env:
#   DRY_RUN=1     With --submit: print sbatch commands without running them.
#   SKIP_PIXY=1   With --run:    skip the pixy section (e.g. for a stats-only re-run).
#   PIXY_ENV=...  Override conda env used for pixy (default: same as CONDA_ENV).

set -euo pipefail

# =====================================================================
#     Configuration
# =====================================================================
RUNS_ROOT="$(pwd -P)/runs"
POST_ROOT="$(pwd -P)/postprocess"

SLURM_ACCOUNT="lu2026-2-62"
PP_CORES=8
PP_MEM="32G"
PP_TIME="08:00:00"

# Sliding-window sizes for pixy + vcftools (space-separated, in bp).
# Each becomes one set of windowed outputs.  Whole-genome and per-chromosome
# aggregations are always emitted in addition to these.
PIXY_WINDOW_SIZES="10000 50000"

CONDA_BASE="$HOME/miniconda3"
# One env with all post-processing tools (bcftools, vcftools, multiqc, tabix, pixy).
# Install with:
#   mamba install -p $CONDA_ENV -c bioconda -c conda-forge bcftools vcftools multiqc tabix pixy
CONDA_ENV="${POSTPROC_ENV:-/lunarc/nobackup/projects/lepeu-lisbon/shared/conda/lepeu-pixy}"
PIXY_ENV="$CONDA_ENV"                                                            # pixy is in the same env
VCFTOOLS_MODULE="${VCFTOOLS_MODULE:-}"                                           # leave blank — vcftools is in CONDA_ENV; set to module names only if you need to fall back
# =====================================================================

usage() {
    sed -n '3,28p' "$0" >&2
}

done_marker() { echo "$1/calling/genotyped-all.vcf.gz.done"; }
input_vcf()   { echo "$1/calling/genotyped-all.vcf.gz"; }
is_finished() { [[ -f "$(done_marker "$1")" && -f "$(input_vcf "$1")" ]]; }

list_finished_runs() {
    for d in "$RUNS_ROOT"/*/; do
        [[ -d "$d" ]] || continue
        is_finished "$d" && basename "$d"
    done
}

# Parse a quoted YAML scalar at `^  <key>:` from a config file.
yaml_value() {
    local file="$1" key="$2"
    grep -E "^  ${key}:" "$file" | head -1 | \
        sed -E 's/^  [^:]+: *"?([^"]*)"?[[:space:]]*$/\1/'
}

# =====================================================================
#     Mode dispatch
# =====================================================================
cmd="${1:-}"
shift || true

case "$cmd" in
  --list)
    printf '%-6s  %s\n' STATUS RUN
    for d in "$RUNS_ROOT"/*/; do
        [[ -d "$d" ]] || continue
        run=$(basename "$d")
        if is_finished "$d"; then
            if [[ -f "$POST_ROOT/$run/.done" ]]; then
                printf '%-6s  %s\n' DONE "$run"
            else
                printf '%-6s  %s\n' READY "$run"
            fi
        else
            printf '%-6s  %s\n' WAIT "$run"
        fi
    done
    exit 0
    ;;

  --submit)
    runs=("$@")
    if [[ ${#runs[@]} -eq 0 ]]; then
        mapfile -t runs < <(list_finished_runs)
    fi
    SCRIPT_PATH="$(readlink -f "$0")"
    SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
    submitted=0; skipped=0
    for run in "${runs[@]}"; do
        run_dir="$RUNS_ROOT/$run"
        if ! is_finished "$run_dir"; then
            echo "SKIP $run (calling not done)"
            skipped=$((skipped+1)); continue
        fi
        pp_dir="$POST_ROOT/$run"
        if [[ -f "$pp_dir/.done" ]]; then
            echo "SKIP $run (already post-processed)"
            skipped=$((skipped+1)); continue
        fi
        mkdir -p "$pp_dir"

        payload="set -euo pipefail
source '$CONDA_BASE/etc/profile.d/conda.sh'
conda activate '$CONDA_ENV'
cd '$SCRIPT_DIR'
'$SCRIPT_PATH' --run '$run'"

        sbatch_args=(
            --account="$SLURM_ACCOUNT"
            --time="$PP_TIME"
            --cpus-per-task="$PP_CORES"
            --mem="$PP_MEM"
            --job-name="pp-$run"
            --output="$pp_dir/slurm-%j.out"
            --error="$pp_dir/slurm-%j.err"
            --wrap="$payload"
        )
        if [[ "${DRY_RUN:-0}" == "1" ]]; then
            echo "DRY_RUN: pp-$run"
            echo "  sbatch ${sbatch_args[*]}"
        else
            echo "Submitting pp-$run"
            sbatch "${sbatch_args[@]}"
        fi
        submitted=$((submitted+1))
    done
    echo
    echo "Submitted: $submitted   Skipped: $skipped"
    exit 0
    ;;

  --multiqc)
    OUT="$POST_ROOT/multiqc"
    mkdir -p "$OUT"
    if ! command -v multiqc >/dev/null 2>&1; then
        echo "ERROR: multiqc not on PATH. Activate the grenepipe env first:" >&2
        echo "  source $CONDA_BASE/etc/profile.d/conda.sh && conda activate $CONDA_ENV" >&2
        exit 1
    fi
    # MultiQC auto-detects bcftools_stats files anywhere under POST_ROOT
    multiqc -f -o "$OUT" -n "lepeu-postprocess" \
        --module bcftools \
        "$POST_ROOT"
    echo "Report: $OUT/lepeu-postprocess.html"
    exit 0
    ;;

  --isec)
    # Cross-run analysis: bcftools isec across all post-processed runs of each species.
    # Input is each run's qual20.biallelic.snps.vcf.gz (normalized + split + biallelic SNPs).
    #
    # Usage:
    #   ./postprocess.sh --isec                  # all species detected from $POST_ROOT
    #   ./postprocess.sh --isec Picarus Pnapi    # only those species
    species_args=("$@")
    if [[ ${#species_args[@]} -eq 0 ]]; then
        mapfile -t species_args < <(
            ls -d "$POST_ROOT"/*/ 2>/dev/null \
              | awk -F'/' '{print $(NF-1)}' \
              | grep -Ev '^(isec|multiqc)$' \
              | awk -F'_' '{print $1}' \
              | sort -u
        )
    fi
    if [[ ${#species_args[@]} -eq 0 ]]; then
        echo "ERROR: no post-processed run dirs found under $POST_ROOT" >&2
        exit 1
    fi

    if ! command -v bcftools >/dev/null 2>&1; then
        echo "ERROR: bcftools not on PATH" >&2
        exit 1
    fi

    ISEC_ROOT="$POST_ROOT/isec"
    mkdir -p "$ISEC_ROOT"

    for species in "${species_args[@]}"; do
        echo
        echo "=== ISEC: $species ==="
        isec_dir="$ISEC_ROOT/$species"

        # Collect biallelic-SNPs VCFs for this species
        vcfs=()
        labels=()
        callers=()
        aligners=()
        for d in "$POST_ROOT/${species}"_*/; do
            [[ -d "$d" ]] || continue
            run=$(basename "$d")
            vcf="$d/filtered/qual20.biallelic.snps.vcf.gz"
            if [[ ! -f "$vcf" ]]; then
                echo "  skip $run (no biallelic SNPs VCF — has the post-process finished?)"
                continue
            fi
            # Parse: <species>_<aligner>_<caller...possibly_with_underscores>
            rest=${run#*_}
            aligner=${rest%%_*}
            caller=${rest#*_}
            vcfs+=("$vcf")
            labels+=("$run")
            callers+=("$caller")
            aligners+=("$aligner")
        done

        n=${#vcfs[@]}
        if (( n < 2 )); then
            echo "  WARN: only $n post-processed VCF(s) for $species — need >=2 for isec; skipping"
            continue
        fi

        echo "  $n VCFs:"
        for i in "${!labels[@]}"; do
            printf "    %2d  %-50s  aligner=%-10s caller=%s\n" "$i" "${labels[i]}" "${aligners[i]}" "${callers[i]}"
        done

        # Fresh prefix dir
        rm -rf "$isec_dir"
        echo "  running bcftools isec..."
        bcftools isec \
            --nfiles +1 \
            --collapse none \
            --output-type z \
            --prefix "$isec_dir" \
            "${vcfs[@]}" \
            2> "$isec_dir.isec.log" || {
                echo "  ERROR: bcftools isec failed for $species (see ${isec_dir}.isec.log)" >&2
                continue
            }
        mv "$isec_dir.isec.log" "$isec_dir/bcftools_isec.log" 2>/dev/null || true

        # Write the input -> run_id/aligner/caller mapping (0-based index = bitmask position L->R)
        {
            printf 'index\trun_id\taligner\tcaller\n'
            for i in "${!labels[@]}"; do
                printf '%d\t%s\t%s\t%s\n' "$i" "${labels[i]}" "${aligners[i]}" "${callers[i]}"
            done
        } > "$isec_dir/run_index.tsv"

        # One awk pass over sites.txt + run_index.tsv to produce all summaries
        awk -F'\t' \
            -v out_run="$isec_dir/per_run_counts.tsv" \
            -v out_conc="$isec_dir/concordance_histogram.tsv" \
            -v out_caller="$isec_dir/per_caller_counts.tsv" \
            -v out_aligner="$isec_dir/per_aligner_counts.tsv" \
            -v out_pat="$isec_dir/pattern_counts.tsv" '
        # Phase 1: run_index.tsv
        NR==FNR {
            if (FNR == 1) next
            idx = $1 + 0
            run[idx] = $2; aligner[idx] = $3; caller[idx] = $4
            n_runs = idx + 1
            if (!seen_caller[$4]++)   caller_list[++n_callers] = $4
            if (!seen_aligner[$3]++)  aligner_list[++n_aligners] = $3
            next
        }
        # Phase 2: sites.txt — columns: CHROM POS REF ALT bitmask
        {
            mask = $5
            pattern_count[mask]++

            ones = 0
            for (j=1; j<=length(mask); j++) ones += (substr(mask,j,1)=="1")
            concordance[ones]++

            # Per-run: count variants where this input has the variant, and (separately) where ONLY this input has it.
            for (i=0; i<n_runs; i++) {
                if (substr(mask, i+1, 1) == "1") {
                    in_set[i]++
                    if (ones == 1) unique[i]++
                }
            }

            # Per-caller: variant is present in any input belonging to this caller
            for (k=1; k<=n_callers; k++) {
                c = caller_list[k]; seen=0
                for (i=0; i<n_runs; i++) {
                    if (caller[i] == c && substr(mask, i+1, 1) == "1") { seen=1; break }
                }
                if (seen) caller_count[c]++
            }

            # Per-aligner: same idea
            for (k=1; k<=n_aligners; k++) {
                a = aligner_list[k]; seen=0
                for (i=0; i<n_runs; i++) {
                    if (aligner[i] == a && substr(mask, i+1, 1) == "1") { seen=1; break }
                }
                if (seen) aligner_count[a]++
            }
        }
        END {
            print "index\trun_id\taligner\tcaller\tin_set\tunique_to_run" > out_run
            for (i=0; i<n_runs; i++) {
                printf "%d\t%s\t%s\t%s\t%d\t%d\n", i, run[i], aligner[i], caller[i], in_set[i]+0, unique[i]+0 > out_run
            }

            print "n_runs_with_variant\tvariant_count" > out_conc
            for (k=1; k<=n_runs; k++) print k"\t"concordance[k]+0 > out_conc

            print "caller\tvariant_count" > out_caller
            for (k=1; k<=n_callers; k++) print caller_list[k]"\t"caller_count[caller_list[k]]+0 > out_caller

            print "aligner\tvariant_count" > out_aligner
            for (k=1; k<=n_aligners; k++) print aligner_list[k]"\t"aligner_count[aligner_list[k]]+0 > out_aligner

            print "pattern\tvariant_count" > out_pat
            for (p in pattern_count) print p"\t"pattern_count[p] > out_pat
        }' "$isec_dir/run_index.tsv" "$isec_dir/sites.txt"

        # Sort pattern_counts descending by count for human consumption
        ( head -1 "$isec_dir/pattern_counts.tsv";
          tail -n +2 "$isec_dir/pattern_counts.tsv" | sort -k2 -rn
        ) > "$isec_dir/pattern_counts.sorted.tsv" && \
            mv "$isec_dir/pattern_counts.sorted.tsv" "$isec_dir/pattern_counts.tsv"

        total_variants=$(wc -l < "$isec_dir/sites.txt")
        all_n=$(awk -v n=$n 'BEGIN{s=""; for (i=0;i<n;i++) s=s"1"; print s}')
        intersection=$(awk -v m="$all_n" '$5==m' "$isec_dir/sites.txt" | wc -l)
        echo "  -> $total_variants variants in union; $intersection in ALL inputs ($(awk "BEGIN{printf \"%.1f\", 100*$intersection/$total_variants}")%)"
        echo "  -> $isec_dir/{per_run_counts,per_caller_counts,per_aligner_counts,concordance_histogram,pattern_counts}.tsv"
    done
    exit 0
    ;;

  --run)
    RUN="${1:-}"
    [[ -z "$RUN" ]] && { echo "ERROR: --run requires a run_id" >&2; exit 1; }
    # fall through to the work section below
    ;;

  -h|--help|"")
    usage; exit 0 ;;

  *)
    echo "ERROR: unknown command '$cmd'" >&2
    usage; exit 1 ;;
esac

# =====================================================================
#     Per-run work (--run)
# =====================================================================
RUN_DIR="$RUNS_ROOT/$RUN"
IN_VCF="$RUN_DIR/calling/genotyped-all.vcf.gz"
OUT_DIR="$POST_ROOT/$RUN"

is_finished "$RUN_DIR" || { echo "ERROR: $RUN calling not done" >&2; exit 1; }

mkdir -p "$OUT_DIR"/{filtered,stats,sfs,per_chrom,pixy,logs}

# --- read reference + samples table from the per-run config snapshot ---
REF=$(yaml_value "$RUN_DIR/config.yaml" reference-genome)
SAMPLES_TSV=$(yaml_value "$RUN_DIR/config.yaml" samples-table)
FAI="${REF}.fai"
[[ -f "$REF"         ]] || { echo "ERROR: reference $REF not found"     >&2; exit 1; }
[[ -f "$FAI"         ]] || { echo "ERROR: index $FAI not found"          >&2; exit 1; }
[[ -f "$SAMPLES_TSV" ]] || { echo "ERROR: samples $SAMPLES_TSV not found" >&2; exit 1; }

# --- pixy populations.tsv (sample <TAB> population from samples_table col 3) ---
# pixy expects NO header — first row is treated as a sample. Likewise vcftools'
# --weir-fst-pop reads one sample per line, no header.
POP_TSV="$OUT_DIR/populations.tsv"
awk -F'\t' 'NR>1 && $1!="" && $3!="" && !seen[$1]++ {print $1"\t"$3}' \
    "$SAMPLES_TSV" > "$POP_TSV"
if [[ ! -s "$POP_TSV" ]]; then
    echo "WARN [$RUN]: populations.tsv empty (no population column in $SAMPLES_TSV?) — pixy will be skipped" >&2
    SKIP_PIXY=1
fi

echo "[$RUN] start $(date -Iseconds)"
echo "       input   : $IN_VCF"
echo "       ref     : $REF"
echo "       samples : $SAMPLES_TSV"
echo "       pops    : $(wc -l < "$POP_TSV") entries"

# =====================================================================
#     Step 0+1+2: normalize, then whole-genome filters
# =====================================================================
# All downstream steps (filters, per-chrom, pixy, vcftools-fst, isec, annotation)
# operate on a NORMALIZED VCF. This is essential for fair cross-caller comparison:
# different callers (HC, DV, bcftools) can emit the same variant with different
# REF/ALT representations (e.g. one row with multi-ALT vs. multiple rows), and
# pixy/isec would treat those as different sites otherwise.
#
# `norm --multiallelics -any` splits multi-allelic rows into one row per ALT.
# `-f REF --check-ref w` left-aligns indels against the reference and warns
# about REF mismatches (doesn't fail). Grenepipe's bcftools caller already ran
# a no-split norm; this step is a near no-op there but essential for HC and DV
# (which don't auto-normalize).
NORM_VCF="$OUT_DIR/filtered/normalized.vcf.gz"
QUAL_VCF="$OUT_DIR/filtered/qual20.vcf.gz"
BIALL_VCF="$OUT_DIR/filtered/qual20.biallelic.snps.vcf.gz"
IOTH=$(( PP_CORES > 2 ? PP_CORES - 2 : 1 ))

# Step 0: normalize the raw VCF — invariant + variant + indels all kept
if [[ ! -f "$NORM_VCF" ]]; then
    echo "[$RUN] norm + split multi-allelics on raw VCF..."
    bcftools norm \
        --multiallelics -any \
        --fasta-ref "$REF" \
        --check-ref w \
        --threads "$IOTH" \
        -O z -o "$NORM_VCF" \
        "$IN_VCF" \
        2> "$OUT_DIR/logs/norm.log"
    bcftools index --tbi --threads 2 "$NORM_VCF"
fi

# Step 1: QUAL>=20 — invariant sites (QUAL=.) are NOT excluded because '.' is not <20
if [[ ! -f "$QUAL_VCF" ]]; then
    echo "[$RUN] filter QUAL>=20..."
    bcftools view --threads "$IOTH" -e 'QUAL<20' -Oz -o "$QUAL_VCF" "$NORM_VCF"
    bcftools index --tbi --threads 2 "$QUAL_VCF"
fi

# Step 2: biallelic SNPs only — norm has already split multi-allelics, so a
# simple `view --types snps -m2 -M2` cleanly drops indels and keeps only the
# SNP rows. Records like REF=A ALT=AT become single-allele AT-insertion rows
# that --types snps drops.
if [[ ! -f "$BIALL_VCF" ]]; then
    echo "[$RUN] filter biallelic SNPs (indels excluded)..."
    bcftools view \
        --types snps \
        --min-alleles 2 --max-alleles 2 \
        --threads "$IOTH" \
        -O z -o "$BIALL_VCF" \
        "$QUAL_VCF"
    bcftools index --tbi --threads 2 "$BIALL_VCF"
fi

# =====================================================================
#     Step 3: per-chromosome subsets of both filtered VCFs
# =====================================================================
mapfile -t CHROMS < <(awk '{print $1}' "$FAI")
echo "[$RUN] per-chrom subsets across ${#CHROMS[@]} chromosomes..."

for label in qual20 qual20.biallelic.snps; do
    src="$OUT_DIR/filtered/${label}.vcf.gz"
    pcdir="$OUT_DIR/per_chrom/${label}"
    mkdir -p "$pcdir"
    for chrom in "${CHROMS[@]}"; do
        # Sanitize chrom name for filename (replace / and : just in case)
        safe=${chrom//\//_}; safe=${safe//:/_}
        out="$pcdir/${safe}.vcf.gz"
        if [[ ! -f "$out" ]]; then
            bcftools view --threads 2 -r "$chrom" -Oz -o "$out" "$src" 2>/dev/null || {
                # No records for this chrom is OK; clean up the empty output
                rm -f "$out"
            }
            [[ -f "$out" ]] && bcftools index --tbi "$out" 2>/dev/null || true
        fi
    done
done

# =====================================================================
#     Step 4: bcftools stats on whole-genome filtered VCFs
# =====================================================================
echo "[$RUN] bcftools stats..."
for label in qual20 qual20.biallelic.snps; do
    src="$OUT_DIR/filtered/${label}.vcf.gz"
    out="$OUT_DIR/stats/${label}.bcftools_stats.txt"
    [[ -f "$out" ]] || bcftools stats -s - "$src" > "$out"
done
# Unfiltered for comparison in MultiQC
out="$OUT_DIR/stats/unfiltered.bcftools_stats.txt"
[[ -f "$out" ]] || bcftools stats -s - "$IN_VCF" > "$out"

# =====================================================================
#     Step 5: SFS + singleton/doubleton counts on biallelic SNPs
# =====================================================================
echo "[$RUN] SFS / singletons..."
SFS_DIR="$OUT_DIR/sfs"

if [[ ! -f "$SFS_DIR/unfolded_sfs.tsv" ]]; then
    {
        printf 'AC\tcount\n'
        bcftools query -f '%AC\n' "$BIALL_VCF" | \
            awk '$1!="." && $1!~/,/ {n[$1]++} END {for (k in n) print k"\t"n[k]}' | \
            sort -n -k1,1
    } > "$SFS_DIR/unfolded_sfs.tsv"
fi

if [[ ! -f "$SFS_DIR/folded_sfs.tsv" ]]; then
    {
        printf 'MAC\tcount\n'
        bcftools query -f '%AC\t%AN\n' "$BIALL_VCF" | \
            awk '$1!="." && $2!="." && $1!~/,/ {
                ac=$1; an=$2; mac=(ac<=an/2)?ac:an-ac; n[mac]++
            } END {
                for (k in n) print k"\t"n[k]
            }' | sort -n -k1,1
    } > "$SFS_DIR/folded_sfs.tsv"
fi

# Summary numbers (singletons = AC==1 or AC==(AN-1) on biallelic sites)
{
    printf 'metric\tvalue\n'
    printf 'biallelic_snps_total\t%s\n' "$(bcftools view -H "$BIALL_VCF" | wc -l)"
    printf 'singletons\t%s\n'           "$(bcftools view -H -i '(AC==1 || AC==(AN-1))' "$BIALL_VCF" | wc -l)"
    printf 'doubletons\t%s\n'           "$(bcftools view -H -i '(AC==2 || AC==(AN-2))' "$BIALL_VCF" | wc -l)"
} > "$SFS_DIR/summary.tsv"

# =====================================================================
#     Step 6: pixy (pi, fst, dxy) on three filter levels
# =====================================================================
# Pixy is run on each of:
#   unfiltered           — raw VCF from grenepipe (has invariant sites)
#   qual20               — QUAL>=20 filter (has invariant sites; pixy's "recommended" input)
#   qual20.biallelic.snps — biallelic SNPs only (NO invariant sites!)
#
# IMPORTANT: pixy on a SNPs-only VCF gives BIASED pi/dxy because invariant sites
# are missing from the denominator (count_comparisons in pixy's output). The
# resulting pi will be inflated relative to the true value. Useful only for
# *relative* comparison (e.g. "filter X reduces apparent pi by N%"); not as an
# absolute estimate.  Take the qual20 numbers as your canonical pi/dxy.
PIXY_DIR="$OUT_DIR/pixy"
CHROM_BED="$PIXY_DIR/chroms.bed"
mkdir -p "$PIXY_DIR"
awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' "$FAI" > "$CHROM_BED"

if [[ "${SKIP_PIXY:-0}" == "1" ]]; then
    echo "[$RUN] pixy: SKIPPED (SKIP_PIXY=1)"
else
    if [[ "$PIXY_ENV" != "$CONDA_ENV" ]]; then
        echo "[$RUN] activating pixy env: $PIXY_ENV"
        conda activate "$PIXY_ENV"
    fi

    if ! command -v pixy >/dev/null 2>&1; then
        echo "WARN [$RUN]: pixy not found in $CONDA_DEFAULT_ENV — skipping pixy section." >&2
    else
        # All three inputs are now normalized (Step 0). The "unfiltered" label
        # means "no QUAL or biallelic filter applied" — but still normalized so
        # cross-caller pi/dxy comparisons are fair.
        declare -A PIXY_INPUTS=(
            [unfiltered]="$NORM_VCF"
            [qual20]="$QUAL_VCF"
            [qual20.biallelic.snps]="$BIALL_VCF"
        )

        for label in unfiltered qual20 qual20.biallelic.snps; do
            src="${PIXY_INPUTS[$label]}"
            dst="$PIXY_DIR/$label"
            mkdir -p "$dst"

            # Biallelic-SNPs-only VCFs have no invariant sites by construction.
            # pixy 2.0 refuses to run on these without an explicit bypass; results
            # are intentionally biased (see Step 6 warning in the script header),
            # included for relative cross-run comparison only.
            # NOTE: --bypass_invariant_check is a boolean flag in pixy 2.0 — passing
            # a value (e.g. `yes`) gives "unrecognized arguments: yes".
            pixy_extra=()
            if [[ "$label" == "qual20.biallelic.snps" ]]; then
                pixy_extra=(--bypass_invariant_check)
            fi

            # Sliding-window pixy — Hudson FST (--fst_type hudson) for each
            # configured window size (10kb + 50kb by default).
            for win_bp in $PIXY_WINDOW_SIZES; do
                # Human-readable label: 10000 -> 10kb, 50000 -> 50kb, 1000000 -> 1Mb
                win_lbl=$(awk -v n=$win_bp 'BEGIN{
                    if (n>=1000000) printf "%dMb", n/1000000
                    else if (n>=1000) printf "%dkb", n/1000
                    else printf "%dbp", n
                }')
                prefix="windows_${win_lbl}"
                if [[ ! -f "$dst/${prefix}_pi.txt" ]]; then
                    echo "[$RUN] pixy $label / $win_lbl..."
                    pixy --stats pi fst dxy \
                         --fst_type hudson \
                         "${pixy_extra[@]}" \
                         --vcf "$src" \
                         --populations "$POP_TSV" \
                         --window_size "$win_bp" \
                         --n_cores "$PP_CORES" \
                         --output_folder "$dst" \
                         --output_prefix "$prefix" \
                         > "$OUT_DIR/logs/pixy_${label}_${win_lbl}.log" 2>&1 || \
                         echo "WARN [$RUN]: pixy $win_lbl on $label failed; see logs/pixy_${label}_${win_lbl}.log" >&2
                fi
            done

            # Per-chromosome (one full-length BED entry per chrom -> one row per chrom per stat)
            if [[ ! -f "$dst/per_chrom_pi.txt" ]]; then
                echo "[$RUN] pixy $label / per-chrom..."
                pixy --stats pi fst dxy \
                     --fst_type hudson \
                     "${pixy_extra[@]}" \
                     --vcf "$src" \
                     --populations "$POP_TSV" \
                     --bed_file "$CHROM_BED" \
                     --n_cores "$PP_CORES" \
                     --output_folder "$dst" \
                     --output_prefix per_chrom \
                     > "$OUT_DIR/logs/pixy_${label}_perchrom.log" 2>&1 || \
                     echo "WARN [$RUN]: pixy per-chrom on $label failed; see logs/pixy_${label}_perchrom.log" >&2
            fi

            # Whole-genome pi/dxy: aggregate per-chrom (exact for these stats).
            # Uses simple variable-index awk syntax instead of $col["..."] field
            # references — easier to read and avoids any shell-parsing
            # ambiguities around `$col[...]` inside the awk literal.
            for stat in pi dxy; do
                src_f="$dst/per_chrom_${stat}.txt"
                dst_f="$dst/whole_genome_${stat}.tsv"
                [[ -f "$src_f" && ! -f "$dst_f" ]] || continue
                awk -F'\t' '
                NR == 1 {
                    # Discover the column indices we need from the header.
                    for (i = 1; i <= NF; i++) col_idx[$i] = i
                    has_pair = ("pop1" in col_idx)
                    if (has_pair) {
                        p1_col = col_idx["pop1"]
                        p2_col = col_idx["pop2"]
                    } else {
                        p_col  = col_idx["pop"]
                    }
                    diffs_col = col_idx["count_diffs"]
                    comps_col = col_idx["count_comparisons"]
                    print (has_pair ? "pop1\tpop2" : "pop") "\tcount_diffs\tcount_comparisons\tvalue"
                    next
                }
                {
                    key = has_pair ? ($p1_col "\t" $p2_col) : $p_col
                    diffs[key] += $diffs_col
                    comps[key] += $comps_col
                }
                END {
                    for (key in diffs) {
                        v = (comps[key] > 0) ? diffs[key] / comps[key] : "NA"
                        print key "\t" diffs[key] "\t" comps[key] "\t" v
                    }
                }' "$src_f" > "$dst_f"
            done

            # Whole-genome FST aggregation (APPROXIMATE — weighted mean by no_snps).
            # Unlike pi/dxy, FST cannot be exactly aggregated from per-window/per-
            # chrom values because pixy doesn't expose the per-site numerator and
            # denominator components. We emit a weighted mean by SNP count, which
            # is the most common back-of-envelope aggregation. The proper
            # ratio-of-averages value is in vcftools_fst/fst_summary.tsv
            # (column `weighted_fst`); this file is here for sanity-checking that
            # pixy and vcftools roughly agree per pop-pair.
            src_f="$dst/per_chrom_fst.txt"
            dst_f="$dst/whole_genome_fst.tsv"
            if [[ -f "$src_f" && ! -f "$dst_f" ]]; then
                awk -F'\t' '
                NR == 1 {
                    for (i = 1; i <= NF; i++) col_idx[$i] = i
                    # pixy 2.0 keeps the column name `avg_wc_fst` for both W&C and
                    # Hudson — but allow `avg_hudson_fst` too in case that changes.
                    fst_col  = ("avg_hudson_fst" in col_idx) ? col_idx["avg_hudson_fst"] \
                                                             : col_idx["avg_wc_fst"]
                    snps_col = col_idx["no_snps"]
                    p1_col   = col_idx["pop1"]
                    p2_col   = col_idx["pop2"]
                    print "pop1\tpop2\tweighted_mean_fst\ttotal_snps\testimator\tnote"
                    next
                }
                {
                    fst = $fst_col; n = $snps_col
                    if (fst == "NA" || n == "NA" || n + 0 == 0) next
                    key = $p1_col "\t" $p2_col
                    num[key] += fst * n
                    den[key] += n
                }
                END {
                    for (key in den) {
                        v = (den[key] > 0) ? num[key] / den[key] : "NA"
                        # Estimator label: pixy was invoked with --fst_type hudson
                        # in Step 6; column name is still avg_wc_fst per pixy 2.0
                        # backwards-compat. Document this explicitly here.
                        print key "\t" v "\t" den[key] "\tHudson\tweighted_mean_by_no_snps_APPROX_compare_with_vcftools_weighted_fst"
                    }
                }' "$src_f" > "$dst_f"
            fi
        done
    fi
fi

# =====================================================================
#     Step 7: vcftools Weir & Cockerham FST (proper ratio-of-averages)
# =====================================================================
echo "[$RUN] vcftools FST..."
VCFT_DIR="$OUT_DIR/vcftools_fst"
mkdir -p "$VCFT_DIR/pops"

# vcftools lives in a Lunarc lmod module (GCC + VCFtools). Load it if it's not
# already on PATH. Lmod is normally exported into sbatch jobs via BASH_ENV; if
# that's not the case, source the init script as a fallback.
if [[ -n "$VCFTOOLS_MODULE" ]] && ! command -v vcftools >/dev/null 2>&1; then
    if ! command -v module >/dev/null 2>&1 && [[ -f /usr/share/lmod/lmod/init/bash ]]; then
        source /usr/share/lmod/lmod/init/bash
    fi
    if command -v module >/dev/null 2>&1; then
        module load $VCFTOOLS_MODULE 2>/dev/null || \
            echo "WARN [$RUN]: failed: module load $VCFTOOLS_MODULE" >&2
    fi
fi

if ! command -v vcftools >/dev/null 2>&1; then
    echo "WARN [$RUN]: vcftools not on PATH (tried module: $VCFTOOLS_MODULE) — skipping Step 7." >&2
elif [[ ! -s "$POP_TSV" ]] || [[ $(awk -F'\t' '{print $2}' "$POP_TSV" | sort -u | wc -l) -lt 2 ]]; then
    echo "WARN [$RUN]: fewer than 2 populations — skipping Step 7 (FST needs >=2 pops)." >&2
else
    # Write per-population sample lists (one sample per line). POP_TSV has no header.
    awk -F'\t' -v outdir="$VCFT_DIR/pops" \
        '$1!="" && $2!="" {print $1 > outdir"/"$2".txt"}' "$POP_TSV"

    mapfile -t POPS < <(ls "$VCFT_DIR/pops"/ | sed 's/\.txt$//' | sort)
    n_pairs=$(( ${#POPS[@]} * (${#POPS[@]} - 1) / 2 ))
    echo "[$RUN] vcftools: ${#POPS[@]} pops -> $n_pairs pairwise FST(s)"

    # Helper: extract overall mean+weighted FST from a vcftools log
    parse_log() {
        local log="$1"
        local mean weighted
        mean=$(grep -m1 'Weir and Cockerham mean Fst estimate'     "$log" 2>/dev/null | awk -F': ' '{print $NF}')
        weighted=$(grep -m1 'Weir and Cockerham weighted Fst estimate' "$log" 2>/dev/null | awk -F': ' '{print $NF}')
        printf '%s\t%s' "${mean:-NA}" "${weighted:-NA}"
    }

    SUMMARY="$VCFT_DIR/fst_summary.tsv"
    {
        printf 'pop1\tpop2\tregion\tn_snps\tmean_fst\tweighted_fst\n'

        for ((i=0; i<${#POPS[@]}; i++)); do
            for ((j=i+1; j<${#POPS[@]}; j++)); do
                p1="${POPS[i]}"; p2="${POPS[j]}"
                pair="${p1}_vs_${p2}"

                # 7a. Whole-genome (no window) — single overall FST per pair
                outpre="$VCFT_DIR/${pair}_whole_genome"
                if [[ ! -f "${outpre}.log" ]]; then
                    vcftools --gzvcf "$BIALL_VCF" \
                        --weir-fst-pop "$VCFT_DIR/pops/${p1}.txt" \
                        --weir-fst-pop "$VCFT_DIR/pops/${p2}.txt" \
                        --out "$outpre" \
                        > /dev/null 2> "${outpre}.log"
                fi
                n=$(awk 'NR>1{c++} END{print c+0}' "${outpre}.weir.fst" 2>/dev/null)
                printf '%s\t%s\twhole_genome\t%s\t%s\n' "$p1" "$p2" "${n:-0}" "$(parse_log "${outpre}.log")"

                # 7b. Sliding windows — one set per configured window size
                for win_bp in $PIXY_WINDOW_SIZES; do
                    win_lbl=$(awk -v n=$win_bp 'BEGIN{
                        if (n>=1000000) printf "%dMb", n/1000000
                        else if (n>=1000) printf "%dkb", n/1000
                        else printf "%dbp", n
                    }')
                    outpre="$VCFT_DIR/${pair}_windows_${win_lbl}"
                    if [[ ! -f "${outpre}.windowed.weir.fst" ]]; then
                        vcftools --gzvcf "$BIALL_VCF" \
                            --weir-fst-pop "$VCFT_DIR/pops/${p1}.txt" \
                            --weir-fst-pop "$VCFT_DIR/pops/${p2}.txt" \
                            --fst-window-size "$win_bp" \
                            --out "$outpre" \
                            > /dev/null 2> "${outpre}.log"
                    fi
                done

                # 7c. Per-chromosome — reuse the per-chrom biallelic VCFs from Step 3
                for chrom in "${CHROMS[@]}"; do
                    safe=${chrom//\//_}; safe=${safe//:/_}
                    chrom_vcf="$OUT_DIR/per_chrom/qual20.biallelic.snps/${safe}.vcf.gz"
                    [[ -f "$chrom_vcf" ]] || continue
                    outpre="$VCFT_DIR/${pair}_chr_${safe}"
                    if [[ ! -f "${outpre}.log" ]]; then
                        vcftools --gzvcf "$chrom_vcf" \
                            --weir-fst-pop "$VCFT_DIR/pops/${p1}.txt" \
                            --weir-fst-pop "$VCFT_DIR/pops/${p2}.txt" \
                            --out "$outpre" \
                            > /dev/null 2> "${outpre}.log"
                    fi
                    n=$(awk 'NR>1{c++} END{print c+0}' "${outpre}.weir.fst" 2>/dev/null)
                    printf '%s\t%s\t%s\t%s\t%s\n' "$p1" "$p2" "$chrom" "${n:-0}" "$(parse_log "${outpre}.log")"
                done
            done
        done
    } > "$SUMMARY"

    echo "[$RUN] vcftools FST summary: $SUMMARY"
fi

# =====================================================================
#     Step 8: bcftools csq annotation (uses species-matched GFF)
# =====================================================================
# Annotates the biallelic-SNPs VCF with:
#   (a) BCSQ — consequence terms from bcftools csq (missense, stop_gained, ...)
#   (b) IMPACT — derived SnpEff/VEP-style severity bucket (HIGH/MODERATE/LOW/MODIFIER),
#       set to the most severe consequence among all BCSQ entries for each site.
# Records that bcftools csq doesn't annotate (no transcript overlap) get IMPACT=MODIFIER.
#
# The GFF is matched to the reference by GCA accession parsed from the ref
# filename (e.g. GCA_937595015.1 in GCA_937595015.1_ilPolIcar1.1_genomic.fna).
echo "[$RUN] annotation (bcftools csq)..."
ANNOT_DIR="$OUT_DIR/annotated"
mkdir -p "$ANNOT_DIR"

REF_DIR=$(dirname "$REF")
ANN_DIR="$REF_DIR/lepeu_annotations"
REF_ACC=$(basename "$REF" | grep -oP 'GCA_[0-9]+\.[0-9]+' || true)
GFF=""
if [[ -n "$REF_ACC" && -d "$ANN_DIR" ]]; then
    GFF=$(ls "$ANN_DIR"/*"$REF_ACC"*.gff 2>/dev/null | head -1)
fi

if [[ -z "$GFF" || ! -f "$GFF" ]]; then
    echo "WARN [$RUN]: no GFF found for $REF_ACC in $ANN_DIR — skipping Step 8." >&2
else
    ANNOT_VCF="$ANNOT_DIR/qual20.biallelic.snps.annotated.vcf.gz"
    BCSQ_VCF="$ANNOT_DIR/qual20.biallelic.snps.bcsq.vcf.gz"   # intermediate (BCSQ only, no IMPACT)
    IMPACT_TAB="$ANNOT_DIR/impact.tsv.gz"                       # CHROM POS REF ALT IMPACT

    # Treat the BCSQ intermediate as "complete" only when both the .vcf.gz AND
    # the .tbi exist — a crash mid-write of the BCSQ output leaves a truncated
    # .vcf.gz without an index, which would otherwise pass a naive existence
    # check and feed corrupted data into Step 8b.
    bcsq_complete=0
    if [[ -f "$BCSQ_VCF" && -f "${BCSQ_VCF}.tbi" ]]; then
        bcsq_complete=1
    elif [[ -f "$BCSQ_VCF" ]]; then
        echo "WARN [$RUN]: partial BCSQ_VCF (no .tbi) — discarding and re-running Step 8a." >&2
        rm -f "$BCSQ_VCF" "${BCSQ_VCF}.tbi"
    fi

    # --- Step 8a: run bcftools csq -> intermediate VCF with BCSQ ---
    # --force keeps csq running past VCF chromosomes that aren't in the GFF
    # (the AnnEvo GFFs only cover annotated scaffolds; variants on unannotated
    # scaffolds / unplaced contigs just pass through without a BCSQ entry).
    if (( ! bcsq_complete )) && [[ ! -f "$ANNOT_VCF" ]]; then
        echo "[$RUN] running bcftools csq with $(basename "$GFF")..."
        bcftools csq \
            --fasta-ref "$REF" \
            --gff-annot "$GFF" \
            --phase a \
            --force \
            --threads "$IOTH" \
            -O z -o "$BCSQ_VCF" \
            "$BIALL_VCF" \
            > "$OUT_DIR/logs/csq.log" 2>&1 || \
            echo "WARN [$RUN]: bcftools csq failed; see logs/csq.log" >&2
        if [[ -f "$BCSQ_VCF" ]]; then
            bcftools index --tbi --threads 2 "$BCSQ_VCF" 2>/dev/null || {
                # If indexing fails the BCSQ file is truncated — discard so the
                # next run regenerates cleanly.
                echo "WARN [$RUN]: BCSQ index failed (truncated csq output) — removing" >&2
                rm -f "$BCSQ_VCF" "${BCSQ_VCF}.tbi"
            }
            [[ -f "${BCSQ_VCF}.tbi" ]] && bcsq_complete=1
        fi
    fi

    # --- Step 8b: derive IMPACT from BCSQ, write a CHROM/POS/REF/ALT/IMPACT lookup ---
    # SnpEff/VEP-style severity mapping for bcftools csq consequence vocabulary.
    # Compound prefixes (* for compound consequence, @<pos> position marker) are stripped first.
    # Same truncation guard on IMPACT_TAB: only consider it done if both .gz and .tbi are present.
    if [[ -f "$IMPACT_TAB" && ! -f "${IMPACT_TAB}.tbi" ]]; then
        echo "WARN [$RUN]: partial impact.tsv.gz (no .tbi) — discarding and re-deriving." >&2
        rm -f "$IMPACT_TAB" "${IMPACT_TAB}.tbi"
    fi
    if (( bcsq_complete )) && [[ ! -f "$IMPACT_TAB" ]]; then
        echo "[$RUN] deriving IMPACT from BCSQ..."
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/BCSQ\n' "$BCSQ_VCF" | \
        awk -F'\t' 'BEGIN {
            OFS="\t"
            # consequence -> impact bucket
            split("stop_gained stop_lost start_lost splice_acceptor splice_donor frameshift", a, " ")
            for (i in a) high[a[i]] = 1
            split("missense inframe_insertion inframe_deletion", a, " ")
            for (i in a) moderate[a[i]] = 1
            split("synonymous start_retained stop_retained splice_region", a, " ")
            for (i in a) low[a[i]] = 1
            # everything else -> MODIFIER
            rank["HIGH"]=4; rank["MODERATE"]=3; rank["LOW"]=2; rank["MODIFIER"]=1
        }
        function cons_to_impact(c) {
            sub(/^\*/, "", c)            # strip compound-consequence prefix
            sub(/@[0-9]+$/, "", c)       # strip @<pos> marker
            if (c in high)     return "HIGH"
            if (c in moderate) return "MODERATE"
            if (c in low)      return "LOW"
            return "MODIFIER"
        }
        function most_severe(bcsq,    n, i, parts, fields, best, cur) {
            if (bcsq == "." || bcsq == "") return "MODIFIER"
            n = split(bcsq, parts, ",")
            best = "MODIFIER"
            for (i=1; i<=n; i++) {
                split(parts[i], fields, "|")
                cur = cons_to_impact(fields[1])
                if (rank[cur] > rank[best]) best = cur
            }
            return best
        }
        {
            print $1, $2, $3, $4, most_severe($5)
        }' | bgzip > "$IMPACT_TAB"
        tabix -s 1 -b 2 -e 2 "$IMPACT_TAB"
    fi

    # --- Step 8c: bcftools annotate the IMPACT field onto the BCSQ VCF ---
    # Truncation guard on ANNOT_VCF: missing .tbi means a previous run crashed
    # mid-write — discard and regenerate.
    if [[ -f "$ANNOT_VCF" && ! -f "${ANNOT_VCF}.tbi" ]]; then
        echo "WARN [$RUN]: partial ANNOT_VCF (no .tbi) — discarding and re-running Step 8c." >&2
        rm -f "$ANNOT_VCF" "${ANNOT_VCF}.tbi"
    fi
    if (( bcsq_complete )) && [[ -f "$IMPACT_TAB" && ! -f "$ANNOT_VCF" ]]; then
        echo "[$RUN] adding IMPACT INFO field..."
        bcftools annotate \
            -a "$IMPACT_TAB" \
            -c CHROM,POS,REF,ALT,INFO/IMPACT \
            -h <(echo '##INFO=<ID=IMPACT,Number=1,Type=String,Description="SnpEff/VEP-style severity bucket (HIGH/MODERATE/LOW/MODIFIER), derived from BCSQ. Most severe of all consequences per record.">') \
            --threads "$IOTH" \
            "$BCSQ_VCF" \
            -O z -o "$ANNOT_VCF"
        bcftools index --tbi --threads 2 "$ANNOT_VCF"
        # The BCSQ-only intermediate is now redundant — clean it up
        rm -f "$BCSQ_VCF" "$BCSQ_VCF.tbi"
    fi

    # --- Step 8d: consequence-level and impact-level count summaries ---
    if [[ -f "$ANNOT_VCF" && ! -f "$ANNOT_DIR/csq_consequence_counts.tsv" ]]; then
        {
            printf 'consequence\tcount\n'
            bcftools query -f '%INFO/BCSQ\n' "$ANNOT_VCF" 2>/dev/null \
                | tr ',' '\n' \
                | awk -F'|' 'NF>1 {c=$1; sub(/^\*/, "", c); sub(/@[0-9]+$/, "", c); print c}' \
                | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2 -rn
        } > "$ANNOT_DIR/csq_consequence_counts.tsv"
        echo "[$RUN] consequence summary: $ANNOT_DIR/csq_consequence_counts.tsv"
    fi

    if [[ -f "$ANNOT_VCF" && ! -f "$ANNOT_DIR/csq_impact_counts.tsv" ]]; then
        {
            printf 'impact\tcount\n'
            bcftools query -f '%INFO/IMPACT\n' "$ANNOT_VCF" 2>/dev/null \
                | awk '{print ($1=="."?"MODIFIER":$1)}' \
                | sort | uniq -c | awk '{print $2"\t"$1}' \
                | awk 'BEGIN{r["HIGH"]=4;r["MODERATE"]=3;r["LOW"]=2;r["MODIFIER"]=1} {print r[$1]"\t"$0}' \
                | sort -k1,1 -rn | cut -f2-
        } > "$ANNOT_DIR/csq_impact_counts.tsv"
        echo "[$RUN] impact summary: $ANNOT_DIR/csq_impact_counts.tsv"
    fi
fi

touch "$OUT_DIR/.done"
echo "[$RUN] DONE $(date -Iseconds)"
echo "       output: $OUT_DIR"
