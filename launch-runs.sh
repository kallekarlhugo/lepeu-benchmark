#!/usr/bin/env bash
#
# launch-runs.sh — submit grenepipe runs for the LepEU benchmarking matrix.
#
# All-in-one-sbatch strategy (HPC admins asked us not to use per-rule SLURM
# submission):
#   * one sbatch per pipeline run,
#   * snakemake runs with --executor local inside the allocation,
#   * --cores controls per-rule parallelism inside the allocation.
#
# CPU runs go to the default (lu48) partition with 48 cores / 240G.
# DeepVariant runs route to the shared GPU partition gpua100i (32 cores,
# 1x A100-40G), see the "GPU routing" block in the loop.
#
# Two-phase strategy:
#
#   ./launch-runs.sh --anchors-only       # Phase 1: one job per species*aligner.
#                                         # Each runs trim + map + dedup + bcftools_cohort.
#
#   <wait for those to finish: squeue -u $USER>
#
#   ./launch-runs.sh --dependents-only    # Phase 2: everything else.
#                                         # Each reuses the anchor's BAMs via mappings-table.
#
# All runs get:
#   keep-invariant-sites: "true"   — emit invariant + variant where supported.
#   filter-variants:      "none"   — unfiltered output; filter downstream.
#
#   DRY_RUN=1 ./launch-runs.sh --anchors-only      # show sbatch commands without submitting

set -euo pipefail

# =====================================================================
#     Edit these for your Lunarc environment
# =====================================================================
# Path to the modded grenepipe (with NGM/DV patches + invariant-sites flag fix).
# Moved out of $HOME so other workshop attendees can read it.
GRENEPIPE_DIR="/lunarc/nobackup/projects/lepeu-lisbon/kalle_temp/LepEU/grenepipe"
RUNS_ROOT="$(pwd -P)/runs"                # -P resolves symlinks; Singularity binds need real paths
PROFILE_REL="profiles/lunarc-slurm"       # relative to GRENEPIPE_DIR

# Master sbatch wrapper resources — the master IS the worker now
# (snakemake --executor local), so this is the full envelope of one
# pipeline run.  CPU defaults below; DeepVariant runs override these
# via the GPU-routing block in the loop.
SLURM_ACCOUNT="lu2026-2-62"
MASTER_CORES=48
MASTER_MEM="240G"
MASTER_TIME="96:00:00"

# Cap on concurrent rules per master (snakemake --jobs).  With --executor local
# this just gates parallelism inside the allocation, not SLURM submission.
JOBS_PER_MASTER=9

# Conda
CONDA_BASE="$HOME/miniconda3"
CONDA_ENV="/lunarc/nobackup/projects/lepeu-lisbon/shared/conda/grenepipe"

# Set to a module name only if your shell does not already have singularity on PATH.
SINGULARITY_MODULE=""
# =====================================================================

usage() {
    cat >&2 <<EOF
Usage:
  $(basename "$0") --anchors-only       [runs.tsv]
  $(basename "$0") --dependents-only    [runs.tsv]

Run anchors first; once they all finish, run dependents.

Environment:
  DRY_RUN=1       Print sbatch commands without submitting.
EOF
}

PHASE=""
RUNS_TSV="runs.tsv"
for arg in "$@"; do
    case "$arg" in
        --anchors-only)    PHASE="anchors" ;;
        --dependents-only) PHASE="dependents" ;;
        -h|--help)         usage; exit 0 ;;
        *.tsv)             RUNS_TSV="$arg" ;;
        *)                 echo "ERROR: unknown arg '$arg'" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$PHASE" ]]; then
    echo "ERROR: pass --anchors-only or --dependents-only" >&2
    usage
    exit 1
fi

if [[ ! -d "$GRENEPIPE_DIR" ]]; then
    echo "ERROR: GRENEPIPE_DIR does not exist: $GRENEPIPE_DIR" >&2
    exit 1
fi
if [[ ! -d "$GRENEPIPE_DIR/$PROFILE_REL" ]]; then
    echo "ERROR: profile not found at $GRENEPIPE_DIR/$PROFILE_REL" >&2
    exit 1
fi
if [[ ! -f "$RUNS_TSV" ]]; then
    echo "ERROR: runs file not found: $RUNS_TSV" >&2
    exit 1
fi

mkdir -p "$RUNS_ROOT"

submitted=0
skipped=0

# Read runs.tsv: tab-separated, columns:
#   run_id  reference  samples_table  mapping_tool  calling_tool  calling_mode  bams_from
while IFS=$'\t' read -r run_id reference samples_table mapping_tool calling_tool calling_mode bams_from _rest; do
    case "$run_id" in
        ""|"run_id"|"#"*) continue ;;
    esac

    # Anchor = empty bams_from. Dependent = has a bams_from value (an upstream run_id).
    is_anchor="no"
    if [[ -z "${bams_from:-}" || "$bams_from" == "-" ]]; then
        is_anchor="yes"
    fi

    if [[ "$PHASE" == "anchors" && "$is_anchor" == "no" ]]; then
        skipped=$((skipped + 1))
        continue
    fi
    if [[ "$PHASE" == "dependents" && "$is_anchor" == "yes" ]]; then
        skipped=$((skipped + 1))
        continue
    fi

    run_dir="$RUNS_ROOT/$run_id"
    mkdir -p "$run_dir"

    # Copy default config and override the per-run fields.
    cp "$GRENEPIPE_DIR/config/config.yaml" "$run_dir/config.yaml"
    sed -i \
        -e "s|^  reference-genome:.*|  reference-genome: \"$reference\"|" \
        -e "s|^  samples-table:.*|  samples-table: \"$samples_table\"|" \
        -e "s|^  mapping-tool:.*|  mapping-tool: \"$mapping_tool\"|" \
        -e "s|^  calling-tool:.*|  calling-tool: \"$calling_tool\"|" \
        -e 's|^  keep-invariant-sites:.*|  keep-invariant-sites: "true"|' \
        -e 's|^  filter-variants:.*|  filter-variants: "none"|' \
        "$run_dir/config.yaml"

    # Bcftools mode override (only meaningful when calling-tool == bcftools).
    case "$calling_mode" in
        individual|combined)
            sed -i -e "s|^    mode: \".*\"|    mode: \"$calling_mode\"|" \
                "$run_dir/config.yaml"
            ;;
        combined_hwe)
            # Run bcftools in combined mode but with per-population HWE prior.
            # We do this without patching grenepipe: bcftools call's free-form
            # extra-arg string (params.bcftools.call) already gets passed through,
            # so we append `-G <populations.tsv>` to it.
            sed -i -e "s|^    mode: \".*\"|    mode: \"combined\"|" \
                "$run_dir/config.yaml"
            awk -F'\t' 'NR>1 && $1!="" && !seen[$1]++ {if ($3!="") print $1"\t"$3}' \
                "$samples_table" > "$run_dir/populations.tsv"
            if [[ ! -s "$run_dir/populations.tsv" ]]; then
                echo "ERROR [$run_id]: populations.tsv empty; samples table needs a non-empty population column" >&2
                exit 1
            fi
            # Append -G inside the bcftools block only.
            sed -i "/^  bcftools:/,/^  [a-zA-Z]/{
                s|^\\(    call: \".*\\)\"|\\1 -G $run_dir/populations.tsv\"|
            }" "$run_dir/config.yaml"
            ;;
        combined_nohwe)
            # Run bcftools in combined mode but explicitly without the HWE prior.
            # `-G -` tells bcftools to treat each sample as its own group, which
            # effectively disables the HWE assumption in the multiallelic caller
            # (a "group" of one sample has no HWE to enforce).
            sed -i -e "s|^    mode: \".*\"|    mode: \"combined\"|" \
                "$run_dir/config.yaml"
            sed -i "/^  bcftools:/,/^  [a-zA-Z]/{
                s|^\\(    call: \".*\\)\"|\\1 -G -\"|
            }" "$run_dir/config.yaml"
            ;;
        ""|default)
            ;;
        *)
            echo "ERROR [$run_id]: unknown calling_mode '$calling_mode'" >&2
            exit 1
            ;;
    esac

    # ----- Dependent setup: point at anchor's BAMs via mappings-table -----
    if [[ "$is_anchor" == "no" ]]; then
        anchor_dir="$RUNS_ROOT/$bams_from"
        if [[ ! -d "$anchor_dir" ]]; then
            echo "WARNING [$run_id]: anchor dir does not exist yet: $anchor_dir" >&2
            echo "         (Run --anchors-only first and wait for them to finish.)" >&2
        fi

        # initialize-bam.smk needs two columns: sample\tbam.
        {
            printf 'sample\tbam\n'
            awk -F'\t' 'NR>1 && $1 != "" && !seen[$1]++ {print $1}' "$samples_table" |
            while read -r sample; do
                printf '%s\t%s\n' "$sample" "$anchor_dir/mapping/final/$sample.bam"
            done
        } > "$run_dir/mappings.tsv"

        sed -i -e "s|^  mappings-table:.*|  mappings-table: \"$run_dir/mappings.tsv\"|" \
            "$run_dir/config.yaml"
    fi

    # =====================================================================
    # GPU routing — MUST come BEFORE the payload heredoc.
    # `set -u` makes unbound $job_cores fatal at heredoc-evaluation time, so
    # these vars need values *before* we build `payload="..."` below.
    # =====================================================================
    if [[ "$calling_tool" == "deepvariant" ]]; then
        # gpua100i: shared GPU partition, 32 cores, 1x A100-40G per job.
        job_partition="gpua100i"
        job_gres="gpu:1"
        job_cores=32
        job_mem="120G"
        job_time="24:00:00"
        # Raise per-user process cap inside the payload: make_examples forks
        # heavily and hits RLIMIT_NPROC on default Lunarc settings.
        dv_ulimit="ulimit -u \$(ulimit -Hu) || true"
    else
        job_partition=""      # default partition (lu48)
        job_gres=""
        job_cores=$MASTER_CORES
        job_mem=$MASTER_MEM
        job_time=$MASTER_TIME
        dv_ulimit=":"
    fi

    # ----- Build the master sbatch payload -----
    if [[ "$is_anchor" == "yes" ]]; then
        # Anchors: full pipeline, then materialize mapping/final/<sample>.bam
        # symlinks + mappings.tsv that dependents read.
        payload="set -euo pipefail
source '$CONDA_BASE/etc/profile.d/conda.sh'
conda activate '$CONDA_ENV'
[[ -n '$SINGULARITY_MODULE' ]] && module load '$SINGULARITY_MODULE' || true
$dv_ulimit
cd '$GRENEPIPE_DIR'
snakemake \\
    --profile '$PROFILE_REL' \\
    --executor local \\
    --restart-times 3 \\
    --rerun-incomplete \\
    --jobs $JOBS_PER_MASTER \\
    --cores $job_cores \\
    --directory '$run_dir' \\
    --configfile '$run_dir/config.yaml'
snakemake \\
    all_bams \\
    --profile '$PROFILE_REL' \\
    --executor local \\
    --restart-times 3 \\
    --rerun-incomplete \\
    --jobs $JOBS_PER_MASTER \\
    --cores $job_cores \\
    --directory '$run_dir' \\
    --configfile '$run_dir/config.yaml'
# Symlink .bai next to each mapping/final/<sample>.bam.
for bam in '$run_dir'/mapping/final/*.bam; do
    real_bam=\$(readlink -f \"\$bam\")
    if [[ -f \"\${real_bam}.bai\" ]]; then
        ln -sf \"\${real_bam}.bai\" \"\${bam}.bai\"
    fi
done
# Write mappings.tsv from the final/ symlinks.
{
    printf 'sample\\tbam\\n'
    for bam in '$run_dir'/mapping/final/*.bam; do
        sample=\$(basename \"\$bam\" .bam)
        printf '%s\\t%s\\n' \"\$sample\" \"\$bam\"
    done
} > '$run_dir/mappings.tsv'
echo \"Anchor $run_id done. mappings.tsv has \$(wc -l < '$run_dir/mappings.tsv') lines.\""
    else
        # Dependent: mappings-table makes grenepipe skip trim/map/dedup.
        payload="set -euo pipefail
source '$CONDA_BASE/etc/profile.d/conda.sh'
conda activate '$CONDA_ENV'
[[ -n '$SINGULARITY_MODULE' ]] && module load '$SINGULARITY_MODULE' || true
$dv_ulimit
cd '$GRENEPIPE_DIR'
snakemake \\
    --profile '$PROFILE_REL' \\
    --executor local \\
    --restart-times 3 \\
    --rerun-incomplete \\
    --jobs $JOBS_PER_MASTER \\
    --cores $job_cores \\
    --directory '$run_dir' \\
    --configfile '$run_dir/config.yaml'"
    fi

    sbatch_args=(
        --account="$SLURM_ACCOUNT"
        --time="$job_time"
        --cpus-per-task="$job_cores"
        --mem="$job_mem"
        --job-name="gp-$run_id"
        --output="$run_dir/slurm-%j.out"
        --error="$run_dir/slurm-%j.err"
    )
    [[ -n "$job_partition" ]] && sbatch_args+=( --partition="$job_partition" )
    [[ -n "$job_gres" ]] && sbatch_args+=( --gres="$job_gres" )
    sbatch_args+=( --wrap="$payload" )

    label="anchor"
    [[ "$is_anchor" == "no" ]] && label="dependent of $bams_from"

    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        echo "=== $run_id ($label; aligner=$mapping_tool, caller=$calling_tool, mode=$calling_mode; partition=${job_partition:-default}, cores=$job_cores, mem=$job_mem) ==="
        echo "sbatch ${sbatch_args[*]}"
        echo
    else
        echo "Submitting $run_id ($label; aligner=$mapping_tool, caller=$calling_tool, mode=$calling_mode; partition=${job_partition:-default}, cores=$job_cores, mem=$job_mem)"
        sbatch "${sbatch_args[@]}"
    fi
    submitted=$((submitted + 1))
done < "$RUNS_TSV"

echo
echo "Phase: $PHASE. Submitted: $submitted. Skipped: $skipped (other phase)."
echo "Monitor with: squeue -u \$USER"
echo "Per-run outputs under: $RUNS_ROOT/"
if [[ "$PHASE" == "anchors" ]]; then
    echo
    echo "Next: wait for all anchor jobs to finish (squeue -u \$USER), then run:"
    echo "  $(basename "$0") --dependents-only"
fi
