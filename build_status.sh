#!/usr/bin/env bash
#
# build_status.sh — generate status.tsv: per-run progress across all pipeline steps.
#
# One row per run in runs.tsv. Columns are Y/N markers for each stage / step.
# Open status.tsv in Excel / LibreOffice / pandas etc., or just inspect with
# `column -t -s$'\t' status.tsv | less -S` on the cluster.
#
# Usage:
#   ./build_status.sh                  # reads runs.tsv, writes status.tsv
#   ./build_status.sh <runs.tsv>       # custom runs file
#   ./build_status.sh <runs.tsv> <output.tsv>

set -euo pipefail

RUNS_TSV="${1:-runs.tsv}"
OUT="${2:-status.tsv}"

RUNS_ROOT="$(pwd)/runs"
POST_ROOT="$(pwd)/postprocess"

[[ -f "$RUNS_TSV" ]] || { echo "ERROR: $RUNS_TSV not found" >&2; exit 1; }

# ----- helpers ---------------------------------------------------------------
# File-exists check
chk() { [[ -e "$1" ]] && echo "Y" || echo "N"; }

# Directory-has-VCFs check (for per-chrom subsets)
chkdir() {
    if [[ -d "$1" ]] && compgen -G "$1/*.vcf.gz" > /dev/null; then
        echo "Y"
    else
        echo "N"
    fi
}

# Stats file — handle both old (.bcftools_stats.txt) and new (.txt) naming
chk_stats() {
    local base="$1"
    [[ -f "${base}.txt" || -f "${base}.bcftools_stats.txt" ]] && echo Y || echo N
}

# ----- header ---------------------------------------------------------------
# Step abbreviations:
#   S1_calling           Stage 1 calling: runs/<run>/calling/genotyped-all.vcf.gz.done
#   S0_norm              Step 0: filtered/normalized.vcf.gz (post-norm + split multi-allelics)
#   S1_qual20            Step 1: filtered/qual20.vcf.gz (QUAL>=20)
#   S2_biallelic         Step 2: filtered/qual20.biallelic.snps.vcf.gz (biallelic SNPs only)
#   S3_pchrom_q20        Step 3: per-chrom VCFs for qual20
#   S3_pchrom_bi         Step 3: per-chrom VCFs for biallelic
#   S4_stats             Step 4: bcftools stats on biallelic VCF
#   S5_sfs               Step 5: SFS + singleton/doubleton counts
#   S6_pixy_unfilt       Step 6: pixy per-chrom on unfiltered (canonical existence marker)
#   S6_pixy_q20          Step 6: pixy per-chrom on qual20 (THE canonical pi/dxy level)
#   S6_pixy_bi           Step 6: pixy per-chrom on biallelic.snps (biased pi/dxy, kept for comparison)
#   S6_50kb_q20          Step 6: pixy 50 kb windowed at qual20 (recent addition)
#   S6_wg_fst_q20        Step 6: pixy whole_genome_fst.tsv at qual20 (approx weighted-mean Hudson FST)
#   S7_fst               Step 7: vcftools fst_summary.tsv (proper W&C ratio-of-averages FST)
#   S8_annot             Step 8: bcftools csq annotated VCF (.vcf.gz.tbi exists -> not truncated)
#   S8_impact            Step 8: IMPACT lookup TSV (impact.tsv.gz + .tbi)
#   PP_done              postprocess/<run>/.done sentinel (touched only on full success)
COLS=(
    run_id
    species aligner caller mode type
    S1_calling
    S0_norm S1_qual20 S2_biallelic
    S3_pchrom_q20 S3_pchrom_bi
    S4_stats S5_sfs
    S6_pixy_unfilt S6_pixy_q20 S6_pixy_bi
    S6_50kb_q20 S6_wg_fst_q20
    S7_fst
    S8_annot S8_impact
    PP_done
)
IFS=$'\t'; echo "${COLS[*]}" > "$OUT"; IFS=$' \t\n'

# ----- per-run rows ---------------------------------------------------------
# Skip header, blank lines, comments (#), and the literal "run_id" line.
awk -F'\t' 'NR>1 && $1!="" && $1!~/^#/ && $1!="run_id"' "$RUNS_TSV" | \
while IFS=$'\t' read -r run_id reference samples_table mapping_tool calling_tool calling_mode bams_from _rest; do
    # Parse run_id pieces (species_aligner_caller-possibly-with-underscores)
    species="${run_id%%_*}"
    rest="${run_id#*_}"
    aligner="${rest%%_*}"
    caller="${rest#*_}"

    # Anchor vs dependent
    if [[ -z "${bams_from:-}" || "$bams_from" == "-" ]]; then
        type="anchor"
    else
        type="dependent"
    fi

    r="$RUNS_ROOT/$run_id"
    p="$POST_ROOT/$run_id"

    # Stage 1
    s1_calling=$(chk "$r/calling/genotyped-all.vcf.gz.done")

    # Stage 2 per-step checks
    s0=$(chk "$p/filtered/normalized.vcf.gz")
    s1=$(chk "$p/filtered/qual20.vcf.gz")
    s2=$(chk "$p/filtered/qual20.biallelic.snps.vcf.gz")
    s3a=$(chkdir "$p/per_chrom/qual20")
    s3b=$(chkdir "$p/per_chrom/qual20.biallelic.snps")
    s4=$(chk_stats "$p/stats/qual20.biallelic.snps")
    s5=$(chk "$p/sfs/summary.tsv")
    s6a=$(chk "$p/pixy/unfiltered/per_chrom_pi.txt")
    s6b=$(chk "$p/pixy/qual20/per_chrom_pi.txt")
    s6c=$(chk "$p/pixy/qual20.biallelic.snps/per_chrom_pi.txt")
    s6_50kb=$(chk "$p/pixy/qual20/windows_50kb_pi.txt")
    s6_wgfst=$(chk "$p/pixy/qual20/whole_genome_fst.tsv")
    s7=$(chk "$p/vcftools_fst/fst_summary.tsv")
    # Step 8: annotated VCF needs both .gz AND .tbi to be considered complete
    if [[ -f "$p/annotated/qual20.biallelic.snps.annotated.vcf.gz" && \
          -f "$p/annotated/qual20.biallelic.snps.annotated.vcf.gz.tbi" ]]; then
        s8=Y
    else
        s8=N
    fi
    # IMPACT lookup: TSV + tabix index
    if [[ -f "$p/annotated/impact.tsv.gz" && -f "$p/annotated/impact.tsv.gz.tbi" ]]; then
        s8i=Y
    else
        s8i=N
    fi
    pp_done=$(chk "$p/.done")

    # Emit row (tab-separated)
    row=(
        "$run_id"
        "$species" "$aligner" "$caller" "${calling_mode:-default}" "$type"
        "$s1_calling"
        "$s0" "$s1" "$s2"
        "$s3a" "$s3b"
        "$s4" "$s5"
        "$s6a" "$s6b" "$s6c"
        "$s6_50kb" "$s6_wgfst"
        "$s7"
        "$s8" "$s8i"
        "$pp_done"
    )
    IFS=$'\t'; echo "${row[*]}"; IFS=$' \t\n'
done >> "$OUT"

# ----- stdout summary --------------------------------------------------------
n_rows=$(($(wc -l < "$OUT") - 1))

echo
echo "Wrote $OUT  ($n_rows runs)"
echo
echo "Completion by step (out of $n_rows runs):"
awk -F'\t' '
NR == 1 {
    for (i = 1; i <= NF; i++) hdr[i] = $i
    ncol = NF
    next
}
{
    for (i = 7; i <= NF; i++) if ($i == "Y") cnt[i]++
    total++
}
END {
    for (i = 7; i <= ncol; i++) {
        bar = ""
        bar_n = int(20 * cnt[i] / total)
        for (j = 0; j < bar_n; j++) bar = bar "#"
        for (j = bar_n; j < 20; j++) bar = bar "."
        printf "  %-20s [%s] %d/%d\n", hdr[i], bar, cnt[i]+0, total
    }
}' "$OUT"

echo
echo "Per-run progress (steps done / total):"
awk -F'\t' '
NR == 1 { ncol = NF; next }
{
    done = 0; total = 0
    for (i = 7; i <= ncol; i++) {
        total++
        if ($i == "Y") done++
    }
    status = (done == total ? "DONE" : (done == 0 ? "....." : "..."))
    printf "  %-45s %2d/%-2d  %s\n", $1, done, total, status
}' "$OUT"

echo
echo "Quick views:"
echo "  Runs fully done:        awk -F'\\t' 'NR==1 || \$NF==\"Y\"' $OUT | column -t -s\$'\\t' | less -S"
echo "  Runs still pending:     awk -F'\\t' 'NR==1 || \$NF==\"N\"' $OUT | column -t -s\$'\\t' | less -S"
echo "  Where 50kb is missing:  awk -F'\\t' 'NR==1 || \$18==\"N\"' $OUT | column -t -s\$'\\t' | less -S"
echo "  Open in spreadsheet:    scp \$USER@cosmos.lunarc.lu.se:.../LepEU/$OUT ."
