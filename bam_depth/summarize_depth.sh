#!/bin/bash
#
# summarize_depth.sh - compile per-sample mosdepth results into one table.
#
# Usage: summarize_depth.sh OUTDIR
#   Reads   OUTDIR/tmp_depth/*_on.mosdepth.summary.txt (and *_off... if present)
#   Writes  OUTDIR/depth_summary.tsv  with columns: sample_id, mean_genome, mean_on, mean_off
#   mean_off is NA for samples without an off-target result (e.g. --on-target-only runs).
#   If any OUTDIR/tmp_depth/*.gene_depth.tsv.gz exist (runs with --gene-bed), also writes
#   OUTDIR/gene_depth.tsv (long format): sample_id, gene, ensg, n_bases, mean_depth
#   [, frac_10x, frac_20x, ... if thresholds were used].

set -euo pipefail

OUTDIR="${1:?Usage: summarize_depth.sh OUTDIR}"
RESULTS_DIR="$OUTDIR/tmp_depth"
[[ -d $RESULTS_DIR ]] || { echo "Error: $RESULTS_DIR not found" >&2; exit 1; }

summary_file="$OUTDIR/depth_summary.tsv"
printf 'sample_id\tmean_genome\tmean_on\tmean_off\n' > "$summary_file"

shopt -s nullglob
on_files=("$RESULTS_DIR"/*_on.mosdepth.summary.txt)
(( ${#on_files[@]} > 0 )) || { echo "Error: no *_on.mosdepth.summary.txt files in $RESULTS_DIR" >&2; exit 1; }

n=0; n_off=0
for on_file in "${on_files[@]}"; do
    base="$(basename "$on_file")"
    sample="${base%_on.mosdepth.summary.txt}"
    off_file="$RESULTS_DIR/${sample}_off.mosdepth.summary.txt"

    # 'total' = whole-genome mean; 'total_region' = mean over the BED regions (4th column)
    mean_genome="$(awk '$1=="total"        {print $4}' "$on_file")"
    mean_on="$(awk     '$1=="total_region" {print $4}' "$on_file")"
    mean_off="NA"
    if [[ -f $off_file ]]; then
        mean_off="$(awk '$1=="total_region" {print $4}' "$off_file")"
        n_off=$((n_off + 1))
    fi

    printf '%s\t%s\t%s\t%s\n' "$sample" "${mean_genome:-NA}" "${mean_on:-NA}" "${mean_off:-NA}" >> "$summary_file"
    n=$((n + 1))
done

echo "Summary written: $summary_file"
echo "Samples: $n (with off-target: $n_off)"

# ---- per-gene table (only if bam_depth.sh was run with --gene-bed) ----
gene_files=("$RESULTS_DIR"/*.gene_depth.tsv.gz)
if (( ${#gene_files[@]} > 0 )); then
    gene_file="$OUTDIR/gene_depth.tsv"
    { printf 'sample_id\t'; zcat "${gene_files[0]}" | head -n 1; } > "$gene_file"
    for f in "${gene_files[@]}"; do
        base="$(basename "$f")"
        zcat "$f" | awk -F'\t' -v OFS='\t' -v s="${base%.gene_depth.tsv.gz}" 'NR > 1 {print s, $0}' >> "$gene_file"
    done
    echo "Gene table written: $gene_file (samples: ${#gene_files[@]})"
fi
