#!/bin/bash
#SBATCH --job-name=bam_depth
#SBATCH --array=1-1257%100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=logs/bam_depth_%A_%a.out
#SBATCH --error=logs/bam_depth_%A_%a.err
#
# c4 SLURM wrapper around bam_depth.sh. One array task = one line of the BAM list.
#
# Usage (from the bam_depth/bam_depth directory, after `mkdir -p logs`):
#   sbatch --array=1-N%100 c4_wrapper.sh ../configs/c4_b38.config bams.txt /path/to/outdir [--on-target-only]
#
# Set the array range to the number of lines in the BAM list (wc -l bams.txt).
# Then run:  bash summarize_depth.sh /path/to/outdir

set -euo pipefail

CONFIG="${1:?config file required}"
BAM_LIST="${2:?BAM list required}"
OUTDIR="${3:?output dir required}"
shift 3

# sbatch runs a spooled copy of this script, so locate the repo from the submit dir
# (or set BAM_DEPTH_DIR explicitly).
SCRIPT_DIR="${BAM_DEPTH_DIR:-${SLURM_SUBMIT_DIR:?run via sbatch from the bam_depth/bam_depth directory, or set BAM_DEPTH_DIR}}"

: "${SLURM_ARRAY_TASK_ID:?this wrapper must be submitted as an array job}"

module load bedtools2 samtools
# mosdepth must be on PATH (module or conda env)

bash "$SCRIPT_DIR/bam_depth.sh" \
    --config "$CONFIG" \
    --bam-list "$BAM_LIST" \
    --outdir "$OUTDIR" \
    --line "$SLURM_ARRAY_TASK_ID" \
    --threads "${SLURM_CPUS_PER_TASK:-2}" \
    "$@"
