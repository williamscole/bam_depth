#!/bin/bash
#
# bam_depth.sh - on-/off-target read depth with mosdepth.
#
# Usage:
#   bam_depth.sh --config CONFIG --bam-list LIST --outdir DIR [--line N] [--on-target-only] [--gene-bed BED [--gene-thresholds LIST]] [--threads T]
#
#   --config          bash config file (BUILD, EXOME_TARGET_BED)
#   --bam-list        text file, one BAM path per line
#   --outdir          results go to DIR/tmp_depth/ (summarize_depth.sh writes DIR/depth_summary.tsv)
#   --line N          only process line N (1-based) of the BAM list (for array jobs).
#                     If omitted, every BAM in the list is processed in turn.
#   --on-target-only  skip the off-target BED and off-target mosdepth run
#   --gene-bed BED    also compute mean depth per gene. BED (.bed or .bed.gz) with one interval
#                     per exon/CDS: chrom, start, end, gene [, ... , 7th col = Ensembl gene ID].
#                     Writes DIR/tmp_depth/<sample>.gene_depth.tsv.gz. Independent of --on-target-only.
#   --gene-thresholds LIST  with --gene-bed: comma-separated depths N for which to also report the
#                     fraction of each gene's bases covered at >= Nx (default 10,20,30; 'none' = skip).
#   --threads T       mosdepth threads (default: $SLURM_CPUS_PER_TASK or 2)

set -euo pipefail

usage() { sed -n '3,20p' "$0" | sed 's/^# \{0,1\}//'; }
die()   { echo "Error: $*" >&2; exit 1; }

CONFIG=""; BAM_LIST=""; OUTDIR=""; LINE=""; ON_ONLY="false"; GENE_BED=""; GENE_THR="10,20,30"
THREADS="${SLURM_CPUS_PER_TASK:-2}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)         CONFIG="${2:-}";   shift 2 ;;
        --bam-list)       BAM_LIST="${2:-}"; shift 2 ;;
        --outdir)         OUTDIR="${2:-}";   shift 2 ;;
        --line)           LINE="${2:-}";     shift 2 ;;
        --threads)        THREADS="${2:-}";  shift 2 ;;
        --on-target-only) ON_ONLY="true";    shift ;;
        --gene-bed)       GENE_BED="${2:-}"; shift 2 ;;
        --gene-thresholds) GENE_THR="${2:-}"; shift 2 ;;
        -h|--help)        usage; exit 0 ;;
        *)                usage; die "unknown argument: $1" ;;
    esac
done

[[ -n $CONFIG && -n $BAM_LIST && -n $OUTDIR ]] || { usage; die "--config, --bam-list and --outdir are required"; }
[[ -f $CONFIG ]]   || die "config not found: $CONFIG"
[[ -f $BAM_LIST ]] || die "BAM list not found: $BAM_LIST"
[[ -z $GENE_BED || -f $GENE_BED ]] || die "gene BED not found: $GENE_BED"
[[ $GENE_THR == none || $GENE_THR =~ ^[0-9]+(,[0-9]+)*$ ]] || die "--gene-thresholds must be comma-separated integers (e.g. 10,20,30) or 'none', got '$GENE_THR'"
if [[ -n $LINE ]]; then
    [[ $LINE =~ ^[1-9][0-9]*$ ]] || die "--line must be a positive integer, got '$LINE'"
fi

# ---- config ----
# shellcheck disable=SC1090
source "$CONFIG"
: "${BUILD:?BUILD not set in $CONFIG}"
: "${EXOME_TARGET_BED:?EXOME_TARGET_BED not set in $CONFIG}"
[[ -f $EXOME_TARGET_BED ]] || die "EXOME_TARGET_BED not found: $EXOME_TARGET_BED"
# STRIP_CHR / ADD_CHR / REF_FASTA from older configs are no longer used: contig names and lengths
# are taken from each BAM's own header (see prepare_regions).

command -v mosdepth >/dev/null || die "mosdepth not on PATH"
command -v bedtools >/dev/null || die "bedtools not on PATH"
command -v samtools >/dev/null || die "samtools not on PATH"

RESULTS_DIR="$OUTDIR/tmp_depth"
mkdir -p "$RESULTS_DIR"

# Per-run scratch space (node-local $TMPDIR if set), removed on exit
WORK="$(mktemp -d "${TMPDIR:-/tmp}/bam_depth.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

export LC_ALL=C

# ---- region files in a canonical contig style (no "chr" prefix, mitochondrion = MT) ----
# Per BAM, these are renamed to that BAM's style, so a cohort can mix "1" and "chr1" BAMs.
canon_contigs() { awk 'BEGIN{FS=OFS="\t"} {sub(/^chr/, "", $1); if ($1 == "M") $1 = "MT"; print}'; }
read_bed()      { if [[ $1 == *.gz ]]; then zcat "$1"; else cat "$1"; fi; }

read_bed "$EXOME_TARGET_BED" | awk -F'\t' -v OFS='\t' '/^(#|track|browser)/ || NF < 3 {next} {print $1, $2, $3}' \
    | canon_contigs | sort -k1,1 -k2,2n > "$WORK/on_canon.bed"
[[ -s $WORK/on_canon.bed ]] || die "no intervals read from EXOME_TARGET_BED: $EXOME_TARGET_BED"

# gene intervals: 4-column BED (chrom, start, end, "gene|ENSG")
if [[ -n $GENE_BED ]]; then
    read_bed "$GENE_BED" \
        | awk -F'\t' -v OFS='\t' '
            /^(#|track|browser)/ || NF == 0 {next}
            NF < 4 {bad = 1; exit}
            {key = $4; if (NF >= 7 && $7 != "") key = key "|" $7; print $1, $2, $3, key}
            END {if (bad) exit 3}' \
        | canon_contigs | sort -k1,1 -k2,2n > "$WORK/gene_canon.bed" \
        || die "gene BED needs at least 4 tab-separated columns (chrom, start, end, gene): $GENE_BED"
    [[ -s $WORK/gene_canon.bed ]] || die "no intervals read from gene BED: $GENE_BED"
fi

# Rename canonical BED $1 to this BAM's contig style, keep only contigs present in the BAM header, write $2.
# Needs BAM_STYLE and $WORK/bam_contigs.txt (set by prepare_regions).
to_bam_style() {
    awk -F'\t' -v OFS='\t' -v style="$BAM_STYLE" '
        NR == FNR {ok[$1] = 1; next}
        {c = $1; if (style == "chr") c = (c == "MT" ? "chrM" : "chr" c)
         if (c in ok) {$1 = c; print}}' "$WORK/bam_contigs.txt" "$1" | sort -k1,1 -k2,2n > "$2"
}

# Per BAM: read contig names/lengths from the header, detect "chr" vs no-"chr" naming, and write
# on_target.bed [, off_target.bed] [, gene.bed] in that BAM's naming.
prepare_regions() {
    local bam="$1" label="$2" n_in n_out
    samtools view -H "$bam" | awk -F'\t' -v OFS='\t' '
        $1 == "@SQ" {n = ""; l = ""
            for (i = 2; i <= NF; i++) {if ($i ~ /^SN:/) n = substr($i, 4); else if ($i ~ /^LN:/) l = substr($i, 4)}
            if (n != "" && l != "") print n, l}' | sort -k1,1 > "$WORK/bam_contigs.txt"
    [[ -s $WORK/bam_contigs.txt ]] || die "$label: no @SQ lines in the BAM header (is it a BAM with a header?)"
    if   awk '$1 == "chr1" {f = 1} END {exit !f}' "$WORK/bam_contigs.txt"; then BAM_STYLE="chr"
    elif awk '$1 == "1"    {f = 1} END {exit !f}' "$WORK/bam_contigs.txt"; then BAM_STYLE="plain"
    else die "$label: cannot tell the contig naming style; expected 1,2,... or chr1,chr2,... but the header starts with: $(head -n 3 "$WORK/bam_contigs.txt" | cut -f1 | paste -sd' ')"
    fi

    to_bam_style "$WORK/on_canon.bed" "$WORK/on_target.bed"
    n_in=$(wc -l < "$WORK/on_canon.bed"); n_out=$(wc -l < "$WORK/on_target.bed")
    [[ $n_out -gt 0 ]] || die "$label: none of the target regions are on contigs in the BAM header (naming: $BAM_STYLE)"
    (( n_out == n_in )) || echo "Warning: $label: $((n_in - n_out)) of $n_in target intervals are on contigs not in the BAM header (skipped)" >&2

    if [[ $ON_ONLY != "true" ]]; then
        bedtools complement -i "$WORK/on_target.bed" -g "$WORK/bam_contigs.txt" > "$WORK/off_target.bed"
    fi

    if [[ -n $GENE_BED ]]; then
        to_bam_style "$WORK/gene_canon.bed" "$WORK/gene.bed"
        n_in=$(wc -l < "$WORK/gene_canon.bed"); n_out=$(wc -l < "$WORK/gene.bed")
        [[ $n_out -gt 0 ]] || die "$label: none of the gene regions are on contigs in the BAM header (naming: $BAM_STYLE)"
        (( n_out == n_in )) || echo "Warning: $label: $((n_in - n_out)) of $n_in gene intervals are on contigs not in the BAM header (skipped)" >&2
    fi
}

# ---- per-BAM work ----
process_bam() {
    local bam="$1" sample
    [[ -f $bam ]] || die "BAM not found: $bam"
    sample="$(basename "$bam" .bam)"
    echo "[$(date +%T)] Processing $sample ($bam)"
    local t0=$SECONDS
    prepare_regions "$bam" "$sample"
    echo "  contig naming: $BAM_STYLE"

    # -n: skip per-base output (not used downstream; much less I/O)
    mosdepth -n -t "$THREADS" -b "$WORK/on_target.bed" "$WORK/${sample}_on" "$bam"
    cp "$WORK/${sample}_on.mosdepth.summary.txt" "$WORK/${sample}_on.mosdepth.global.dist.txt" "$RESULTS_DIR/"

    if [[ $ON_ONLY != "true" ]]; then
        mosdepth -n -t "$THREADS" -b "$WORK/off_target.bed" "$WORK/${sample}_off" "$bam"
        cp "$WORK/${sample}_off.mosdepth.summary.txt" "$WORK/${sample}_off.mosdepth.global.dist.txt" "$RESULTS_DIR/"
    fi

    if [[ -n $GENE_BED ]]; then
        local thr_args=()
        if [[ $GENE_THR != none ]]; then
            thr_args=(--thresholds "$GENE_THR")
        else
            printf '#chrom\tstart\tend\tregion\n' | gzip > "$WORK/${sample}_gene.thresholds.bed.gz"
        fi
        mosdepth -n -t "$THREADS" -b "$WORK/gene.bed" "${thr_args[@]}" "$WORK/${sample}_gene" "$bam"
        # per-gene results, weighted by interval length:
        #   mean_depth = sum(mean_i * len_i) / sum(len_i)      (mosdepth rounds mean_i to 2 decimals)
        #   frac_Nx    = sum(bases_i at >= Nx) / sum(len_i)
        # input 1: thresholds.bed.gz (header + bases >= N per interval); input 2: regions.bed.gz (mean per interval)
        awk -F'\t' -v OFS='\t' '
            FNR == 1 {f++}
            f == 1 {
                if ($0 ~ /^#/) {nt = NF - 4; for (i = 5; i <= NF; i++) tn[i - 4] = tolower($i); next}
                k = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
                for (i = 1; i <= nt; i++) thr[k, i] = $(i + 4)
                next
            }
            {
                len = $3 - $2; if (len <= 0) next
                g = $4; bases[g] += len; cov[g] += $5 * len
                k = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
                for (i = 1; i <= nt; i++) tb[g, i] += thr[k, i]
            }
            END {
                h = "gene\tensg\tn_bases\tmean_depth"
                for (i = 1; i <= nt; i++) h = h "\tfrac_" tn[i]
                print h
                for (g in bases) {
                    split(g, a, "|"); ensg = (a[2] == "" ? "NA" : a[2])
                    line = sprintf("%s\t%s\t%d\t%.4f", a[1], ensg, bases[g], cov[g] / bases[g])
                    for (i = 1; i <= nt; i++) line = line sprintf("\t%.4f", tb[g, i] / bases[g])
                    print line
                }
            }' <(zcat "$WORK/${sample}_gene.thresholds.bed.gz") <(zcat "$WORK/${sample}_gene.regions.bed.gz") \
            | { IFS= read -r header; echo "$header"; sort -k1,1 -k2,2; } \
            | gzip > "$RESULTS_DIR/${sample}.gene_depth.tsv.gz"
    fi

    rm -f "$WORK/${sample}_"*
    echo "[$(date +%T)] Done $sample in $((SECONDS - t0))s"
}

if [[ -n $LINE ]]; then
    bam="$(sed -n "${LINE}p" "$BAM_LIST")"
    [[ -n $bam ]] || die "line $LINE of $BAM_LIST is empty or past the end of the file"
    process_bam "$bam"
else
    while IFS= read -r bam || [[ -n $bam ]]; do
        [[ -z $bam || $bam == \#* ]] && continue
        process_bam "$bam"
    done < "$BAM_LIST"
fi
