#!/bin/bash
#
# bam_depth.sh - on-/off-target read depth with mosdepth.
#
# Usage:
#   bam_depth.sh --config CONFIG --bam-list LIST --outdir DIR [--line N] [--on-target-only] [--gene-bed BED [--gene-thresholds LIST]] [--threads T]
#
#   --config          bash config file (REF_FASTA, EXOME_TARGET_BED, STRIP_CHR/ADD_CHR, BUILD)
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
STRIP_CHR="false"; ADD_CHR="false"
# shellcheck disable=SC1090
source "$CONFIG"
: "${BUILD:?BUILD not set in $CONFIG}"
: "${REF_FASTA:?REF_FASTA not set in $CONFIG}"
: "${EXOME_TARGET_BED:?EXOME_TARGET_BED not set in $CONFIG}"
[[ -f $REF_FASTA ]]        || die "REF_FASTA not found: $REF_FASTA"
[[ -f ${REF_FASTA}.fai ]]  || die "index not found: ${REF_FASTA}.fai (run: samtools faidx $REF_FASTA)"
[[ -f $EXOME_TARGET_BED ]] || die "EXOME_TARGET_BED not found: $EXOME_TARGET_BED"
[[ $STRIP_CHR == "true" && $ADD_CHR == "true" ]] && die "STRIP_CHR and ADD_CHR cannot both be true (check $CONFIG)"

command -v mosdepth >/dev/null || die "mosdepth not on PATH"
command -v bedtools >/dev/null || die "bedtools not on PATH"

RESULTS_DIR="$OUTDIR/tmp_depth"
mkdir -p "$RESULTS_DIR"

# Per-run scratch space (node-local $TMPDIR if set), removed on exit
WORK="$(mktemp -d "${TMPDIR:-/tmp}/bam_depth.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

# ---- build regions (rebuilt every run on purpose, so parallel tasks never clash) ----
export LC_ALL=C

cut -f1,2 "${REF_FASTA}.fai" | sort -k1,1 > "$WORK/genome.txt"
cut -f1 "$WORK/genome.txt" > "$WORK/valid_chroms.txt"

if [[ $EXOME_TARGET_BED == *.gz ]]; then
    zcat "$EXOME_TARGET_BED" | cut -f1-3 > "$WORK/on_raw.bed"
else
    cut -f1-3 "$EXOME_TARGET_BED" > "$WORK/on_raw.bed"
fi

if [[ $STRIP_CHR == "true" ]]; then
    sed 's/^chr//' "$WORK/on_raw.bed" > "$WORK/on_nochr.bed"
    mv "$WORK/on_nochr.bed" "$WORK/on_raw.bed"
fi

if [[ $ADD_CHR == "true" ]]; then
    awk 'BEGIN{OFS="\t"} $1 !~ /^chr/ {$1="chr"$1} {print}' "$WORK/on_raw.bed" > "$WORK/on_chr.bed"
    mv "$WORK/on_chr.bed" "$WORK/on_raw.bed"
fi

# keep only contigs present in the reference, then sort
awk 'NR==FNR{ok[$1]=1; next} ($1 in ok)' "$WORK/valid_chroms.txt" "$WORK/on_raw.bed" \
    | sort -k1,1 -k2,2n > "$WORK/on_target.bed"

[[ -s $WORK/on_target.bed ]] || die "no target regions left after filtering to reference contigs (check STRIP_CHR / ADD_CHR / REF_FASTA contig names in $CONFIG)"

# gene intervals: 4-column BED (chrom, start, end, "gene|ENSG"), same contig handling as the target BED
if [[ -n $GENE_BED ]]; then
    if [[ $GENE_BED == *.gz ]]; then zcat "$GENE_BED"; else cat "$GENE_BED"; fi \
        | awk -F'\t' -v OFS='\t' '
            /^(#|track|browser)/ || NF == 0 {next}
            NF < 4 {bad = 1; exit}
            {key = $4; if (NF >= 7 && $7 != "") key = key "|" $7; print $1, $2, $3, key}
            END {if (bad) exit 3}' > "$WORK/gene_raw.bed" \
        || die "gene BED needs at least 4 tab-separated columns (chrom, start, end, gene): $GENE_BED"
    if [[ $STRIP_CHR == "true" ]]; then
        sed 's/^chr//' "$WORK/gene_raw.bed" > "$WORK/gene_nochr.bed"; mv "$WORK/gene_nochr.bed" "$WORK/gene_raw.bed"
    fi
    if [[ $ADD_CHR == "true" ]]; then
        awk 'BEGIN{OFS="\t"} $1 !~ /^chr/ {$1="chr"$1} {print}' "$WORK/gene_raw.bed" > "$WORK/gene_chr.bed"
        mv "$WORK/gene_chr.bed" "$WORK/gene_raw.bed"
    fi
    awk 'NR==FNR{ok[$1]=1; next} ($1 in ok)' "$WORK/valid_chroms.txt" "$WORK/gene_raw.bed" \
        | sort -k1,1 -k2,2n > "$WORK/gene.bed"
    [[ -s $WORK/gene.bed ]] || die "no gene regions left after filtering to reference contigs (check STRIP_CHR / ADD_CHR in $CONFIG)"
    n_in=$(wc -l < "$WORK/gene_raw.bed"); n_kept=$(wc -l < "$WORK/gene.bed")
    (( n_kept == n_in )) || echo "Warning: dropped $((n_in - n_kept)) of $n_in gene intervals on contigs not in the reference" >&2
fi

if [[ $ON_ONLY != "true" ]]; then
    bedtools complement -i "$WORK/on_target.bed" -g "$WORK/genome.txt" > "$WORK/off_target.bed"
fi

# ---- per-BAM work ----
process_bam() {
    local bam="$1" sample
    [[ -f $bam ]] || die "BAM not found: $bam"
    sample="$(basename "$bam" .bam)"
    echo "[$(date +%T)] Processing $sample ($bam)"
    local t0=$SECONDS

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
