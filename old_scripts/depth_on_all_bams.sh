#!/bin/bash
#SBATCH --job-name=mosdepth
#SBATCH --array=1-1257%100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=/zivlab/data3/colew/Proj_GLIMPSE/_Output/logs/mosdepth_idx%a_%A.out
#SBATCH --error=/zivlab/data3/colew/Proj_GLIMPSE/_Output/logs/mosdepth_idx%a_%A.err

BUILD=${1}

if [[ $BUILD != "b37" && $BUILD != "b38" ]]; then
    echo "Error: BUILD must be either 'b37' or 'b38'"
    echo "Usage: sbatch mosdepth.sh [b37|b38]"
    exit 1
fi

rootdir="/zivlab/data3"
projdir="$rootdir/colew/Proj_GLIMPSE"
outputdir="$projdir/_Output/depth/${BUILD}/mosdepth"

mkdir -p $outputdir

module load bcftools bedtools2 samtools

if [[ $BUILD == "b38" ]]
then
    ref="$rootdir/colew/datasets/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    exome_target="$rootdir/geri_wes/hg38_exome_v2.0.2_targets_sorted_validated.annotated.bed.gz"
    bamfile_list="$projdir/_Data/b38_geri_bamfiles.txt"
    suffix="_b38.DeDup.bam"
else
    ref="$rootdir/colew/datasets/ref/human_g1k_v37.fasta"
    exome_target="$rootdir/geri_wes/hg19_exome_v2.0.2_merged_probes_sorted_validated.annotated.bed"
    bamfile_list="$projdir/_Data/b37_geri_bamfiles.txt"
    suffix="_b37.DeDup.bam"
fi

# Create genome file from reference index
if [[ $BUILD == "b38" ]]
then
    # b38 ref already has chr prefix
    cut -f1,2 ${ref}.fai | sort -k1,1 > $TMPDIR/genome.txt
else
    # b37 ref doesn't have chr prefix
    cut -f1,2 ${ref}.fai | sort -k1,1 > $TMPDIR/genome.txt
fi

# Create on-target BED file
if [[ $exome_target == *.gz ]]
then
    zcat $exome_target | cut -f1-3 | sort -k1,1 -k2,2n > $TMPDIR/on_target_raw.bed
else
    cut -f1-3 $exome_target | sort -k1,1 -k2,2n > $TMPDIR/on_target_raw.bed
fi

# Remove chr prefix for b37 to match BAM chromosomes
if [[ $BUILD == "b37" ]]
then
    sed 's/^chr//' $TMPDIR/on_target_raw.bed > $TMPDIR/on_target_temp.bed
    mv $TMPDIR/on_target_temp.bed $TMPDIR/on_target_raw.bed
fi

# Filter BED to only include chromosomes in genome file
cut -f1 $TMPDIR/genome.txt > $TMPDIR/valid_chroms.txt
awk 'NR==FNR{chroms[$1]=1; next} $1 in chroms' $TMPDIR/valid_chroms.txt $TMPDIR/on_target_raw.bed > $TMPDIR/on_target.bed

# Create off-target BED file (complement of on-target regions)
bedtools complement -i $TMPDIR/on_target.bed -g $TMPDIR/genome.txt > $TMPDIR/off_target.bed

# Get the BAM file for this array task
bam=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $bamfile_list)

if [[ -z $bam ]]; then
    echo "Error: No BAM file found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

sample_id=$(basename $bam $suffix)

echo "Processing sample: $sample_id"
echo "BAM file: $bam"

start_time=$(date +%s)

mosdepth -t 2 -b $TMPDIR/on_target.bed $TMPDIR/${sample_id}_on $bam
mosdepth -t 2 -b $TMPDIR/off_target.bed $TMPDIR/${sample_id}_off $bam

end_time=$(date +%s)
elapsed_time=$((end_time - start_time))

# Copy output files
cp $TMPDIR/${sample_id}_*.mosdepth.summary.txt $outputdir
cp $TMPDIR/${sample_id}_*.mosdepth.global.dist.txt $outputdir

echo "Completed in ${elapsed_time}s"