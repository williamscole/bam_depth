BUILD=${1}

if [[ $BUILD != "b37" && $BUILD != "b38" ]]; then
    echo "Error: BUILD must be either 'b37' or 'b38'"
    echo "Usage: bash compile_results.sh [b37|b38]"
    exit 1
fi

rootdir="/zivlab/data3"
projdir="$rootdir/colew/Proj_GLIMPSE"
outputdir="$projdir/_Output/depth/${BUILD}/mosdepth"

# Create output summary file
summary_file="$outputdir/depth_summary.tsv"

# Write header
echo -e "sample_id\tmean_genome\tmean_on\tmean_off" > $summary_file

# Get list of unique sample IDs
sample_ids=$(ls $outputdir/*_on.mosdepth.summary.txt | xargs -n1 basename | sed 's/_on.mosdepth.summary.txt//' | sort -u)

# Process each sample
for sample_id in $sample_ids
do
    on_file="$outputdir/${sample_id}_on.mosdepth.summary.txt"
    off_file="$outputdir/${sample_id}_off.mosdepth.summary.txt"
    
    # Check if both files exist
    if [[ ! -f $on_file ]] || [[ ! -f $off_file ]]; then
        echo "Warning: Missing files for sample $sample_id"
        continue
    fi
    
    # Extract mean depth values (4th column)
    mean_genome=$(grep "^total\s" $on_file | awk '{print $4}')
    mean_on=$(grep "^total_region\s" $on_file | awk '{print $4}')
    mean_off=$(grep "^total_region\s" $off_file | awk '{print $4}')
    
    # Write to summary file
    echo -e "${sample_id}\t${mean_genome}\t${mean_on}\t${mean_off}" >> $summary_file
done

echo "Summary compiled: $summary_file"
echo "Total samples processed: $(tail -n +2 $summary_file | wc -l)"