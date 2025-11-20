#!/bin/bash
#SBATCH --job-name=trimming_fastq_for_RNAseq_script
#SBATCH -p standard
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 28
#SBATCH --mem 64000mb
#SBATCH --time=16:00:00
#SBATCH --mail-user=mkhan35@ur.rochester.edu
#SBATCH --mail-type=END,FAIL

module load salmon

# Path to the Salmon index
salmon_index="/scratch/bio257_2025/Users/group4_ADxTummy/RNAseq_data/GRCh38_salmon_index"

# Set the path to the "fastq" folder
fastq_dir="/scratch/bio257_2025/Users/group4_ADxTummy/RNAseq_data/trimmed_fastq_for_RNAseq"

# set output directory
outdir="/scratch/bio257_2025/Users/group4_ADxTummy/RNAseq_data"
mkdir -p "$outdir"

# Loop through all the directories within the "fastq" folder
for dir in "${fastq_dir}"/SRR*; do
    # Find the R1 and R2 FASTQ files
    r1_file=$(find "$dir" -name "*_1_val_1.fq")
    r2_file=$(find "$dir" -name "*_2_val_2.fq")

    # Extract the sample name
    samp=$(basename "$dir")

    salmon quant -i "$salmon_index" -l A \
        -1 "$r1_file" \
        -2 "$r2_file" \
        -p 28 \
        --validateMappings \
        -o "${outdir}/${samp}_quant"
done

echo "finished all 12"
