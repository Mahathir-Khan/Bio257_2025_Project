#!/bin/bash
#SBATCH --job-name=trimming_fastq_for_RNAseq_script
#SBATCH -p standard
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 32000mb
#SBATCH --time=16:00:00
#SBATCH --mail-user=mkhan35@ur.rochester.edu
#SBATCH --mail-type=END,FAIL

# chat do i hardcode this
# i hardcoded it hehe

mkdir /scratch/bio257_2025/Users/group4_ADxTummy/trimmed_fastq_for_RNAseq

module load trimgalore/0.4.4 cutadapt/b1

pair1_list=("SRR14864859_1.fastq"
"SRR14864860_1.fastq" "SRR14864861_1.fastq" "SRR14864862_1.fastq" "SRR14864863_1.fastq" "SRR14864864_1.fastq" "SRR14864$
)

pair2_list=("SRR14864859_2.fastq"
"SRR14864860_2.fastq" "SRR14864861_2.fastq" "SRR14864862_2.fastq" "SRR14864863_2.fastq" "SRR14864864_2.fastq" "SRR14864$
)

input_dir="/scratch/bio257_2025/Users/group4_ADxTummy/untrimmed_fastq_for_RNAseq"
output_dir="/scratch/bio257_2025/Users/group4_ADxTummy/trimmed_fastq_for_RNAseq"

mkdir -p "$output_dir"

for ((i=0; i<12; i++)); do
    read1="${input_dir}/${pair1_list[$i]}"
    read2="${input_dir}/${pair2_list[$i]}"

    trim_galore --paired --phred33 --suppress_warn \
        -o "$output_dir" \
        "$read1" "$read2" || true

        echo "finished pair $((i+1))"

done

echo "finished all 12"




