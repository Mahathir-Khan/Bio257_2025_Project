#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64000mb
#SBATCH --time=16:00:00
#SBATCH --mail-user=mkhan35@ur.rochester.edu
#SBATCH --mail-type=END,FAIL

module load salmon

workdir="/scratch/bio257_2025/Users/group4_ADxTummy/RNAseq_data"
cd "$workdir"


salmon index \
  -t GRCh38_and_decoys.fa.gz \
  -d decoys.txt \
  -p 16 \
  -i GRCh38_salmon_index \
  --gencode

