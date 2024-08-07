#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="star_index_primary"
#SBATCH --output=star_index_primary.out
#SBATCH --mail-user=n.dubey1@universityofgalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=60G   # memory per cpu-core
#SBATCH -p highmem

module load Anaconda3
source activate alignment

# date written = 31/05/2024 

start_time=$SECONDS

cd /data4/msc23104469/star/
echo "current directory..."
pwd
mkdir -p index_star_primary

STAR --runMode genomeGenerate --genomeDir /data4/msc23104469/star/index_star_primary --genomeFastaFiles /data4/msc23104469/star/primary_sequence/Homo_sapiens.GRCh38.dna.primary_assembly.fa --runThreadN 8

# Calculate and print elapsed time
elapsed=$(( SECONDS - start_time ))
echo "Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"

