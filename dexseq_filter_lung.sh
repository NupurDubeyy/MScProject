#!/bin/bash
#SBATCH --job-name="DEXseq_filter"
#SBATCH --error=error/dexseq_filter_new.error
#SBATCH --output=Dexseq_filter.out
#SBATCH --mail-user=n.dubey1@universityofgalway.ie
#SBATCH --mail-type=ALL
#SBATCH --time=5-16:00:00
#SBATCH -p highmem
#SBATCH -c 8
#SBATCH --mem-per-cpu=60G 

cd /data4/msc23104469/dexseq/dexseq_filter_l

#activate conda environment with R and essential packages installed
module load Anaconda3
conda activate dexseq

Rscript --vanilla dexseq_filter.R
