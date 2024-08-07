#!/bin/bash
#SBATCH --job-name="DEXseq_lung"
#SBATCH --error=error/dexseq_test1.error
#SBATCH --output=Dexseq_lung.out
#SBATCH --mail-user=n.dubey1@universityofgalway.ie
#SBATCH --mail-type=ALL
#SBATCH --time=5-16:00:00
#SBATCH -p normal
#SBATCH -c 8

cd /data4/msc23104469/dexseq/dexseq_lung

#activate conda environment with R and essential packages installed
module load Anaconda3
conda activate dexseq

Rscript --vanilla dexseq_lung.R
