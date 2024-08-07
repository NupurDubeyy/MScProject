#!/bin/bash
#SBATCH --job-name="DEXseq_breast"
#SBATCH --error=error/%A_%a.error
#SBATCH --output=Dexseq_breast.out
#SBATCH --mail-user=n.dubey1@universityofgalway.ie
#SBATCH --mail-type=ALL
#SBATCH --time=5-16:00:00
#SBATCH --mem-per-cpu=60G   # memory per cpu-core
#SBATCH -p highmem
#SBATCH --cpus-per-task=8


cd /data4/msc23104469/dexseq/dexseq_breast

#activate conda environment with R and essential packages installed
module load Anaconda3
conda activate dexseq

Rscript --vanilla dexseq_breast.R
