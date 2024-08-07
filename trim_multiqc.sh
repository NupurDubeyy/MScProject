#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="trim_multiqc"
#SBATCH --output=trim_multiqc.out
#SBATCH --mail-user=n.dubey1@universityofgalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=120G   # memory per cpu-core
#SBATCH -p highmem

module load Anaconda3
conda activate qc

# Date written = 08/05/2024

start_time=$SECONDS

cd /data4/msc23104469/caf-breast/caf-breast/trimmed_reads/multiqc_trimmed_reads/
#echo "current directory..."
#pwd
#files="/home/msc23104469/sample_sheet.tsv"
#cutadapt_path="/data4/msc23104469/ .conda/envs/qc/bin/cutadapt"
#echo "my files"
#cat $files

# Trim reads using Trim Galore
#echo "Trimming reads..."
#mkdir -p trimmed_reads

#while IFS=$'\t' read -r -a myArray
#do
# echo "reading file..."
# echo "${myArray[0]}"
# var=$(date)
# echo "$var"
# echo "reading reverse read..."
# echo "${myArray[1]}"
# trim_galore --fastqc "${myArray[0]}" "${myArray[1]}" --cores 4 --paired --gzip --output_dir trimmed_reads
#done < $files

# Generate MultiQC report on trimmed reads
echo "Generating MultiQC report..."
#files="/data4/msc23104469/caf-lung/outdir_lung/trimmed_reads/multiqc_trimmed_reads"
#echo "my files"
#cat $files

# Define directories
TRIMMED_READS_DIR="/data4/msc23104469/caf-breast/caf-breast/trimmed_reads/multiqc_trimmed_reads"
OUTPUT_DIR="$TRIMMED_READS_DIR/multiqc_report/"

#cd /data4/msc23104469/caf-lung/outdir_lung/trimmed_reads
echo "multiqc directory..."
pwd
multiqc "$TRIMMED_READS_DIR" -o "$OUTPUT_DIR" -n "Trimmed Reads Report"
#done < $files


# Calculate and print elapsed time
elapsed=$(( SECONDS - start_time ))
echo "Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
