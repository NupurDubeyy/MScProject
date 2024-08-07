#!/bin/bash
#SBATCH --time=3-16:00:00
#SBATCH --job-name="fastqc_multiqc"
#SBATCH --output=fastqc_multiqc_highmem.out
#SBATCH --mail-user=n.dubey1@universityofgalway.ie
#SBATCH --mail-type=ALL
#SBATCH --mem-per-cpu=120G   # memory per cpu-core
#SBATCH -p highmem

module load Anaconda3
conda activate qc

# date written = 06/05/2024 

start_time=$SECONDS

cd /data4/msc23104469/
echo 'current directory...'
pwd
#mkdir outdir
files="/home/msc23104469/sample_sheet.tsv"
echo "my files"
cat $files
# run fastqc on all samples
while IFS=$'\t' read -r -a myArray
do
 echo "reading file..."
 echo "${myArray[0]}"
 var=$(date)
 echo "$var"
 echo "reading reverse read..."
 echo "${myArray[1]}"
 fastqc "${myArray[0]}" "${myArray[1]}" -o /data4/msc23104469/caf-breast/multiqc_rerun/
done < $files

#Run MultiQC on the FastQC output files
cd /data4/msc23104469/caf-breast/multiqc_rerun/
echo "multiqc directory..."
pwd
multiqc . -o .

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
