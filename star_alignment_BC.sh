#!/bin/bash
#SBATCH --job-name="BC_star"
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=error/%A_%a.error
#SBATCH -p highmem # either normal or highmem - if you need >60GB ram, go for highmem
#SBATCH --mail-user=n.dubey1@universityofgalway.ie
#SBATCH --mail-type=ALL
#SBATCH --time=5-16:00:00
#SBATCH --array=1-24 # initiallise 24 jobs - running star 24 times
#SBATCH --mem-per-cpu=120gb # memory per cpu core


start_time=$SECONDS

indir="/data4/msc23104469/alignment/alignment_BC"
cd $indir
mkdir -p ${indir}/star_outdir_SAM
# activate conda environment in which STAR is installed
module load Anaconda3
conda activate alignment

RUN=${SLURM_ARRAY_TASK_ID:-1}

echo "Run: ${RUN}"
fq1=$(cat /data4/msc23104469/caf-breast/all_fastqs/trimmed_sample_sheet.csv | cut -f 1 -d "," | sed -n ${SLURM_ARRAY_TASK_ID}p)
fq2=$(cat /data4/msc23104469/caf-breast/all_fastqs/trimmed_sample_sheet.csv | cut -f 2 -d "," | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "Read 1: $fq1"
echo "Read 2: $fq2"

file="$(basename -- $fq1)"
echo "Filename: $file"
sample=${file:0:7}
echo "Sample: $sample"

mkdir -p "${indir}/star_outdir_SAM/${sample}"

cd "${indir}/star_outdir_SAM/${sample}"
echo "STAR running directory..."
pwd

echo "prefix: '${indir}/star_outdir_SAM/${sample}/${sample}_'"

STAR --runThreadN 12 --genomeDir /data4/msc23104469/star/index_star_primary/ --sjdbOverhang 149 --sjdbGTFfile /data4/msc23104469/star/Homo_sapiens.GRCh38.112.gtf --readFilesCommand gunzip -c --readFilesIn $fq1 $fq2 --outFileNamePrefix "${indir}/star_outdir_SAM/${sample}/${sample}_"  --outSAMtype SAM --outSAMunmapped Within --outSAMattributes NH HI AS nM NM MD --outSAMstrandField intronMotif

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
