#!/bin/bash
#SBATCH --job-name="lung_rMATS"
#SBATCH --output=logs/rmats_run.out
#SBATCH --error=error/error.error
#SBATCH -p highmem # either normal or highmem - if you need >60GB ram, go for highmem
#SBATCH --mail-user=n.dubey1@universityofgalway.ie
#SBATCH --mail-type=ALL
#SBATCH --time=5-16:00:00
#SBATCH --mem-per-cpu=120gb # memory per cpu core


start_time=$SECONDS



indir="/data4/msc23104469/rMATS_files/lung_rmats"
cd $indir
mkdir -p ${indir}/rmats_outdir
# activate conda environment in which rMATS is installed
module load singularity/3.4.1
singularity shell /data/containers/rmats_latest.sif

RUN=${SLURM_ARRAY_TASK_ID:-1}

echo "Run: ${RUN}"

singularity run -B /data4/msc23104469/rMATS_files /data/containers/rmats_latest.sif python /rmats/rmats.py --s1 /data4/msc23104469/rMATS_files/lung_rmats/s1.txt --s2 /data4/msc23104469/rMATS_files/lung_rmats/s2.txt --gtf /data4/msc23104469/rMATS_files/Homo_sapiens.GRCh38.112.gtf --bi /data4/msc23104469/rMATS_files/STAR_index -t paired --libType fr-firststrand --readLength 150 --nthread 32 --paired-stats --novelSS --individual-counts --od /data4/msc23104469/rMATS_files/lung_rmats/rmats_outdir --tmp /data4/msc23104469/rMATS_files/lung_rmats/rmats_outdir

elapsed=$(( SECONDS - start_time ))
eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
