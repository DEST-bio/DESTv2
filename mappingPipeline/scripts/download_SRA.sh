#!/usr/bin/env bash
#
#SBATCH -J download_SRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/split_and_run.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

wd=/scratch/aob2x/dest
### grep -E "ES_ba_12|AT_gr_12" /scratch/aob2x/dest/DEST/populationInfo/samps.csv | cut -f1,13 -d',' > /scratch/aob2x/fastq/todl.csv
### run as: sbatch --array=1-$( wc -l < /scratch/aob2x/fastq/todl.csv ) /scratch/aob2x/dest/DEST/mappingPipeline/scripts/download_SRA.sh
### sacct -j 18750292

module load sratoolkit/2.10.5

#SLURM_ARRAY_TASK_ID=1

sranum=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/fastq/todl.csv | cut -f2 -d',' )
sampName=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/fastq/todl.csv | cut -f1 -d',' )

echo $sampName " / " $sranum

prefetch \
-o /scratch/aob2x/fastq/${sranum}.sra \
${sranum}

fasterq-dump \
--split-files \
--split-3 \
--outfile /scratch/aob2x/fastq/${sranum} \
-e 10 \
/scratch/aob2x/fastq/${sranum}.sra

gzip /scratch/aob2x/fastq/${sranum}_1.fastq
gzip /scratch/aob2x/fastq/${sranum}_2.fastq

#rm /scratch/aob2x/fastq/${sranum}.sra
