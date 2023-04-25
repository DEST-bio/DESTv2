#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 7
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00
#SBATCH --mem 60G
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/RunDest.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-16

pwd
echo $SLURM_CPUS_PER_TASK 

### modules
module load singularity

###################################
# Part  0. Get Sample information #
###################################

sampleid=$( cat mapping.data.only.collapsed.txt  | sed '1d' | awk '{ print $1 }' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

numFlies=$( cat mapping.data.only.collapsed.txt  | sed '1d' | awk '{ print $2 }' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

echo $sampleid
echo $numFlies

###################################
# Part  2. Run Docker             #
###################################

singularity run \
/scratch/yey2sn/DEST/dest_freeze1_latest.sif \
${sampleid}/${sampleid}.joint_1.fastq \
${sampleid}/${sampleid}.joint_2.fastq \
${sampleid} \
/project/berglandlab/DEST/dest_mapped/RECENT_OUTPUTS \
--cores $SLURM_CPUS_PER_TASK \
--max-cov 0.95 \
--min-cov 4 \
--base-quality-threshold 25 \
--num-flies ${numFlies} \
--do_poolsnp \
--do_snape

echo "done"
date
