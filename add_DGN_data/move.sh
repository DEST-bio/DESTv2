#!/usr/bin/env bash
#
#SBATCH -J move # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:20:00 ### 0.5 hours
##SBATCH --mem 3G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/move.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol8083

module load parallel htslib

wd="/scratch/aob2x/dest"

# rm -R -v !("/project/berglandlab/DEST/dest_mapped/example_pipeline_output")

#SLURM_ARRAY_TASK_ID=1
pop=$( cat ${wd}/dgn/pops.delim | cut -f3 | sort | uniq | sed "${SLURM_ARRAY_TASK_ID}q;d" )

echo ${pop}



[ ! -d /project/berglandlab/DEST/dest_mapped/${pop}/ ] && mkdir /project/berglandlab/DEST/dest_mapped/${pop}

tabix -f -b 2 -s 1 -e 2 ${wd}/dest/wholeGenomeSyncData/${pop}.sync.gz

mv ${wd}/dest/wholeGenomeSyncData/${pop}.sync.gz ${wd}/dest/wholeGenomeSyncData/${pop}.masked.sync.gz
mv ${wd}/dest/wholeGenomeSyncData/${pop}.sync.gz.tbi ${wd}/dest/wholeGenomeSyncData/${pop}.masked.sync.gz.tbi

rsync ${wd}/dest/wholeGenomeSyncData/${pop}.masked.sync.gz* /project/berglandlab/DEST/dest_mapped/${pop}/


echo $pop
