#!/usr/bin/env bash
#
#SBATCH -J wide2long # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:30:00 ### 30 minutes
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/wide2long.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/wide2long.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


####################
### wide to long ###
####################

###
  #SLURM_ARRAY_TASK_ID=1

  filestem=$( grep -w "^${SLURM_ARRAY_TASK_ID}" /scratch/aob2x/dest/dgn/dgn_wideFiles.delim | cut -f2 )

  mkdir -p /scratch/aob2x/dest/dgn/longData




  sed 's/\(.\)/\1\n/g' \
  /scratch/aob2x/dest/dgn/wideData/${filestem} > \
  /scratch/aob2x/dest/dgn/longData/${filestem}.long
