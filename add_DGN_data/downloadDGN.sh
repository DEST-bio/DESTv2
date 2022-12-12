#!/usr/bin/env bash
#
#SBATCH -J download_DGN # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/download_DGN.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/download_DGN.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


### download datasets
  #SLURM_ARRAY_TASK_ID=4

  name=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | cut -f2 -d' ' )
  url=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | cut -f3 -d' ' )

  echo "Downloading data"
  echo ${name}

  wget -P /scratch/aob2x/dest/dgn/rawData \
  ${url}
