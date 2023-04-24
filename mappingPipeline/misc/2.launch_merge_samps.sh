#!/usr/bin/env bash
#
#SBATCH -J blsmm_env # A single job name for the array
#SBATCH --ntasks-per-node=5 # ten cores
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 ### 20 minutes per job per 10 permutations
#SBATCH --mem 40G
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-16

module load gcc/7.1.0 openmpi/3.1.4 R/4.1.1 gdal proj

Rscript \
--vanilla \
2.Merge_samps_to_Collapse.R \
${SLURM_ARRAY_TASK_ID}


