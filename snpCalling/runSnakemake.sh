#!/bin/bash
#
#SBATCH -J check_fastq_encoding # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 01:00:00 ### 1 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/DESTv2_output/logs/runSnakemake.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DESTv2_output/logs/runSnakemake.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard


module load gcc/9.2.0 openmpi/3.1.6 python/3.7.7 snakemake/6.0.5
cd /scratch/aob2x/DESTv2/snpCalling
snakemake --profile /scratch/aob2x/DESTv2/snpCalling/slurm
