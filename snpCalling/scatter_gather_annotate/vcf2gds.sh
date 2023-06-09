#!/usr/bin/env bash
#
#SBATCH -J vcf2gds # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:02:00  ### 48 hours
#SBATCH --mem 24G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/vcf2gds.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/vcf2gds.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

### sbatch /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.sh
### sacct -j 50241880
### cat /scratch/aob2x/dest/slurmOutput/vcf2gds.22867938

module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3

Rscript --vanilla /scratch/aob2x/DEST_freeze1/snpCalling/scatter_gather_annotate/vcf2gds.R \
/project/berglandlab/DEST/vcf/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz
