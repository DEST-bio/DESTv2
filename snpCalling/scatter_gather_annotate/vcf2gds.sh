#!/usr/bin/env bash
#
#SBATCH -J vcf2gds # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 4:00:00  ### 4 hours
#SBATCH --mem 24G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/vcf2gds.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/vcf2gds.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

### sbatch /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.sh
### sacct -j 5675053
### cat /scratch/aob2x/dest/slurmOutput/vcf2gds.5675053*

module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1

Rscript --vanilla /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.R \
/project/berglandlab/DEST/vcf/dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz \
/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.3May2024.ann.gds



### sneaking a index job here too
# module load samtools/1.17
# tabix -p vcf /project/berglandlab/DEST/vcf/dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz
