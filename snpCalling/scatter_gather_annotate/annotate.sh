#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

### #
### #SBATCH -J manual_annotate # A single job name for the array
### #SBATCH --ntasks-per-node=20 # one core
### #SBATCH -N 1 # on one node
### #SBATCH -t 14:00:00 ### 1 hours
### #SBATCH --mem 40G
### #SBATCH -o /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_annotate.%A_%a.out # Standard output
### #SBATCH -e /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_annotate.%A_%a.err # Standard error
### #SBATCH -p standard
### #SBATCH --account biol4559-aob2x
###
### ### cat /scratch/aob2x/DESTv2_output_SNAPE/logs/runSnakemake.49369837*.err
###
### ### sbatch /scratch/aob2x/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/manual_annotate.sh
### ### sacct -j 49432588
### ### cat /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_annotate*.out

module purge

#module load  htslib/1.10.2 bcftools/1.9 intel/18.0 intelmpi/18.0 parallel/20200322 R/3.6.3 samtools vcftools
module load htslib/1.17  bcftools/1.17 parallel/20200322 gcc/11.4.0 openmpi/4.1.4 python/3.11.4 perl/5.36.0 vcftools/0.1.16 samtools/1.17 R/4.3.1
echo "R_LIBS_USER=~/R/goolf/4.3" > ~/.Renviron


popSet=${1}
method=${2}
maf=${3}
mac=${4}
version=${5}
wd=${6}
snpEffPath=${7}

cd ${wd}

echo "concat"

  ls -d ${wd}/sub_bcf/dest.*.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | grep -E "2L|2R|3L|3R|X" > \
  ${wd}/sub_bcf/vcf_order.genome

  bcftools concat \
  -f ${wd}/sub_bcf/vcf_order.genome \
  -O z \
  -o ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz


  tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

echo "convert to vcf & annotate"
  bcftools view \
  --threads 20 \
  ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | \
  java -jar ~/snpEff/snpEff.jar \
  eff \
  BDGP6.86 - > \
  ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf

echo "make GDS"
  Rscript --vanilla /home/kjl5t/Bergland/snpCalling_v2/snpCalling/scatter_gather_annotate/vcf2gds.R ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf

echo "bgzip & tabix"
 bgzip -@20 -c ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf > ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz
 tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz
















### old

# echo "concat"
#
#   ls -d ${wd}/sub_bcf/dest.*.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz > \
#   ${wd}/sub_bcf/vcf_order.genome
#
#
#   vcf-concat \
#   -f ${wd}/sub_bcf/vcf_order.genome \
#   -s \
#   |  \
#   bgzip -c > ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#
#   tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#
# echo "convert to vcf & annotate"
#   bcftools view \
#   --threads 10 \
#   ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | \
#   java -jar ${snpEffPath}/snpEff.jar \
#   eff \
#   BDGP6.86 - > \
#   ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf
#
#echo "make GDS"
#   Rscript --vanilla /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.R ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf
#
#echo "bgzip & tabix"
#  bgzip -@20 -c ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf > ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz
#  tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz






# ### old
#
# #!/usr/bin/env bash
#
# module purge
#
# module load  htslib/1.10.2 bcftools/1.9 intel/18.0 intelmpi/18.0 parallel/20200322 R/3.6.3
#
#
# popSet=${1}
# method=${2}
# maf=${3}
# mac=${4}
# version=${5}
# wd=${6}
# snpEffPath=${7}
#
# cd ${wd}
#
# echo "index"
#   bcftools index -f ${wd}/sub_bcf/dest.2L.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f ${wd}/sub_bcf/dest.2R.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f ${wd}/sub_bcf/dest.3L.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f ${wd}/sub_bcf/dest.3R.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f  ${wd}/sub_bcf/dest.X.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f  ${wd}/sub_bcf/dest.Y.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f  ${wd}/sub_bcf/dest.4.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#   bcftools index -f ${wd}/sub_bcf/dest.mitochondrion_genome.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#
#
# echo "concat"
#   bcftools concat \
#   -n \
#   -O z \
#   ${wd}/sub_bcf/dest.*.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz \
#   -o ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
#
# echo "convert to vcf & annotate"
#   bcftools view \
#   --threads 10 \
#   ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | \
#   java -jar ${snpEffPath}/snpEff.jar \
#   eff \
#   BDGP6.86 - > \
#   ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf
#
#
# echo "fix header" #this is now fixed in PoolSNP.py
# #  sed -i '0,/CHROM/{s/AF,Number=1/AF,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
# #  sed -i '0,/CHROM/{s/AC,Number=1/AC,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
# #  sed -i '0,/CHROM/{s/AD,Number=1/AD,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
# #  sed -i '0,/CHROM/{s/FREQ,Number=1/FREQ,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
#
# #  bcftools view -h ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf > ${wd}/tmp.header
# #
# #  bcftools reheader --threads 10 -h ${wd}/tmp.header -o ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.header.bcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.bcf
#
# echo "make GDS"
#   Rscript --vanilla ${wd}/../DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.R ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf
#
# echo "bgzip & tabix"
#   bgzip -c ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf > ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz
#   tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf.gz
#
#   rm ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.ann.vcf
#   rm ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf
#
