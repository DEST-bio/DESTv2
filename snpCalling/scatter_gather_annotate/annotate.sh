#!/usr/bin/env bash

module purge

module load  htslib/1.10.2 bcftools/1.9 intel/18.0 intelmpi/18.0 parallel/20200322 R/3.6.3


popSet=${1}
method=${2}
maf=${3}
mac=${4}
version=${5}
wd=${6}
snpEffPath=${7}

cd ${wd}

echo "index"
  bcftools index -f ${wd}/sub_bcf/dest.2L.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
  bcftools index -f ${wd}/sub_bcf/dest.2R.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
  bcftools index -f ${wd}/sub_bcf/dest.3L.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
  bcftools index -f ${wd}/sub_bcf/dest.3R.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
  bcftools index -f ${wd}/sub_bcf/dest.X.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
  bcftools index -f ${wd}/sub_bcf/dest.Y.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
  bcftools index -f ${wd}/sub_bcf/dest.4.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
  bcftools index -f ${wd}/sub_bcf/dest.mitochondrion_genome.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz


echo "concat"
  bcftools concat \
  -n \
  -O z \
  ${wd}/sub_bcf/dest.*.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz \
  -o ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz

echo "convert to vcf & annotate"
  bcftools view \
  --threads 10 \
  ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz | \
  java -jar ${snpEffPath}/snpEff.jar \
  eff \
  BDGP6.86 - > \
  ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf


echo "fix header" #this is now fixed in PoolSNP.py
#  sed -i '0,/CHROM/{s/AF,Number=1/AF,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
#  sed -i '0,/CHROM/{s/AC,Number=1/AC,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
#  sed -i '0,/CHROM/{s/AD,Number=1/AD,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
#  sed -i '0,/CHROM/{s/FREQ,Number=1/FREQ,Number=A/}' ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf

#  bcftools view -h ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf > ${wd}/tmp.header
#
#  bcftools reheader --threads 10 -h ${wd}/tmp.header -o ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.header.bcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.bcf

echo "make GDS"
  Rscript --vanilla ${wd}/../DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.R ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf

echo "bgzip & tabix"
  bgzip -c ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf > ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf.gz
  tabix -p vcf ${wd}/dest.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf.gz

rm ${wd}/tmp.header
