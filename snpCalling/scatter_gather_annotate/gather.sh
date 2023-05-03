#!/usr/bin/env bash

module purge

module load htslib/1.10.2 bcftools/1.9 intel/18.0 intelmpi/18.0 parallel/20200322 vcftools/0.1.16

popSet=${1}
method=${2}
maf=${3}
mac=${4}
version=${5}
wd=${6}
chr=${7}

echo "Chromosome: $chr"

bcf_outdir="${wd}/sub_bcf"
if [ ! -d $bcf_outdir ]; then
    mkdir $bcf_outdir
fi

outdir=$wd/sub_vcfs
cd ${wd}


echo "making list"

ls -d ${outdir}/*.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | \
rev | cut -f1 -d '/' |rev | grep "^${chr}_"| sort -t"_" -k2n,2 -k4g,4 | \
sed "s|^|$outdir|g" > $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort



echo "Concatenating"

vcf-concat \
-f $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort \
-s | \
bgzip -c > $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

echo "tabix'ing"

tabix -p vcf $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
