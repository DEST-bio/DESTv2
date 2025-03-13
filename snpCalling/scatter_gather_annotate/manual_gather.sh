#!/bin/bash
#
#SBATCH -J manual_gather # A single job name for the array
#SBATCH --ntasks-per-node=8 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 14:00:00 ### 1 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_gather.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_gather.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4559-aob2x

### sbatch /scratch/aob2x/CompEvoBio_modules/utils/snpCalling/scatter_gather_annotate/manual_gather.sh
### sacct -j 53544100
### cat /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_gather.53544098.err
### cat /scratch/aob2x/compBio_SNP_25Sept2023/logs/manual_gather
### cd /scratch/aob2x/compBio_SNP_25Sept2023

load gcc/11.4.0  openmpi/4.1.4 python/3.11.4

#module load htslib bcftools parallel intel/18.0 intelmpi/18.0 mvapich2/2.3.1 R/3.6.3 python/3.6.6 vcftools/0.1.16
#module load htslib/1.10.2 bcftools/1.9 parallel/20200322 intel/18.0 intelmpi/18.0 R/3.6.3 python/3.6.6 vcftools/0.1.16
module load htslib/1.17  bcftools/1.17 parallel/20200322 gcc/11.4.0 openmpi/4.1.4 python/3.11.4 perl/5.36.0 vcftools/0.1.16

concatVCF() {

  popSet=PoolSeq
  method=PoolSNP
  maf=001
  mac=50
  version=28Sept2024_ExpEvo
  wd=/scratch/aob2x/compBio_SNP_28Sept2024
  # chr=2L

  chr=${1}

  echo "Chromosome: $chr"

  bcf_outdir="${wd}/sub_bcf"
  if [ ! -d $bcf_outdir ]; then
      mkdir $bcf_outdir
  fi

  outdir=$wd/sub_vcfs
  cd ${wd}

  echo "generate list"
  #ls -d *.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | grep '^${chr}_' | sort -t"_" -k2n,2 -k4g,4 \
  #> $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort


  ls -d ${outdir}/*.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | \
  rev | cut -f1 -d '/' |rev | grep -E "^${chr}_" | sort -t"_" -k2n,2 -k4g,4 | \
  sed "s|^|$outdir/|g" > $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort

  # less -S $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort


  echo "Concatenating"

  bcftools concat \
  -f $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort \
  -O z \
  -n \
  -o $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

  # vcf-concat \
  # -f $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort \
  # -s | \
  # bgzip -c > $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

  tabix -p vcf $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz


}
export -f concatVCF

parallel -j8 concatVCF ::: 2L 2R 3L 3R 4 mitochondrion X Y
parallel -j8 concatVCF ::: 3L
