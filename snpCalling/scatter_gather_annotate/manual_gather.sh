#!/bin/bash
#
#SBATCH -J manual_gather # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1:00:00 ### 1 hours
#SBATCH --mem 20G
#SBATCH -o /scratch/aob2x/DESTv2_output_26April2023/logs/manual_gather.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DESTv2_output_26April2023/logs/manual_gather.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

### cat /scratch/aob2x/DESTv2_output_SNAPE/logs/runSnakemake.49369837*.err

### sbatch /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/manual_gather.sh
### sacct -j 49411622
### cat /scratch/aob2x/DESTv2_output_26April2023/logs/manual_gather.49411616*.out

module purge

module load htslib/1.10.2 bcftools/1.9 intel/18.0 intelmpi/18.0 parallel/20200322

concatVCF() {

  popSet=all
  method=PoolSNP
  maf=001
  mac=50
  version=26April2023
  wd=/scratch/aob2x/DESTv2_output_26April2023


  chr=${1}

  echo "Chromosome: $chr"

  bcf_outdir="${wd}/sub_bcf"
  if [ ! -d $bcf_outdir ]; then
      mkdir $bcf_outdir
  fi
#
  outdir=$wd/sub_vcfs
  cd ${outdir}
#
  echo "generate list"
  ls -d *.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz | grep ${chr} | sort -t"_" -k2n,2 -k4g,4 \
  > $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort
#
  echo "Concatenating"
#
  bcftools concat \
  -f $outdir/vcfs_order.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.sort \
  -O z \
  -n \
  -o $bcf_outdir/dest.${chr}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz
}
export -f concatVCF

parallel concatVCF ::: 2L 2R 3L 3R 4 mitochondrion X Y
