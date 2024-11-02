#!/usr/bin/env bash
#
#SBATCH -J updateDESTv2_24Aug2024 # A single job name for the array
#SBATCH -c 10 ### 10 cores
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00
#SBATCH --mem 5G
#SBATCH -o /scratch/aob2x/logs/updateDESTv2_24Aug2024.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/updateDESTv2_24Aug2024.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run as: sbatch /scratch/aob2x/DESTv2/populationInfo/scripts/update_DataFiles_26August2024.sh
### sacct -j 65726654
### cat /scratch/aob2x/logs/updateDESTv2_24Aug2024.*.err

  module load bcftools/1.17
  module load samtools/1.17
  module load gcc/11.4.0  openmpi/4.1.4 R/4.3.1; echo "R_LIBS_USER=~/R/goolf/4.3" > ~/.Renviron

### define things

  prevVersion_VCF=dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz
  newVersion_VCF=dest.all.PoolSNP.001.50.24Aug2024.ann.vcf.gz

  wd=/project/berglandlab/DEST/vcf

### pull header
  # bcftools view -h ${wd}/${prevVersion_VCF} > /scratch/aob2x/header_dest2.txt

### modify header
  #sed '1d' /scratch/aob2x/DESTv2/populationInfo/scripts/update_DataFiles_26August2024_conversionTable.csv > \
  #/scratch/aob2x/update_DataFiles_26August2024_conversionTable.noHeader.csv

  #cp /scratch/aob2x/header_dest2.txt /scratch/aob2x/header_dest2.new.txt

  #while read p; do
  #  newName=$( echo ${p} | cut -d',' -f1 )
  #  oldName=$( echo ${p} | cut -d',' -f2 )

  #  echo ${newName}
  #  sed -i "s/${oldName}/${newName}/g" /scratch/aob2x/header_dest2.new.txt

  #done < /scratch/aob2x/update_DataFiles_26August2024_conversionTable.noHeader.csv

  #cmp --silent /scratch/aob2x/header_dest2.txt /scratch/aob2x/header_dest2.new.txt || echo "files are different"

### rename
  #bcftools reheader \
  #-h /scratch/aob2x/header_dest2.new.txt \
  #-o ${wd}/${newVersion_VCF} \
  #--threads 10 \
  #${wd}/${prevVersion_VCF}

  #tabix ${wd}/${newVersion_VCF}

### convert new VCF to GDS calls the script below
  gdsfn=$( echo ${newVersion_VCF} | sed 's/vcf.gz/gds/g' )
  echo ${gdsfn}

  Rscript --vanilla \
  /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.R \
  /project/berglandlab/DEST/vcf/${newVersion_VCF} \
  /project/berglandlab/DEST/gds/${gdsfn} \
