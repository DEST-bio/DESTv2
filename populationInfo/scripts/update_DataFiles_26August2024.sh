###
  ijob -A berglandlab_standard -c10 -p standard

  module load bcftools/1.17
  module load samtools/1.17

### define things

  prevVersion_VCF=dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz
  newVersion_VCF=dest.all.PoolSNP.001.50.24Aug2024.ann.vcf.gz

  wd=/project/berglandlab/DEST/vcf

### pull header
  bcftools view -h ${wd}/${prevVersion_VCF} > /scratch/aob2x/header_dest2.txt


  bcftools reheader \
  -h /scratch/aob2x/header_dest2.new.txt \
  -o ${wd}/${newVersion_VCF} \
  --threads 10 \
  ${wd}/${prevVersion_VCF}

  tabix ${wd}/${newVersion_VCF}

### convert new VCF to GDS calls the script below
  `DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.sh`
