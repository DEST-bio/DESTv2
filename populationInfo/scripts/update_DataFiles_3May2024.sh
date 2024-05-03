###
ijob -A berglandlab_standard -c4 -p dev

module load bcftools/1.17
module load samtools/1.17

### define things
  prevVersion_VCF=dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz
  newVersion_VCF=dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz

  wd=/project/berglandlab/DEST/vcf

### pull header
  bcftools view -h ${wd}/${prevVersion_VCF} > /scratch/aob2x/header_dest2.txt

  sed 's/IT_Sas_Rec_1_2018-10-18/PT_Por_Rec_1_2018-10-18/g' /scratch/aob2x/header_dest2.txt > /scratch/aob2x/header_dest2.new.txt
  grep "IT_Sas_Rec_1_2018-10-18" /scratch/aob2x/header_dest2.txt
  grep "IT_Sas_Rec_1_2018-10-18" /scratch/aob2x/header_dest2.new.txt
  grep "PT_Por_Rec_1_2018-10-18" /scratch/aob2x/header_dest2.new.txt

  bcftools reheader \
  -h /scratch/aob2x/header_dest2.new.txt \
  -o ${wd}/${newVersion_VCF} \
  --threads 4 \
  ${wd}/${prevVersion_VCF}

  tabix ${wd}/${newVersion_VCF}

### convert new VCF to GDS
  `DESTv2/snpCalling/scatter_gather_annotate/vcf2gds.sh`
