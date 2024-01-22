#!/usr/bin/env bash
#
#SBATCH -J remap_fastq_individual # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 15:00:00 ### 6 hours
#SBATCH --mem 75G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/remap.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/remap.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4559-aob2x

wd=/scratch/aob2x/dest
### run as: sbatch --array=2 /scratch/aob2x/DESTv2/mappingPipeline/misc/remap_for_unmapped.sh/remap_individual.sh
### sacct -j 57080425
### cat /scratch/aob2x/dest/slurmOutput/remap.57072677_2.out | tail

###   samtools idxstats /project/berglandlab/DEST/dest_mapped/Cville/US_Vir_Cha_1_2016-07-08/US_Vir_Cha_1_2016-07-08.original.bam | grep -vE "2L|2R|3L|3R|4|X|Y|mitochondrion_genome|sim_2L|sim_2R|sim_3L|sim_3R|sim_4|sim_X|sim_mtDNA" | cut -f1,2 | awk '{print $1"\t"1"\t"$2}' > /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed
###   sed -i '$d' /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed


module load gcc/11.4.0 sratoolkit/3.0.3 samtools/1.17 openmpi/4.1.4 bwa/0.7.17 picard/2.23.4 cutadapt/3.4
threads=10

#SLURM_ARRAY_TASK_ID=528
#SLURM_ARRAY_TASK_ID=2
## SLURM_ARRAY_TASK_ID=200

### combine keep lists:
### four columns: "sranum", "sampleId", "type" (iso or wild), "set" (2016 or OW)

  # cat /standard/vol186/bergland-lab/alyssa/2016filters/keepList_v4.txt | awk '{print $1","$1",CartersMountainIsofemale,2016"}' > /scratch/aob2x/dest/individual_data/samples2use.csv
  # grep "CM" /project/berglandlab/Dmel_Single_Individuals/Overwintering_2018_2019/Individuals_Metadata_OW_n5_n7_corrected.csv | \
  # cut -f20,28,22 -d',' | sed 's/ //g' | sed 's/\"//g' | awk -F',' '{print $1","$3","$2",OW"}' >> /scratch/aob2x/dest/individual_data/samples2use.csv

### 2016 flies: 4 sets of PE files per sample
  sranum=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dest/individual_data/samples2use.csv | cut -f1 -d',')
  sample=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dest/individual_data/samples2use.csv | cut -f2 -d',')
  set=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dest/individual_data/samples2use.csv | cut -f4 -d',')

  echo $sample " / " $sranum " / " $set

### standardize input files
  if [[ "$set" == "2016" ]]; then
    echo "2016 sample"
    # ls -lh /standard/vol186/bergland-lab/alyssa/2016fastq/*_1_${sranum}*.gz
    # ls -lh /standard/vol186/bergland-lab/alyssa/2016fastq/*_2_${sranum}*.gz

    cat /standard/vol186/bergland-lab/alyssa/2016fastq/*_1_${sample}*.gz > /scratch/aob2x/dest/individual_data/fastq/${sranum}_1.fastq.gz
    cat /standard/vol186/bergland-lab/alyssa/2016fastq/*_2_${sample}*.gz > /scratch/aob2x/dest/individual_data/fastq/${sranum}_2.fastq.gz

    r1_filename=/scratch/aob2x/dest/individual_data/fastq/${sranum}_1.fastq.gz
    r2_filename=/scratch/aob2x/dest/individual_data/fastq/${sranum}_2.fastq.gz

  else
    echo "OW sample"
    r1_filename=/project/berglandlab/Dmel_Single_Individuals/Overwintering_2018_2019/raw_data/usftp21.novogene.com/raw_data/${sranum}/*1.fq.gz
    r2_filename=/project/berglandlab/Dmel_Single_Individuals/Overwintering_2018_2019/raw_data/usftp21.novogene.com/raw_data/${sranum}/*2.fq.gz
  fi

### trim adapters
  if [ ! -f /scratch/aob2x/dest/individual_data/fastq/${sranum}.trimmed1.fq.gz ]; then
    echo "cutadapt start"
    cutadapt \
    -q 18 \
    --minimum-length 25 \
    -o /scratch/aob2x/dest/individual_data/fastq/${sranum}.trimmed1.fq.gz \
    -p /scratch/aob2x/dest/individual_data/fastq/${sranum}.trimmed2.fq.gz \
    -b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
    -B CAAGCAGAAGACGGCATACGAGAT \
    -O 15 \
    -n 3 \
    --cores=$threads \
    ${r1_filename} \
    ${r2_filename}
  fi

  echo "cutadapt done"

### merge reads
  if [ ! -f /scratch/aob2x/dest/fastq/${sranum}.1_un.fq ]; then
    echo "bbmerge start"
    /scratch/aob2x/dest/bbmap/bbmerge.sh \
    in1=/scratch/aob2x/dest/individual_data/fastq/${sranum}.trimmed1.fq.gz \
    in2=/scratch/aob2x/dest/individual_data/fastq/${sranum}.trimmed2.fq.gz \
    out=/scratch/aob2x/dest/individual_data/fastq/${sranum}.merged.fq.gz \
    outu1=/scratch/aob2x/dest/individual_data/fastq/${sranum}.1_un.fq.gz \
    outu2=/scratch/aob2x/dest/individual_data/fastq/${sranum}.2_un.fq.gz
  fi

  echo "bbmerge end"
  ls -lh /scratch/aob2x/dest/fastq//scratch/aob2x/dest/bam/*${sample}*

### remap
  if [ ! -f /scratch/aob2x/dest/bam/${sample}.merged.bam ]; then
    echo "remap start"
    bwa mem -t $threads -M -R "@RG\tID:${sranum}\tSM:${sample}\tPL:illumina\tLB:lib1" \
    /scratch/aob2x/dest/remap_for_unmapped/ref/holo_dmel_6.12.fa \
    /scratch/aob2x/dest/individual_data/fastq/${sranum}.1_un.fq.gz \
    /scratch/aob2x/dest/individual_data/fastq/${sranum}.2_un.fq.gz | \
    samtools view -@ $threads -Sbh - > /scratch/aob2x/dest/individual_data/bam/${sample}.merged_un.bam

    bwa mem -t $threads -M -R "@RG\tID:${sranum}\tSM:${sample}\tPL:illumina\tLB:lib1" \
    /scratch/aob2x/dest/remap_for_unmapped/ref/holo_dmel_6.12.fa \
    /scratch/aob2x/dest/individual_data/fastq/${sranum}.merged.fq.gz | \
    samtools view -@ $threads -Sbh - > /scratch/aob2x/dest/individual_data/bam/${sample}.merged.bam
  fi
  echo "mapping done"

### merge
  java -jar $EBROOTPICARD/picard.jar MergeSamFiles \
  I=/scratch/aob2x/dest/individual_data/bam/${sample}.merged.bam \
  I=/scratch/aob2x/dest/individual_data/bam/${sample}.merged_un.bam \
  SO=coordinate \
  USE_THREADING=true \
  O=/scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.bam

  echo "merging done"
  ls -lh /scratch/aob2x/dest/individual_data/bam/*

### get reads mapping to mitochondria
  # inputFile=/project/berglandlab/DEST/dest_mapped/Cville/US_Vir_Cha_1_2018-11-29/US_Vir_Cha_1_2018-11-29.original.bam
  # samtools idxstats ${inputFile} | grep -E "mitochondrion_genom|sim_mtDNA" | cut -f1,2 | awk '{print $1"\t"1"\t"$2}' > /scratch/aob2x/DESTv2_unmapped_reads/mitoGenome.bed
  #sed -i '$d' /scratch/aob2x/DESTv2_unmapped_reads/mitoGenome.bed
  echo "getting mito reads"

  samtools view -@ 20 -L /scratch/aob2x/DESTv2_unmapped_reads/mitoGenome.bed /scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.bam -b > \
  /scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.mito.bam

### get unmapped reads
  echo "get unmapped"
  samtools view -h -@ 20 -b -f 12 /scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.bam  > \
  /scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.unmapped.bam

### get reads mapping to non-Drosophila genomes
  # samtools idxstats ${inputFile} | grep -vE "2L|2R|3L|3R|4|X|Y|mitochondrion_genome|sim_2L|sim_2R|sim_3L|sim_3R|sim_4|sim_X|sim_mtDNA" | cut -f1,2 | awk '{print $1"\t"1"\t"$2}' > /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed
  # sed -i '$d' /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed
  echo "get unmapped"

  samtools view -@ 20 -L /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed /scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.bam -b > \
  /scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.nonDros.bam

### index
  samtools index /scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.mito.bam
  samtools index /scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.nonDros.bam
  samtools index /scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.unmapped.bam
  samtools index /scratch/aob2x/dest/individual_data/bam/${sample}.sorted_merged.bam

  rm /scratch/aob2x/dest/individual_data/bam/${sample}.merged.bam
  rm /scratch/aob2x/dest/individual_data/bam/${sample}.merged_un.bam

  echo "remap done"
