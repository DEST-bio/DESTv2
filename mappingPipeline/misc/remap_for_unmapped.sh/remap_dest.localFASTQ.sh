#!/usr/bin/env bash
#
#SBATCH -J remap_fastq # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 15:00:00 ### 6 hours
#SBATCH --mem 50G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/remap.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/remap.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4559-aob2x

wd=/scratch/aob2x/dest
### grep -E "ES_ba_12|AT_gr_12" /scratch/aob2x/dest/DEST/populationInfo/samps.csv | cut -f1,13 -d',' > /scratch/aob2x/fastq/todl.csv
### run as: sbatch --array=2-20%5 /scratch/aob2x/DESTv2/mappingPipeline/misc/remap_for_unmapped.sh/remap_dest.localFASTQ.sh

#sbatch --array=23,25,26,35,42,43,45,46,48,49,50,51,52,53,54,55,56,57,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,119,120,130,131,135,136,137,138,139,140,144,145,146,147,148,149,150,151,152,153,156,157,158,159,161,162,163,164,165,166,167,168,169,170,171%20 /scratch/aob2x/DESTv2/mappingPipeline/misc/remap_for_unmapped.sh/remap_dest.localFASTQ.sh
### sacct -j 55855458
### cat /scratch/aob2x/dest/slurmOutput/remap.55835711_42.err

### sacct -j 55835711 | grep "FAIL" | grep "standard" | cut -f1 -d' ' | cut -f2 -d'_' | tr '\n' ','

module load gcc/11.4.0 sratoolkit/3.0.3 samtools/1.17 openmpi/4.1.4 bwa/0.7.17 picard/2.23.4 cutadapt/3.4
threads=10

#SLURM_ARRAY_TASK_ID=10

### get sample
  # sranum=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dros3.all.mapping.guide.txt | cut -f2 )
  # sample=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dros3.all.mapping.guide.txt | cut -f1 )
  sranum=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dest/missingSamples.DrosEU_3.delim | cut -f2 )
  sample=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dest/missingSamples.DrosEU_3.delim | cut -f1 )

  echo $sample " / " $sranum

### copy over

### trim
  if [ ! -f /scratch/aob2x/dest/fastq/${sranum}.trimmed1.fq ]; then
    echo "cutadapt start"

    cutadapt \
    -q 18 \
    --minimum-length 25 \
    -o /scratch/aob2x/dest/fastq/${sranum}.trimmed1.fq.gz \
    -p /scratch/aob2x/dest/fastq/${sranum}.trimmed2.fq.gz \
    -b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
    -B CAAGCAGAAGACGGCATACGAGAT \
    -O 15 \
    -n 3 \
    --cores=$threads \
    /project/berglandlab/DEST/raw_reads/DrosEU_3_Jan2023/${sranum}_1.fastq.gz \
    /project/berglandlab/DEST/raw_reads/DrosEU_3_Jan2023/${sranum}_2.fastq.gz

    gunzip /scratch/aob2x/dest/fastq/${sranum}.trimmed1.fq.gz
    gunzip /scratch/aob2x/dest/fastq/${sranum}.trimmed2.fq.gz
  fi

  echo "cutadapt done"
  ls -lh /scratch/aob2x/dest/fastq/*

### merge reads
  if [ ! -f /scratch/aob2x/dest/fastq/${sranum}.1_un.fq ]; then
    echo "bbmerge start"
    /scratch/aob2x/dest/bbmap/bbmerge.sh \
    in1=/scratch/aob2x/dest/fastq/${sranum}.trimmed1.fq \
    in2=/scratch/aob2x/dest/fastq/${sranum}.trimmed2.fq \
    out=/scratch/aob2x/dest/fastq/${sranum}.merged.fq \
    outu1=/scratch/aob2x/dest/fastq/${sranum}.1_un.fq \
    outu2=/scratch/aob2x/dest/fastq/${sranum}.2_un.fq
  fi

  echo "bbmerge end"
  ls -lh /scratch/aob2x/dest/fastq/*

### remap
  if [ ! -f /scratch/aob2x/dest/bam/${sample}.merged.bam ]; then
    echo "remap start"
    bwa mem -t $threads -M -R "@RG\tID:${sranum}\tSM:${sample}\tPL:illumina\tLB:lib1" \
    /scratch/aob2x/dest/remap_for_unmapped/ref/holo_dmel_6.12.fa \
    /scratch/aob2x/dest/fastq/${sranum}.1_un.fq \
    /scratch/aob2x/dest/fastq/${sranum}.2_un.fq | \
    samtools view -@ $threads -Sbh - > /scratch/aob2x/dest/bam/${sample}.merged_un.bam

    bwa mem -t $threads -M -R "@RG\tID:${sranum}\tSM:${sample}\tPL:illumina\tLB:lib1" \
    /scratch/aob2x/dest/remap_for_unmapped/ref/holo_dmel_6.12.fa \
    /scratch/aob2x/dest/fastq/${sranum}.merged.fq | \
    samtools view -@ $threads -Sbh - > /scratch/aob2x/dest/bam/${sample}.merged.bam
  fi
  echo "mapping done"


  java -jar $EBROOTPICARD/picard.jar MergeSamFiles \
  I=/scratch/aob2x/dest/bam/${sample}.merged.bam \
  I=/scratch/aob2x/dest/bam/${sample}.merged_un.bam \
  SO=coordinate \
  USE_THREADING=true \
  O=/scratch/aob2x/dest/bam/${sample}.sorted_merged.bam

  echo "remap done"
  ls -lh /scratch/aob2x/dest/bam/*

### get unmapped reads
  echo "get unmapped"
  samtools view -h -@ 20 -b -f 12 /scratch/aob2x/dest/bam/${sample}.sorted_merged.bam  > \
  /scratch/aob2x/dest/bam/${sample}.sorted_merged.unmapped.bam

### get reads mapping to non-Drosophila genomes
  # samtools idxstats ${inputFile} | grep -vE "2L|2R|3L|3R|4|X|Y|mitochondrion_genome|sim_2L|sim_2R|sim_3L|sim_3R|sim_4|sim_X|sim_mtDNA" | cut -f1,2 | awk '{print $1"\t"1"\t"$2}' > /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed
  # sed -i '$d' /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed

  samtools view -@ 20 -L /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed /scratch/aob2x/dest/bam/${sample}.sorted_merged.bam -b > \
  /scratch/aob2x/dest/bam/${sample}.sorted_merged.nonDros.bam

### index
  samtools index /scratch/aob2x/dest/bam/${sample}.sorted_merged.unmapped.bam
  samtools index /scratch/aob2x/dest/bam/${sample}.sorted_merged.nonDros.bam

  echo "remap done"
  ls -lh /scratch/aob2x/dest/bam/*
  ls -lh /scratch/aob2x/dest/bam/${sample}.sorted_merged.nonDros.bam

### clean up
  rm /scratch/aob2x/fastq/${sranum}.sra
  rm /scratch/aob2x/dest/fastq/${sranum}.trimmed1.fq
  rm /scratch/aob2x/dest/fastq/${sranum}.trimmed2.fq
  rm /scratch/aob2x/dest/fastq/${sranum}_1.fastq
  rm /scratch/aob2x/dest/fastq/${sranum}_2.fastq
  rm /scratch/aob2x/dest/fastq/${sranum}.merged.fq
  rm /scratch/aob2x/dest/fastq/${sranum}.1_un.fq
  rm /scratch/aob2x/dest/fastq/${sranum}.2_un.fq
  rm /scratch/aob2x/dest/bam/${sample}.merged.bam
  rm /scratch/aob2x/dest/bam/${sample}.merged_un.bam
  rm /scratch/aob2x/dest/bam/${sample}.sorted_merged.bam
