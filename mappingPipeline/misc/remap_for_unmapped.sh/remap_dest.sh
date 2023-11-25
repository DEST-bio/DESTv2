#!/usr/bin/env bash
#
#SBATCH -J remap_fastq # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 10:00:00 ### 6 hours
#SBATCH --mem 50G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/remap.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/remap.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4559-aob2x

wd=/scratch/aob2x/dest
### grep -E "ES_ba_12|AT_gr_12" /scratch/aob2x/dest/DEST/populationInfo/samps.csv | cut -f1,13 -d',' > /scratch/aob2x/fastq/todl.csv
### run as: sbatch --array=2-754%10 /scratch/aob2x/DESTv2/mappingPipeline/misc/remap_for_unmapped.sh/remap_dest.sh
### sacct -j 55233826
### cat /scratch/aob2x/dest/slurmOutput/remap.55114839_2.err

module load sratoolkit/2.10.5 samtools/1.9 gcc/9.2.0 bwa/0.7.17 picard/2.23.4 cutadapt/3.4
threads=10

#SLURM_ARRAY_TASK_ID=2

### get sample
  sranum=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dest/dest_v2.samps_8Jun2023.csv | cut -f31 -d',' )
  sample=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dest/dest_v2.samps_8Jun2023.csv | cut -f1 -d',' )

  echo $sample " / " $sranum

### download if necesary
  if [ ! -f "/scratch/aob2x/dest/fastq/${sranum}.sra" ]; then
    prefetch \
    -o /scratch/aob2x/dest/fastq/${sranum}.sra \
    ${sranum}
  fi

  if [ ! -f "/scratch/aob2x/dest/fastq/${sranum}_1.fastq" ]; then

    fasterq-dump \
    --split-files \
    --split-3 \
    --outfile /scratch/aob2x/dest/fastq/${sranum} \
    -e 10 \
    /scratch/aob2x/dest/fastq/${sranum}.sra

    # gzip /scratch/aob2x/dest/fastq/${sranum}_1.fastq
    # gzip /scratch/aob2x/dest/fastq/${sranum}_2.fastq

    # rm /scratch/aob2x/fastq/${sranum}.sra
  fi

### trim
  if [ ! -f /scratch/aob2x/dest/fastq/${sranum}.trimmed1.fq ]; then
    echo "cutadapt start"

    cutadapt \
    -q 18 \
    --minimum-length 25 \
    -o /scratch/aob2x/dest/fastq/${sranum}.trimmed1.fq \
    -p /scratch/aob2x/dest/fastq/${sranum}.trimmed2.fq \
    -b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
    -B CAAGCAGAAGACGGCATACGAGAT \
    -O 15 \
    -n 3 \
    --cores=$threads \
    /scratch/aob2x/dest/fastq/${sranum}_1.fastq \
    /scratch/aob2x/dest/fastq/${sranum}_2.fastq
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
