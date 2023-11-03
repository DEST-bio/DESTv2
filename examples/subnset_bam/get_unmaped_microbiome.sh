#!/usr/bin/env bash
#
#SBATCH -J fst # A single job name for the array
#SBATCH -c 20 ### 10 cores
#SBATCH -N 1 # on one node
#SBATCH -t 0:90:00
#SBATCH --mem 20G
#SBATCH -o /scratch/aob2x/logs/dest_fst.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/logs/dest_fst.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol4559-aob2x

### Nov 2 2023
### run as: sbatch --array=1-743 /scratch/aob2x/DESTv2/examples/subnset_bam/get_unmaped_microbiome.sh

### sacct -j 54804087
### cat /scratch/aob2x/logs/dest_fst.54723955_46.out
### SLURM_ARRAY_TASK_ID=2

### modules
  module load samtools
  #module load intel/18.0  intelmpi/18.0 R/4.1.1

### get file list
  inputFile=$( ls -d /project/berglandlab/DEST/dest_mapped/*/*/*.original.bam | sed "${SLURM_ARRAY_TASK_ID}q;d" )
  ls -lh $inputFile
  fileStem=$( echo $inputFile | rev | cut -f1 -d'/' | rev )

### get unmapped reads
  samtools view -h -@ 20 -b -f 4 $inputFile  > /scratch/aob2x/DESTv2_unmapped_reads_v2/unmapped.${fileStem}

### get reads mapping to non-Drosophila genomes
  # samtools idxstats ${inputFile} | grep -vE "2L|2R|3L|3R|4|X|Y|mitochondrion_genome|sim_2L|sim_2R|sim_3L|sim_3R|sim_4|sim_X|sim_mtDNA" | cut -f1,2 | awk '{print $1"\t"1"\t"$2}' > /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed
  # sed -i '$d' /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed

  samtools view -@ 20 -L /scratch/aob2x/DESTv2_unmapped_reads/nonDrosGenome.bed $inputFile -b > /scratch/aob2x/DESTv2_unmapped_reads_v2/nonDros.${fileStem}

### index
  samtools index /scratch/aob2x/DESTv2_unmapped_reads_v2/unmapped.${fileStem}
  samtools index /scratch/aob2x/DESTv2_unmapped_reads_v2/nonDros.${fileStem}

### log
  ls -lh /scratch/aob2x/DESTv2_unmapped_reads_v2/*
