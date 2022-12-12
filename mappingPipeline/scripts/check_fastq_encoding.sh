#!/bin/bash
#
#SBATCH -J check_fastq_encoding # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 00:30:00 ### 1 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/check_fastq_encoding.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/check_fastq_encoding.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


wd=/scratch/aob2x/dest

chmod +x ${wd}/DEST_freeze1/mappingPipeline/scripts/guess_encoding.py

if test -f ${wd}/fastq/qualEncodings.delim; then
    rm ${wd}/fastq/qualEncodings.delim
    touch ${wd}/fastq/qualEncodings.delim
else
    touch ${wd}/fastq/qualEncodings.delim
fi

for f in ${wd}/fastq/*fastq.gz; do
  ### f=/scratch/aob2x/dest/fastq/SRX2885350_2.fastq.gz
  zcat $f | head -n 5000 | awk 'NR % 4 == 0' | \
  ${wd}/DEST/mappingPipeline/scripts/guess_encoding.py | grep -v "#" |
  sed -E "s|^|${f}\t|g" >> ${wd}/fastq/qualEncodings.delim
done
