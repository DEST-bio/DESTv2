#!/usr/bin/env bash
#
#SBATCH -J liftover_r5_to_r6 # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:15:00 ### 0.5 hours
#SBATCH --mem 6G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/liftover_r5_to_r6.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/liftover_r5_to_r6.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account biol8083

module load htslib

### note: this script uses the ram disk to write temporary files.

wd="/scratch/aob2x/dest"

#SLURM_ARRAY_TASK_ID=10; SLURM_JOB_ID=4
tmpdir=/dev/shm/$USER/
[ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}


pop=$( grep  "^${SLURM_ARRAY_TASK_ID}[[:space:]]" ${wd}/dgn/pops.delim | cut -f3 )
chr=$( grep  "^${SLURM_ARRAY_TASK_ID}[[:space:]]" ${wd}/dgn/pops.delim | cut -f2 )

zcat ${wd}/dest/wholeGenomeSyncData/${pop}_Chr${chr}.sync.gz |
awk '{
  print "chr"$1"\t"$2"\t"$2+1"\t"$3","$4",dm5_"$1"_"$2
}' > /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/${pop}_Chr${chr}.dm3.bed

### do liftover
~/liftOver \
/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/${pop}_Chr${chr}.dm3.bed \
${wd}/dgn/liftoverChains/dm3ToDm6.over.chain \
/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/${pop}_Chr${chr}.dm6.bed \
/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/${pop}_Chr${chr}.dm6.unmapped.bed

### get chromosome lengths from 6.12
### head -n6 /scratch/aob2x/dest/referenceGenome/r6/holo_dmel_6.12.fa.fai | cut -f1,2


if [ ${chr} = "2L" ]; then
  maxLen=23513712
elif [ ${chr} = "2R" ]; then
  maxLen=25286936
elif [ ${chr} = "3L" ]; then
  maxLen=28110227
elif [ ${chr} = "3R" ]; then
  maxLen=32079331
elif [ ${chr} = "X" ]; then
  maxLen=23542271
fi



cat /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/${pop}_Chr${chr}.dm6.bed | \
grep -E "^chr${chr}[[:space:]]" | \
sort -n -k2 - > ~/sort.bed | \
awk -v chr=${chr} -v chrLen=${maxLen} '{

  if($2==p) {
    next
  }
  if($2!=p+1) {
    for(i=p+1; i<=$2-1; i++) {
      print chr"\t"i"\tN\t0:0:0:0:0:0"
    }
    p=$2-1
  }
  if($2==p+1) {
    split($4, sp, ",")
    print chr"\t"$2"\t"sp[1]"\t"sp[2]

  }
}
{p=$2}
END{
  if(p<chrLen) {
    for(i=p+1; i<=chrLen; i++) {
      print chr"\t"i"\tN\t0:0:0:0:0:0"
    }
  }
}' | bgzip -c > ${wd}/dest/wholeGenomeSyncData/${pop}_Chr${chr}.gSYNC.gz

tabix -f -b 2 -s 1 -e 2 ${wd}/dest/wholeGenomeSyncData/${pop}_Chr${chr}.gSYNC.gz


### checks

#zcat ${wd}/dest/wholeGenomeSyncData/${pop}_Chr${chr}.gSYNC.gz | wc -l
#echo $maxLen
#
#zcat ${wd}/dest/wholeGenomeSyncData/${pop}_Chr${chr}.gSYNC.gz | head
#zcat ${wd}/dest/wholeGenomeSyncData/${pop}_Chr${chr}.gSYNC.gz | tail
#
#zcat ${wd}/dest/wholeGenomeSyncData/${pop}_Chr${chr}.gSYNC.gz | awk '$2!=p+1{print p"-"$2}{p=$2}'


#zcat ${wd}/dest/wholeGenomeSyncData/${pop}_Chr${chr}.gSYNC.gz | grep -m2 240254
