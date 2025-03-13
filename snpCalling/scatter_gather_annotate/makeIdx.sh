#!/bin/bash
#
#SBATCH -J manual_annotate # A single job name for the array
#SBATCH --ntasks-per-node=20 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 14:00:00 ### 1 hours
#SBATCH --mem 40G
#SBATCH -o /scratch/aob2x/DESTv2_output_26April2023/logs/manual_annotate.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/DESTv2_output_26April2023/logs/manual_annotate.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard

### cat /scratch/aob2x/DESTv2_output_SNAPE/logs/runSnakemake.49369837*.err

### sbatch /scratch/aob2x/DESTv2/snpCalling/scatter_gather_annotate/manual_annotate.sh
### sacct -j 49572492
### cat /scratch/aob2x/DESTv2_output_26April2023/logs/manual_annotate.49572492*.out

module purge

module load samtools/1.17 htslib/1.17

syncs=$( ls /project/berglandlab/DEST/dest_mapped/*/*/*.monomorphic.masked.sync.gz )

for i in /project/berglandlab/DEST/dest_mapped/*/*/*.monomorphic.masked.sync.gz; do
   #echo "$i"
   if [ ! -f $i.tbi ]; then
     echo $i
   fi
   # or do whatever with individual element of the array
done


gunzip /project/berglandlab/DEST/dest_mapped/DEST_SA/CO_Cun_LaV_1_2019-02-21/CO_Cun_LaV_1_2019-02-21.SNAPE.monomorphic.masked.sync.gz
bgzip /project/berglandlab/DEST/dest_mapped/DEST_SA/CO_Cun_LaV_1_2019-02-21/CO_Cun_LaV_1_2019-02-21.SNAPE.monomorphic.masked.sync

tabix -s 1 -b 2 -e 2 /project/berglandlab/DEST/dest_mapped/DEST_SA/CO_Cun_LaV_1_2019-02-21/CO_Cun_LaV_1_2019-02-21.SNAPE.monomorphic.masked.sync.gz
