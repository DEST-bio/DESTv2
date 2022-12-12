#!/usr/bin/env bash
#
#SBATCH -J makePopSync # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0:60:00 ### 5 minutes
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/makePopSync.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/makePopSync.%A_%a.err # Standard error
#SBATCH -p dev
#SBATCH --account berglandlab

module load htslib

#####################
### get jobs task ###
#####################

#SLURM_ARRAY_TASK_ID=60
##chr_i=$( echo "${SLURM_ARRAY_TASK_ID}%5+1" | bc )
##pop_i=$( echo "${SLURM_ARRAY_TASK_ID}%35+1" | bc )

pop=$( grep  "^${SLURM_ARRAY_TASK_ID}[[:space:]]" /scratch/aob2x/dest/dgn/pops.delim | cut -f3 )
chr=$( grep  "^${SLURM_ARRAY_TASK_ID}[[:space:]]" /scratch/aob2x/dest/dgn/pops.delim | cut -f2 )


##if [ ${chr_i} == "1" ]; then
##  chr="2L"
##elif [ ${chr_i} == "2" ]; then
##  chr="2R"
##elif [ ${chr_i} == "3" ]; then
##  chr="3L"
##elif [ ${chr_i} == "4" ]; then
##  chr="3R"
##elif [ ${chr_i} == "5" ]; then
##  chr="X"
##fi
##

echo $pop
echo $chr

############################
### paste per chromosome ###
############################

#if [ ! -f /scratch/aob2x/dest/dgn/csvData/${pop}_Chr${chr}.csv ]; then
#
#    paste -d',' /scratch/aob2x/dest/dgn/longData/${pop}_*_Chr${chr}*long  \
#    /scratch/aob2x/dest/dgn/csvData/${pop}_Chr${chr}.csv
#
#fi

samps=$( grep "^"${pop}"," < ${wd}/DEST/populationInfo/dpgp.ind.use.csv | cut -f2 -d',' | sed "s/$/_Chr${chr}/g" | tr '\n' '|' | sed 's/.$//' )

nFiles_obs=$( ls /scratch/aob2x/dest/dgn/longData/${pop}_*_Chr${chr}*long | grep -E "${samps}" | wc -l )

nFiles_exp=$( grep "^"${pop}"," < ${wd}/DEST/populationInfo/dpgp.ind.use.csv | cut -f9 -d',' | head -n1 )

echo ${pop}","${chr}","${nFiles_obs}","${nFiles_exp} >> /scratch/aob2x/dest/dgn/confirm_files

paste -d',' $( ls /scratch/aob2x/dest/dgn/longData/${pop}_*_Chr${chr}*long | grep -E "${samps}" ) > \
/scratch/aob2x/dest/dgn/csvData/${pop}_Chr${chr}.csv


##########################
### csv to sync format ###
##########################

paste -d' ' \
/scratch/aob2x/dest/referenceGenome/r5/${chr}.long \
/scratch/aob2x/dest/dgn/csvData/${pop}_Chr${chr}.csv | \
awk -F' ' -v chr=${chr} '
{
nN=gsub(/N/,"",$2)
nA=gsub(/A/,"",$2)
nT=gsub(/T/,"",$2)
nC=gsub(/C/,"",$2)
nG=gsub(/G/,"",$2)

nObs=nA+nT+nC+nG

print chr"\t"NR"\t"toupper($1)"\t"nA":"nT":"nC":"nG":"nN":0"

}' | bgzip -c > /scratch/aob2x/dest/dest/wholeGenomeSyncData/${pop}_Chr${chr}.sync.gz

#tabix -f -b 2 -s 1 -e 2 /scratch/aob2x/dest/dest/wholeGenomeSyncData/${pop}_Chr${chr}.sync.gz
