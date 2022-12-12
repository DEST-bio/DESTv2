#!/usr/bin/env bash
#
#SBATCH -J unpack # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 1G
#SBATCH -o /scratch/aob2x/dest/slurmOutput/unpack.%A_%a.out # Standard output
#SBATCH -e /scratch/aob2x/dest/slurmOutput/unpack.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### unpack datasets
  #SLURM_ARRAY_TASK_ID=2

  name=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | cut -f2 -d' ' )
  url=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | cut -f3 -d' ' )
  fileName=$( grep ^${SLURM_ARRAY_TASK_ID} /scratch/aob2x/dest/DEST/add_DGN_data/dgn.list | rev | cut -f1 -d'/' | rev )


cd /scratch/aob2x/dest/dgn/wideData

if [ "${SLURM_ARRAY_TASK_ID}" == "1" ]; then
  ### DPGP2

  #tar -zvx \
  #--directory=/scratch/aob2x/dest/dgn/rawData \
  #--file=/scratch/aob2x/dest/dgn/rawData/${fileName}

  cd /scratch/aob2x/dest/dgn/rawData/
  tar -xvf dpgp2_Chr2L_sequences.tar
  tar -xvf dpgp2_Chr2R_sequences.tar
  tar -xvf dpgp2_Chr3L_sequences.tar
  tar -xvf dpgp2_Chr3R_sequences.tar
  tar -xvf dpgp2_ChrX_sequences.tar

  for i in *.seq; do
    pre=$( echo $i |  cut -f1 -d'_' | cut -c1,2 )
    mv $i /scratch/aob2x/dest/dgn/wideData/${pre}_${i}
  done

  #ls *.seq | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/DPGP2_{}

elif [ "${SLURM_ARRAY_TASK_ID}" == "2" ]; then
  ###DPGP3

  tar -zvx \
  --directory=/scratch/aob2x/dest/dgn/rawData \
  --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

  cd dpgp3_sequences

  tar -xvf dpgp3_Chr2L.tar
  tar -xvf dpgp3_Chr2R.tar
  tar -xvf dpgp3_Chr3L.tar
  tar -xvf dpgp3_Chr3R.tar
  tar -xvf dpgp3_ChrX.tar

  ls ZI* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/ZI_{}
  cd dpgp3_ChrX
  ls ZI* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/ZI_{}


elif [ "${SLURM_ARRAY_TASK_ID}" == "3" ]; then
  ### DGRP
  #tar -jvx \
  #--directory=/scratch/aob2x/dest/dgn/rawData \
  #--file=/scratch/aob2x/dest/dgn/rawData/${fileName}

  cd /scratch/aob2x/dest/dgn/rawData/
  tar -xvf dgrp_Chr2L.tar
  tar -xvf dgrp_Chr2R.tar
  tar -xvf dgrp_Chr3L.tar
  tar -xvf dgrp_Chr3R.tar
  tar -xvf dgrp_ChrX.tar

  ls RAL* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/RAL_{}

elif [ "${SLURM_ARRAY_TASK_ID}" == "5" ]; then
    tar -zvx \
    --directory=/scratch/aob2x/dest/dgn/rawData \
    --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences
    tar -xvf CLARK_Chr2L_sequences.tar
    tar -xvf CLARK_Chr2R_sequences.tar
    tar -xvf CLARK_Chr3L_sequences.tar
    tar -xvf CLARK_Chr3R_sequences.tar
    tar -xvf CLARK_ChrX_sequences.tar

    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences/CLARK_Chr2L
    ls B* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/B_{}
    ls I* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/I_{}
    ls N* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/N_{}
    ls T* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/T_{}

    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences/CLARK_Chr2R
    ls B* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/B_{}
    ls I* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/I_{}
    ls N* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/N_{}
    ls T* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/T_{}

    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences/CLARK_Chr3L
    ls B* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/B_{}
    ls I* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/I_{}
    ls N* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/N_{}
    ls T* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/T_{}

    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences/CLARK_Chr3R
    ls B* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/B_{}
    ls I* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/I_{}
    ls N* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/N_{}
    ls T* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/T_{}

    cd /scratch/aob2x/dest/dgn/rawData/CLARK_sequences/CLARK_ChrX
    ls B* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/B_{}
    ls I* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/I_{}
    ls N* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/N_{}
    ls T* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/T_{}

elif [ "${SLURM_ARRAY_TASK_ID}" == "6" ]; then
  tar -zvx \
  --directory=/scratch/aob2x/dest/dgn/rawData \
  --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences
  tar -xvf NUZHDIN_Chr2L_sequences.tar
  tar -xvf NUZHDIN_Chr2R_sequences.tar
  tar -xvf NUZHDIN_Chr3L_sequences.tar
  tar -xvf NUZHDIN_Chr3R_sequences.tar
  tar -xvf NUZHDIN_ChrX_sequences.tar

  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences/NUZHDIN_Chr2L
  ls w* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/W_{}

  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences/NUZHDIN_Chr2R
  ls w* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/W_{}

  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences/NUZHDIN_Chr3L
  ls w* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/W_{}

  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences/NUZHDIN_Chr3R
  ls w* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/W_{}

  cd /scratch/aob2x/dest/dgn/rawData/NUZHDIN_sequences/NUZHDIN_ChrX
  ls w* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/W_{}

elif [ "${SLURM_ARRAY_TASK_ID}" == "7" ]; then
  tar -zvx \
  --directory=/scratch/aob2x/dest/dgn/rawData \
  --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

  cd /scratch/aob2x/dest/dgn/rawData/POOL_sequences
  tar -xvf POOL_Chr2L_Sequences.tar
  tar -xvf POOL_Chr2R_Sequences.tar
  tar -xvf POOL_Chr3L_Sequences.tar
  tar -xvf POOL_Chr3R_Sequences.tar
  tar -xvf POOL_ChrX_Sequences.tar

  cd /scratch/aob2x/dest/dgn/rawData/POOL_sequences/POOL_Chr2L
  for i in *.seq; do
    pre=$( echo $i |  cut -f1 -d'_' | cut -c1,2 )
    mv $i /scratch/aob2x/dest/dgn/wideData/${pre}_${i}
  done

  cd /scratch/aob2x/dest/dgn/rawData/POOL_sequences/POOL_Chr2R
  for i in *.seq; do
    pre=$( echo $i |  cut -f1 -d'_' | cut -c1,2 )
    mv $i /scratch/aob2x/dest/dgn/wideData/${pre}_${i}
  done

  cd /scratch/aob2x/dest/dgn/rawData/POOL_sequences/POOL_Chr3L
  for i in *.seq; do
    pre=$( echo $i |  cut -f1 -d'_' | cut -c1,2 )
    mv $i /scratch/aob2x/dest/dgn/wideData/${pre}_${i}
  done

  cd /scratch/aob2x/dest/dgn/rawData/POOL_sequences/POOL_Chr3R
  for i in *.seq; do
    pre=$( echo $i |  cut -f1 -d'_' | cut -c1,2 )
    mv $i /scratch/aob2x/dest/dgn/wideData/${pre}_${i}
  done

  cd /scratch/aob2x/dest/dgn/rawData/POOL_sequences/POOL_ChrX
  for i in *.seq; do
    pre=$( echo $i |  cut -f1 -d'_' | cut -c1,2 )
    mv $i /scratch/aob2x/dest/dgn/wideData/${pre}_${i}
  done

elif [ "${SLURM_ARRAY_TASK_ID}" == "8" ]; then
  tar -zvx \
  --directory=/scratch/aob2x/dest/dgn/rawData \
  --file=/scratch/aob2x/dest/dgn/rawData/${fileName}

  cd /scratch/aob2x/dest/dgn/rawData/
  ls Simulans* | xargs -t -I{} mv {} /scratch/aob2x/dest/dgn/wideData/SIM_{}
fi
