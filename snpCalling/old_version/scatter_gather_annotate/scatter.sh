#!/usr/bin/env bash

module purge
module load gcc/7.1.0 openmpi/3.1.4
#module load htslib bcftools parallel intel/18.0 intelmpi/18.0 mvapich2/2.3.1 R/3.6.3 python/3.6.6 vcftools/0.1.16
module load htslib/1.10.2 bcftools/1.9 parallel/20200322 intel/18.0 intelmpi/18.0 R/3.6.3 python/3.6.6 vcftools/0.1.16

### r, mvapch, parallel


## Run params
  popSet=${1}
  method=${2}
  maf=${3}
  mac=${4}
  version=${5}
  jobid=${6}
  job=$( echo $jobid | sed 's/_/,/g')
  script_dir=${7}
  wd=${8}
  pipeline_output=${9}

  #### popSet="all"; method="poolSNP"; maf="001"; mac=5; jobs="jobs.csv"; script_dir="/scratch/aob2x/DESTv2/snpCalling"; wd="/scratch/aob2x/DESTv2_output"; SLURM_JOB_ID=1;
  #### pipeline_output="/project/berglandlab/DEST/dest_mapped/"; job="2L,1,10000"

## working & temp directory
  outdir="${wd}/sub_vcfs" #### outdir=${wd}"/sub_vcfs"
    if [ ! -d $outdir ]; then
        mkdir $outdir
    fi

## get list of SNYC files based on popSet & method
### full list
  echo "foo: "${9}
  echo "pipeline_output: "${pipeline_output}
  echo $( ls ${pipeline_output}/*/*/*.sync.gz )

## get job
  cat ${wd}/jobs.csv
#job=$( cat ${wd}/jobs.txt | sed "${SLURM_ARRAY_TASK_ID}q;d" )
#jobid=$( echo ${job} | sed 's/,/_/g' )
  echo "jobid is " $jobid
  echo "job is " $job

## set up RAM disk
  [ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
  [ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
  tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}

  echo "Temp dir is $tmpdir"

## get sub section
  subsection () {
    syncFile=${1}
    job=${2}
    jobid=$( echo ${job} | sed 's/,/_/g' )
    tmpdir=${3}

    pop=$( echo ${syncFile} | rev | cut -f1 -d'/' | rev | sed 's/.masked.sync.gz//g' )

    chr=$( echo $job | cut -f1 -d',' )
    start=$( echo $job | cut -f2 -d',' )
    stop=$( echo $job | cut -f3 -d',' )

    echo ${pop}_${jobid}

    tabix -b 2 -s 1 -e 2 \
    ${syncFile} \
    ${chr}:${start}-${stop} > ${tmpdir}/${pop}_${jobid}

  }
  export -f subsection

  echo "subset"

  if [[ "${method}" == "SNAPE" && "${popSet}" == "PoolSeq" ]]; then
    echo "SNAPE" ${method}
    parallel -j 4 subsection ::: $( ls ${pipeline_output}/*/*/*.masked.sync.gz | tr '  ' '\n' | grep "SNAPE" | grep "monomorphic" | tr '\n' ' ' ) ::: ${job} ::: ${tmpdir}
  elif [[ "${method}" == "PoolSNP" && "${popSet}" == "all" ]]; then
    echo "PoolSNP" ${method}
    parallel -j 4 subsection ::: $( ls ${pipeline_output}/*/*/*.masked.sync.gz | tr '  ' '\n' | grep -v "SNAPE" ) ::: ${job} ::: ${tmpdir}
  elif [[ "${method}" == "PoolSNP" && "${popSet}" == "PoolSeq" ]]; then
    echo "PoolSNP" ${method}
    parallel -j 4 subsection ::: $( ls ${pipeline_output}/*/*/*.masked.sync.gz | tr '  ' '\n' | grep -v "SNAPE" | grep -v "DGN" ) ::: ${job} ::: ${tmpdir}
  fi




### paste function
  echo "paste"
  Rscript --no-save --no-restore ${script_dir}/scatter_gather_annotate/paste.R ${job} ${tmpdir} ${method}

### run through SNP calling
  echo "SNP calling"

  if [[ "${method}" == "SNAPE" ]]; then
    echo $method
    cat ${tmpdir}/allpops.${method}.sites | python ${script_dir}/PoolSNP/PoolSnp.py \
    --sync - \
    --min-cov 4 \
    --max-cov 0.95 \
    --miss-frac 0.5 \
    --min-count 0 \
    --min-freq 0 \
    --posterior-prob 0.9 \
    --SNAPE \
    --names $( cat ${tmpdir}/allpops.${method}.names |  tr '\n' ',' | sed 's/,$//g' )  > ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf

  elif [[ "${method}"=="PoolSNP" ]]; then
    echo $method

    cat ${tmpdir}/allpops.${method}.sites | python ${script_dir}/PoolSNP/PoolSnp.py \
    --sync - \
    --min-cov 4 \
    --max-cov 0.95 \
    --min-count ${mac} \
    --min-freq 0.${maf} \
    --miss-frac 0.5 \
    --names $( cat ${tmpdir}/allpops.${method}.names |  tr '\n' ',' | sed 's/,$//g' )  > ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf
  fi

### compress and clean up
  echo "compress and clean"
  # cp ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf
  # cat ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf | bgzip -c > ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz
   cat ${tmpdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf | vcf-sort | bgzip -c > ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz

  tabix -p vcf ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz


  module load gcc/9.2.0 bedtools/2.29.2

  bedtools intersect -sorted -v -header \
  -b ${script_dir}/scatter_gather_annotate/repeat_bed/repeats.sort.bed.gz \
  -a ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.vcf.gz |
  bgzip -c > \
  ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

  tabix -p vcf ${outdir}/${jobid}.${popSet}.${method}.${maf}.${mac}.${version}.norep.vcf.gz

  #echo "vcf -> bcf "
  #bcftools view -Ou ${tmpdir}/${jobid}.vcf.gz > ${outdir}/${jobid}.bcf

  rm -fr ${tmpdir}

### done
  echo "done"
