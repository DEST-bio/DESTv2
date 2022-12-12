#!/bin/bash

### this variable points to the directory where all of the mapping folders are
mappingDir=/mnt/spicy_2/dest/mapped


getSimCont () {



	rootDir=${1}
	samp=${2}

	(>&2 echo ${samp})

	samtools idxstats \
	${rootDir}/${samp}/mapping/${samp}*.decon.bam-sim.bam 2>/tmp/null | \
	grep -E "^(2L|2R|3L|3R|X){1}[[:space:]]" |
	sed "s/^/${samp}\tsim\t/g"

	samtools idxstats \
	${rootDir}/${samp}/mapping/${samp}*.decon.bam-mel.bam 2>/tmp/null | \
	grep -E "^(2L|2R|3L|3R|X){1}[[:space:]]" |
	sed "s/^/${samp}\tmel\t/g"

}
export -f getSimCont


nohup parallel --gnu -j1 getSimCont ::: ${mappingDir} ::: $( ls ${mappingDir} ) > ~/sim_mel_mappingRates.DrosRTEC.delim &
