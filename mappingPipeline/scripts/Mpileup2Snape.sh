#!/bin/bash

mpileup=$1
output=$2
sample=$3
theta=$4
D=$5
priortype=$6
fold=$7
nflies=$8

cd $output/$sample/

awk '{if (last != $1) close(last); print >> $1; last = $1}' $mpileup

for chr in {2L,2R,3L,3R,4}; do
  snape-pooled -nchr $(($nflies*2)) -theta $theta -D $D -priortype $priortype -fold $fold < ${chr} > ${chr}-$sample-SNAPE.txt
done

for chr in {X,Y}; do
  snape-pooled -nchr $nflies -theta $theta -D $D -priortype $priortype -fold $fold < ${chr} > ${chr}-$sample-SNAPE.txt
done

for chr in {2L,2R,3L,3R,4,X,Y,mitochondrion_genome}; do
  rm $chr
done

cat *-$sample-SNAPE.txt  > ${sample}.SNAPE.output.txt

rm *-$sample-SNAPE.txt
