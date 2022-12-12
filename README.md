# DEST
  > Scripts for mapping and quality control of DEST dataset

## Metadata
  > `populationInfo\`: Has supplemental data from the DrosEU, DrosRTEC, and DPGP files to make a unified meta-datafile. This datafile is
  supplemental table 1 of Kapun et al 2021. The meta-data file includes collection information, library quality filtering, inversion frequencies,
  weather stations and WorldClim environmental data.

## PoolSeq mapping pipeline
  > `mappingPipeline\`: Contains dockerized mapping pipeline. Downloads data, produces bam files, filter files, gSYNC files

## Incorporate Drosophila Genome Nexus data
  > `add_DGN_datda\`: Downloads, formats, lifts-over the Drosophila Genome Nexus data into gSYNC format

## SNP calling
  > `snpCalling\`: SNP calling based on gSYNC files. Runs PoolSNP and SNAPE using snakemake pipeline

## Example Scripts
  > `examples`: Provides a basic script to pull allele frequencies from a GDS table and also calculates the effective coverage
