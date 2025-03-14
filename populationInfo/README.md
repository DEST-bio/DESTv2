# Sample information

## Description for DESTv2 sample metadata

## sample metadata
> `DESTv2/populationInfo/CHANGE_LOG.md` has info on the latest versions and changes between dated versions

## scripts
> Several scripts that generate the sample metadata
>  `DESTv2/populationInfo/scripts/makeJointSampleInfo_DESTv2.R`
> There are also a bunch of update scripts that fix various naming erorrs or add more info










## Description for DESTv1 sample medatdata
>  This directory contains scripts to generate meta-data files for the DEST2 dataset. The 25Feb2023 version is the most recent and matches with the SNP calls that use the same date
>
## Sample metadata
  > `DEST_freeze1/populationInfo/makeJointSampleInfo.R` generates several files:
  > 1. `DEST_freeze1/populationInfo/samps.csv` <br> Contains collection information (locality, date, SRA accession, weather station ID, etc). This is Supplemental Table 1 <br>
  > 2. `DEST_freeze1/populationInfo/dest.worldclim.csv` contains WorldClim data for sampling localities

## Library sequencing statistics
  > `DEST_freeze1/populationInfo/sequencingStats/sequencingSummaryStats.R` generates data that is part of Figure 3 and is found in Supplemental Table 2 <br>
  > 1. `DEST_freeze1/populationInfo/sequencingStats/rd.csv` Average read depth
  > 2. `DEST_freeze1/populationInfo/sequencingStats/pcr.csv` PCR duplicate rate
  > 3. `DEST_freeze1/populationInfo/sequencingStats/simulans.csv` D. simulans contamination rate

## Sample filtering based on pn/ps & the number of private SNPs
  > `DEST_freeze1/populationInfo/sequencingStats/sequencingSummaryStats.R` generates the filter file
  > 1. `DEST_freeze1/populationInfo/sampleFilter/classify_pops.txt` is the filter file

## Inversion frequency estimation
  > `DEST_freeze1/populationInfo/Inversions/inversions.sh` estimates inversion frequencies from the data
  > 1. `DEST_freeze1/populationInfo/Inversions/inversion.af` are the estimates

## World clim
  > `DEST_freeze1/populationInfo/makeJointSampleInfo.R` generates the data internally, but it is save to a file for quicker updates
  > 1. `DEST_freeze1/populationInfo/worldClim/dest.worldclim.csv`
