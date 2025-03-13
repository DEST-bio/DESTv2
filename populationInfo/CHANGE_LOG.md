## Major updates - these are changes to the number of samples, or the names of samples
### 24 August 2024
  • Newest version is now 24August2024
  • Several latitudes and longitudes were incorrect in previous verisons of the DESTv2 meta-data. This led to problems with the automatic identification of provinces and subsequent sample names. Lat & long entries have been manually corrected by the collectors and the sample metadata has been updated.
  • sample metadata fixed: `DESTv2/populationInfo/scripts/update_sampleMetaData_26August2024.R`
  • VCF/GDS fixed: `DESTv2/populationInfo/scripts/update_DataFiles_26August2024.sh`
  • mapping output renamed: `DESTv2/populationInfo/scripts/update_pipelineOutputNames_24August2024.sh`

### 3 May 2024
  • Update the metadata to include standardized baits and country codes
  • Fix problemtic sample name: "IT_Sas_Rec_1_2018-10-18" should be "PR_Sas_Rec_1_2018-10-18"
  • sample metadata fixed: `DESTv2/populationInfo/scripts/update_sampleMetaData_3May2024.R`
  • VCF/GDS fixed: `DESTv2/populationInfo/scripts/update_DataFiles_3May2024.sh`

## Minor updates - these are changes to the metadata
### 13 March 2025
  • The dated version stays the same `24August2024`
  • The updated file is here: `DESTv2/populationInfo/dest_v2.samps_24Aug2024.xa.csv`
  • Added in sex-ratio info. Sex-ratio calculated here: `DESTv2/populationInfo/scripts/update_add_XAratio_24August2024.R`
