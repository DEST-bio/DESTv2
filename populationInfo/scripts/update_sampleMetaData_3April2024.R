### goals:
### • incorporate standardized bait/substrate labels
### • standardize some country codes (US & UK), seq type
### • correct sample idetity: "IT_Sas_Rec_1_2018-10-18" should be "PR_Sas_Rec_1_2018-10-18" and the longitude should be -8.51 (and not 8.41)... province should be "Porto" (not "Sassari") (edited)


### libraries
  library(data.table)

### load in old data
  samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/old_versions/dest_v2.samps_8Jun2023.csv")

### load in standardized bait data
  sb <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/OriginalMetadata/dest_v2.samps_8Jun2023_unifyBaits_v1.tab")
  table(is.na(sb$fruit_type_Currated))
  table(sb$fruit_type)
  table(is.na(sb$Boiled))
  setnames(sb, "fruit_type_Currated", "fruit_type_curated")
  samps <- merge(samps, sb[,c("sampleId", "fruit_type_curated"), with=F])

  str(samps)
