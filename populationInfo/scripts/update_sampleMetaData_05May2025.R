### notes
  ## 1. Switch to Ukranian spelling of some province names. done
  ## 2. Differentiate DrosEU1 & DrosEU2 from DrosEU set
  ## 3. update SRA info for new DrosEU samples
  ### we are keeping the same dated version because the sample names haven't changed

### libraries
  library(data.table)

### load in current version
  samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.xa.csv")

### change province names in Ukraine to Ukrainian spelling
  samps[city=="Kyiv", province:="Kyiv_City"]
  samps[city=="Odesa", province:="Odesa"]

### load in old sample names
  ct <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/scripts/update_DataFiles_24August2024_conversionTable.csv")
  ct <- rbindlist(list(ct,
                        data.table(old_sampleId="IT_Sas_Rec_1_2018-10-18", sampleId="PT_Por_Rec_1_2018-10-18"),
                        data.table(old_sampleId="ES_Bal_Tom_1_2021-10-07", sampleId="ES_Ciu_Tom_1_2021-10-07")),
                  use.names=T)

  samps <- merge(samps, ct, by="sampleId", all.x=T)
  samps[is.na(old_sampleId), old_sampleId:=sampleId]
  table(samps$sampleId==samps$old_sampleId)

### load in SRA info
  sra <- fread("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/OriginalMetadata/table_download.droseu3_sra.tsv")
  str(sra)
  sra <- sra[,c("Title", "Accession")]
  sra[,old_sampleId:=tstrsplit(Title, " ")[[1]]]
  sra <- sra[,c("old_sampleId", "Accession")]

  oldSampleId_newSRAs <- merge(samps[,c("Recommendation", "set", "sampleId", "old_sampleId", "SRA_Accession")], sra, by.x="old_sampleId", by.y="old_sampleId", all=T)

  ### merging on the old sample names
    ### which samples are in the samps table but do not have an accession? There are 10 of them
    table(oldSampleId_newSRAs$set, is.na(oldSampleId_newSRAs$Accession), oldSampleId_newSRAs$Recommendation)
    table(oldSampleId_newSRAs$set, is.na(oldSampleId_newSRAs$Accession))

    oldSampleId_newSRAs[set%like%"DrosEU_3" & is.na(Accession)]

    ### which samples are NOT in the samps but have an accession? There is one
    table(is.na(oldSampleId_newSRAs$set), is.na(oldSampleId_newSRAs$Accession))
    oldSampleId_newSRAs[is.na(set) & !is.na(Accession)]


  ### tack together
    oldSampleId_newSRAs[is.na(SRA_Accession)]
    oldSampleId_newSRAs[is.na(SRA_Accession), SRA_Accession:=Accession]
    table(oldSampleId_newSRAs$set, is.na(oldSampleId_newSRAs$SRA_Accession))

    oldSampleId_newSRAs[set%like%"DrosEU_3" & is.na(SRA_Accession)]
    oldSampleId_newSRAs[is.na(set) & !is.na(SRA_Accession)]

    str(samps)
    samps2 <- merge(samps[,-c("old_sampleId")], oldSampleId_newSRAs[,c("sampleId", "SRA_Accession")], by="sampleId", all.x=T)
    table(is.na(samps2$SRA_Accession.x), is.na(samps2$SRA_Accession.y))
    table(is.na(samps2$SRA_Accession.y), samps2$set)

    samps2 <- samps2[,-c("SRA_Accession.x")]
    setnames(samps2, "SRA_Accession.y", "SRA_Accession")

    samps2[is.na(SRA_Accession) & set!="dgn",c("sampleId", "Recommendation", "library_result")]

    table(samps2$Recommendation, samps2$library_result)

### write
  write.csv(samps2, quote=F, row.names=F, file="~/dest_v2_samps_24Aug2024.xa.sra.csv")


table(samps$set)
dim(sra)




  table(newSampleId_newSRAs$set, is.na(newSampleId_newSRAs$Accession))

  oldSampleId_newSRAs[is.na(sampleId)]
  newSampleId_newSRAs[is.na(Accession)]


  newSRAs[is.na(Accession) & set=="DrosEU_3"]$sampleId
  table(newSRAs$set)
  dim(sra)

  newSRAs[is.na(SRA_Accession) & !is.na(Accession), SRA_Accession:=Accession]
  table(is.na(newSRAs$SRA_Accession))

  #samps_new <- newSRAs[,-c("Accession")]
  table(is.na(samps_new$SRA_Accession), samps_new$set)

  samps_new[is.na(SRA_Accession) & set=="DrosEU_3"]$sampleId


### set
  table(samps_new[set%like%"DrosEU"]$year, samps_new[set%like%"DrosEU"]$set)
  samps_new[year==2014 & set=="DrosEU", set:="DrosEU_1"]
  samps_new[year>2014 & set=="DrosEU", set:="DrosEU_2"]
