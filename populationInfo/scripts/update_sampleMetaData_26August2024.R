### libraries
  library(googlesheets4)
  library(data.table)
  library(countrycode)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(sf)
  library(stringi)

### load in manually corrected data
  droseu <- as.data.table(read_sheet(ss="https://docs.google.com/spreadsheets/d/1boyVdV-_XZ7TMSfExc7RUSlHLEHUVv5FgNYnLCgiXuk/edit?gid=426327939#gid=426327939", sheet="DrosEU"))
  dest_plus <- as.data.table(read_sheet(ss="https://docs.google.com/spreadsheets/d/1boyVdV-_XZ7TMSfExc7RUSlHLEHUVv5FgNYnLCgiXuk/edit?gid=426327939#gid=426327939", sheet="dest_plus"))
  DrosEU_3 <- as.data.table(read_sheet(ss="https://docs.google.com/spreadsheets/d/1boyVdV-_XZ7TMSfExc7RUSlHLEHUVv5FgNYnLCgiXuk/edit?gid=426327939#gid=426327939", sheet="DrosEU_3"))
  DrosRTEC <- as.data.table(read_sheet(ss="https://docs.google.com/spreadsheets/d/1boyVdV-_XZ7TMSfExc7RUSlHLEHUVv5FgNYnLCgiXuk/edit?gid=426327939#gid=426327939", sheet="DrosRTEC"))
  cville <- as.data.table(read_sheet(ss="https://docs.google.com/spreadsheets/d/1boyVdV-_XZ7TMSfExc7RUSlHLEHUVv5FgNYnLCgiXuk/edit?gid=426327939#gid=426327939", sheet="cville"))

  corrected <- rbindlist(list(droseu, dest_plus, DrosEU_3, DrosRTEC, cville), fill=T)

  corrected[is.na(correct_long), correct_long:=reported_long]
  corrected[is.na(correct_lat), correct_lat:=reported_lat]
  corrected[is.na(correct_city), correct_city:=city]
  corrected

### get new political identifiers
  states <- ne_states(returnclass="sf")
  countries <- ne_countries(returnclass="sf")

  DT <- corrected[,c("sampleId", "correct_long", "correct_lat"), with=F]
  setnames(DT, c("correct_long", "correct_lat"), c("longitude", "latitude"))
  DT_sf = st_as_sf(DT, coords = c("longitude", "latitude"),
                crs = 4326, agr = "constant")

  sf_use_s2(FALSE)
  dest_states <- cbind(DT, as.data.table(states)[st_nearest_feature(DT_sf, states)])
  dest_countries <- cbind(DT, as.data.table(countries)[st_nearest_feature(DT_sf, countries)])

  corrected <- merge(corrected, dest_states[,c("sampleId", "iso_a2", "name"), with=F], by="sampleId", all.x=T)
  corrected <- merge(corrected, dest_countries[,c("sampleId", "continent"), with=F], by="sampleId", all.x=T)

  setnames(corrected, "name", "province")
  corrected[,province:=stri_trans_general(str = province,
                                      id = "Latin-ASCII")]

  corrected[province=="L'viv",province:="Lviv"]
  corrected[,continent:=gsub(" ", "_", continent)]
  corrected[,province:=gsub(" ", "_", province)]
  corrected[,locality:=paste(iso_a2, substr(province, 0, 3), substr(city, 0, 3), sep="_")]
  corrected[,ss:=tstrsplit(sampleId, "_")[[4]]]
  corrected[,date:=tstrsplit(sampleId, "_")[[5]]]


  corrected[,new_sampleId:=paste(locality, ss, date, sep="_")]
  corrected[sampleId!=new_sampleId]
  setnames(corrected, "province", "correct_province")
  setnames(corrected, "country", "correct_country")
  setnames(corrected, "locality", "correct_locality")
  setnames(corrected, "new_sampleId", "correct_sampleId")
  corrected_small <- corrected[,c("sampleId", "correct_sampleId", "correct_locality", "correct_lat", "correct_long", "correct_city", "correct_province", "correct_country" ), with=F]

### merge
  samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/old_versions/dest_v2.samps_3May2024.csv")
  sm <- merge(corrected_small, samps, by="sampleId", all.y=T)
  sm[is.na(correct_sampleId), correct_sampleId:=sampleId]
  sm[is.na(correct_locality), correct_locality:=locality]
  sm[is.na(correct_lat), correct_lat:=lat]
  sm[is.na(correct_long), correct_long:=long]
  sm[is.na(correct_city), correct_city:=city]
  sm[is.na(correct_province), correct_province:=province]
  sm[is.na(correct_country), correct_country:=country]

  sm[correct_sampleId!=sampleId]

  smo <- sm[order(correct_sampleId==sampleId)]
  write_sheet(smo, ss="https://docs.google.com/spreadsheets/d/1boyVdV-_XZ7TMSfExc7RUSlHLEHUVv5FgNYnLCgiXuk/edit?gid=2051031941#gid=2051031941", sheet="UPDATED METADATA")

### clean up
  names(sm)
  setnames(sm,
          c("sampleId", "locality", "lat", "long", "city", "province", "country"),
          paste("old", c("sampleId", "locality", "lat", "long", "city", "province", "country"), sep="_"))

  setnames(sm,
          paste("correct", c("sampleId", "locality", "lat", "long", "city", "province", "country"),sep="_"),
          c("sampleId", "locality", "lat", "long", "city", "province", "country"))

  sm[sampleId!=old_sampleId]


### samps
  samps_new <- sm[,-paste("old", c("sampleId", "locality", "lat", "long", "city", "province", "country"), sep="_")]
  samps_new
  write.csv(samps_new, file="/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_26Aug2024.csv", quote=F, row.names=F)
  conversion <- sm[sampleId!=old_sampleId,c("sampleId", "old_sampleId")]

  write.csv(conversion, file="/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/scripts/update_DataFiles_26August2024_conversionTable.csv", quote=F, row.names=F)
