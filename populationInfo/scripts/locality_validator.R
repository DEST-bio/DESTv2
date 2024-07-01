### libraries
  library(data.table)
  library(ggplot2)
  library(foreach)
  library(dplyr)
  library(tidygeocoder)
  library(geosphere)
  library(doMC)
  registerDoMC(4)
  library(googlesheets4)
### load metadata
  samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_3May2024.csv")

### geocode
  dist_out <- foreach(i=1:dim(samps)[1], .combine="rbind")%dopar%{
    # i<-586
    message(paste(i, dim(samps)[1], sep=" / "))

    inferred_lat_long <- geocode(samps[i,c("sampleId", "city", "country"), with=F], city="city", country="country", return_input=F, limit=1000)
    inferred_lat_long <- as.data.table(inferred_lat_long)
    inferred_lat_long[,sampleId:=samps[i]$sampleId]
    inferred_lat_long <- merge(samps[,c("sampleId", "lat", "long"),with=F], inferred_lat_long, by="sampleId")
  #foreach(lat.method=c("orig", "invert", "decimal_reduce", "decimal_enlarge"))%do%{
  #  foreach(long.method=c("orig", "invert",  "decimal_reduce", "decimal_enlarge"))%do%{
  #    latlong <- inferred_lat_long[,c("long.x", "lat.x"), with=F]
  #
  #    if(lat.method=="invert") latlong[,lat.x:=lat.x*-1]
  #    if(lat.method=="decimal_reduce") latlong[,lat.x:=lat.x*.1]
  #    if(lat.method=="decimal_enlarge") latlong[,lat.x:=lat.x*10]

  #    if(long.method=="invert") latlong[,long.x:=long.x*-1]
  #    if(long.method=="decimal_reduce")  latlong[,long.x:=long.x*.1]
  #    if(long.method=="decimal_enlarge") latlong[,long.x:=long.x*10]

  #
  #


  #  }
  #}



    inferred_lat_long[,dist_km:=distVincentyEllipsoid(p1=inferred_lat_long[,c("long.x", "lat.x"), with=F], p2=inferred_lat_long[,c("long.y", "lat.y"), with=F])/1000]
    inferred_lat_long[,nMatches:=dim(inferred_lat_long)[1]]
    inferred_lat_long <- inferred_lat_long[which.min(dist_km)]
    return(inferred_lat_long)
  }

  dist_out <- merge(dist_out, samps[,c("sampleId", "set"), with=F], by="sampleId")
  dist_out <- merge(dist_out, samps[,c("sampleId", "collector"), with=F], by="sampleId")

  ggplot(data=dist_out, aes((dist_km))) + geom_histogram(bins=100)



  options(width=200)

  sheet_write(data=dist_out[dist_km>0][order(-dist_km)][set=="DrosEU"], ss="https://docs.google.com/spreadsheets/d/1boyVdV-_XZ7TMSfExc7RUSlHLEHUVv5FgNYnLCgiXuk/edit?gid=0#gid=0", sheet="DrosEU")
  sheet_write(data=dist_out[dist_km>0][order(-dist_km)][set=="dest_plus"], ss="https://docs.google.com/spreadsheets/d/1boyVdV-_XZ7TMSfExc7RUSlHLEHUVv5FgNYnLCgiXuk/edit?gid=0#gid=0", sheet="dest_plus")
  sheet_write(data=dist_out[dist_km>0][order(-dist_km)][set=="DrosEU_3"], ss="https://docs.google.com/spreadsheets/d/1boyVdV-_XZ7TMSfExc7RUSlHLEHUVv5FgNYnLCgiXuk/edit?gid=0#gid=0", sheet="DrosEU_3")
  sheet_write(data=dist_out[dist_km>0][order(-dist_km)][set=="DrosRTEC"], ss="https://docs.google.com/spreadsheets/d/1boyVdV-_XZ7TMSfExc7RUSlHLEHUVv5FgNYnLCgiXuk/edit?gid=0#gid=0", sheet="DrosRTEC")

  dist_out[dist_km>0][order(-dist_km)][set=="dest_plus"]
  dist_out[dist_km>0][order(-dist_km)][set=="DrosEU_3"]
  dist_out[dist_km>0][order(-dist_km)][set=="DrosRTEC"]

  table(dist_out$dist_km>25, dist_out$set)


### reverse geocode
  inferred_city <- reverse_geocode(samps[,c("sampleId", "lat", "long"), with=F], lat="lat", long="long", return_input=T)
  inferred_city <- as.data.table(inferred_city)
  samps <- merge(samps, inferred_city, by="sampleId")

  table(samps$country.x==samps$country.y)
  table(samps$city.x==samps$city.y)

  samps[sampleId=="ES_Ler_Gim_1_2020-10-26",c("city.x", "city.y", "lat.x", "long.x"), with=F]


  samps[samp]
