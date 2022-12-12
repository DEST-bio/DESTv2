### libraries
  library(data.table)
  library(readxl)
  library(sp)
  library(foreach)
  library(lubridate)
  library(ggplot2)

### set wd
  setwd("/Users/alanbergland/Documents/GitHub/")

### load DESTv1
  dest_v1 <- fread("DESTv2/populationInfo/samps_10Nov2020.csv")

### load new DrosEU (v3)
  ### this is the collection metadata
    droseu_sample <- as.data.table(read_excel("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/OriginalMetadata/DrosEUExtraction_metadata_12 Dec 2022.xlsx"))
    setnames(droseu_sample,
            c("Sample ID", "Plate number", "Plate position", "Location name"),
            c("sampleid", "Plate_number", "Plate_position", "Location_name"))
    droseu_sample <- droseu_sample[!is.na(sampleid)]
    droseu_sample[,lat:=as.numeric(lat)]
    droseu_sample[,long:=as.numeric(long)]

    #droseu_sample[is.na(as.numeric(as.character(droseu_sample$lat)))]
    #droseu_sample[is.na(as.numeric(as.character(droseu_sample$long)))]

    droseu_sample[is.na(long) & !is.na(lat), long:=c(41.15,35.1467,35.1467,56.3204,39.1689,55.9323,55.9324)]
    droseu_sample[Location_name=="Charlottesville", lat:=37.9790]
    droseu_sample[Location_name=="Charlottesville", long:=-78.4897]

  ### this is the DNA_library metadata
    droseu_lib <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/Final Extraction data  _Microgen_DrosEU_2017-2021.xlsx"))
    setnames(droseu_lib,
             c("...1", "Location name"),
             c("sampleid", "Location_name"))

  ### First batch of sequencing includes all samples except these:
    low_qual <- c(32, 50, 76, 77, 83, 87, 92, 96, 107, 120, 121, 127,136, 139, 144, 158, 185, 205, 210, 215, 233, 235, 236, 241, 243, 244, 246, 247, 249, 250, 252, 253)

  ### merge droseu
    droseu <- merge(droseu_sample, droseu_lib, by="sampleid")
    table(droseu$"Location_name.x" == droseu$"Location name.y")
    as.data.frame(droseu[,c("Location_name.x", "Location_name.y"), with=F])

  ### which "locality" does each new sample belong to?
    dest_locality <- na.omit(dest_v1[,list(lat=mean(lat), long=mean(long)), list(locality, continent)])


    locs <- foreach(i=1:dim(droseu)[1], .combine="rbind")%do%{
      # i <- 250

      if(!is.na(droseu[i]$lat)) {
        dists <- spDistsN1(pts=as.matrix(dest_locality[,c("long", "lat"), with=F]),
                            pt=as.matrix(droseu[i, c("long", "lat"), with=F]), longlat=T)

        o <- data.table(sampleid=droseu[i]$sampleid,
                    locality=dest_locality[which.min(dists)]$locality,
                    continent=dest_locality[which.min(dists)]$continent,

                    locality_dist=min(dists))
      } else {
        o <- data.table(sampleid=droseu[i]$sampleid,
                    locality=NA, continent=NA,
                    locality_dist=NA)

      }
      return(o)
    }

    droseu <- merge(droseu, locs, by="sampleid")
    droseu[,date:=ymd(paste(Year, Month, Date.x, sep="-"))]
    droseu[,yday:=yday(date)]
    setnames(droseu, "sampleid", "sampleId")

### Cville data
  cville <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/TableS1.Sample Metadata.xlsx"))
  cville <- cville[source_data=="This Study"][Seq_strategy=="pooled"]
  cville[,yday:=yday(ymd(collectionDate))]
  cville[,year:=year(ymd(collectionDate))]
  cville[,continent:="North_America"]
  cville[, lat:=37.9790]
  cville[, long:=-78.4897]

### merge together
  dest_v2 <- rbindlist(list(
    dest_v1[set%in%c("DrosEU", "DrosRTEC")][,c("sampleId", "year", "yday", "locality", "lat", "long", "continent"), with=F],
    droseu[,c("sampleId", "year", "yday", "locality", "lat", "long", "continent"), with=F],
    cville[,c("sampleId", "year", "yday", "locality", "lat", "long", "continent"), with=F]))

  dest_v2[continent=="NorthAmerica", continent:="North_America"]


### output
  write.csv(dest_v2, )

### basic plot
  sampPlot <- ggplot(dest_v2[year>=2009]) +
              geom_line(aes(x=yday, y=lat, group=locality), linetype="dashed") +
              geom_point(aes(x=yday, y=lat, color=continent)) +
              facet_grid(~year)
sampPlot
