### libraries
  library(data.table)
  library(readxl)
  library(sp)
  library(foreach)

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
    dest_locality <- dest_v1[,list(lat=mean(lat), long=mean(long)), list(locality)]

    locs <- foreach(i=1:dim(droseu)[1], .combine="rbind")%do%{
      # i <- 76

      if(!is.na(droseu[i]$lat)) {
        dists <- spDistsN1(pts=as.matrix(na.omit(dest_locality[,c("long", "lat"), with=F])),
                            pt=as.matrix(droseu[i, c("long", "lat"), with=F]))

        o <- data.table(sampleid=droseu[i]$sampleid,
                    locality=dest_locality[which.min(dists)]$locality,
                    locality_dist=min(dists))
      } else {
        o <- data.table(sampleid=droseu[i]$sampleid,
                    locality=NA,
                    locality_dist=NA)

      }
      return(o)
    }

    droseu <- merge(droseu, )


### Cville data
  cville <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/TableS1.Sample Metadata.xlsx"))
  cville <- cville[source_data=="This Study"][Seq_strategy=="pooled"]
