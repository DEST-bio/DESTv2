### libraries
  library(data.table)
  library(readxl)
  library(sp)
  library(foreach)
  library(lubridate)
  library(ggplot2)
  library(countrycode)
  library(sp)
  library(rworldmap)

### set wd
  setwd("/Users/alanbergland/Documents/GitHub/")

### load DESTv1
  dest_v1 <- fread("DESTv2/populationInfo/samps_10Nov2020.csv")

### load new DrosEU (v3)
  ### this is the collection metadata
    #droseu_sample <- as.data.table(read_excel("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/OriginalMetadata/DrosEUExtraction_metadata_12 Dec 2022.xlsx"))
    droseu_sample <- as.data.table(read_excel("populationInfo/Metadata collection sheets/Finalized docs/Dec 15 2022/JCBN.DrosEUExtraction_metadata_12 Dec 2022.xlsx"))
    setnames(droseu_sample,
            c("SampleID", "Plate number", "Plate position", "Location name"),
            c("sampleId", "Plate_number", "Plate_position", "Location_name"))
    droseu_sample <- droseu_sample[!is.na(sampleId)]
    droseu_sample[,lat:=as.numeric(lat)]
    droseu_sample[,long:=as.numeric(iconv(long, 'utf-8', 'ascii', sub=''))]
    droseu_sample[Altitude=="60 m", Altitude:=60]
    droseu_sample[,altitude:=as.numeric(iconv(Altitude, 'utf-8', 'ascii', sub=''))]

    #droseu_sample[is.na(as.numeric(as.character(droseu_sample$lat)))]
    #droseu_sample[is.na(as.numeric(as.character(droseu_sample$long)))]

    droseu_sample[Location_name=="Charlottesville", lat:=37.9790]
    droseu_sample[Location_name=="Charlottesville", long:=-78.4897]

  ### this is the DNA_library metadata
    #droseu_lib <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/Final Extraction data  _Microgen_DrosEU_2017-2021.xlsx"))
    droseu_lib <- as.data.table(read_excel("populationInfo/2022_DESTv2/Metadata collection sheets/Finalized docs/Dec 15 2022/AOB.JCBN.Final Extraction data _Microgen_DrosEU_2017-2021.xlsx"))

    setnames(droseu_lib,
             c("SampleID", "Sample ID**\r\n(<15 letters)"),
             c("sampleId", "SequencingId"))

  ### First batch of sequencing includes all samples except these:
    low_qual <- c(32, 50, 76, 77, 83, 87, 92, 96, 107, 120, 121, 127,136, 139, 144, 158, 185, 205, 210, 215, 233, 235, 236, 241, 243, 244, 246, 247, 249, 250, 252, 253)

  ### merge droseu
    droseu.full <- merge(droseu_sample, droseu_lib, by="sampleId", all=T)
    table(droseu.full$"Location_name" == droseu$"Location name")
    as.data.frame(droseu.full[,c("Location_name", "Location name"), with=F])

  ### clean up a bit
    droseu <- droseu.full[,
                          c("sampleId", "SequencingId",
                            "Location_name", "Country.x", "Country.y",
                            "Date.x", "Day_apprx", "Month", "Year", "Date_Exact", "Date_Comp",
                            "lat", "long", "altitude",
                            "Extraction",
                            "Sampling strategy", "Wild/F1", "Fruit type"), with=F]


  ### which "locality" does each new sample belong to?
    dest_locality <- na.omit(dest_v1[,list(lat=mean(lat), long=mean(long), lat_var=var(lat), long_var=var(long)),
                                      list(city_dest=city,
                                           country_dest=country,
                                           locality_dest=locality,
                                           continent)])

    locs <- foreach(i=1:dim(droseu)[1], .combine="rbind")%do%{
      # i <- 76

      if(!is.na(droseu[i]$lat)) {
        dists <- spDistsN1(pts=as.matrix(dest_locality[,c("long", "lat"), with=F]),
                            pt=as.matrix(droseu[i, c("long", "lat"), with=F]), longlat=T)

        o <- data.table(sampleId=droseu[i]$sampleId,
                    locality_dest=dest_locality[which.min(dists)]$locality_dest,
                    continent_dest=dest_locality[which.min(dists)]$continent,
                    city_dest=dest_locality[which.min(dists)]$city_dest,
                    country_dest=dest_locality[which.min(dists)]$country_dest,
                    locality_dist=min(dists),
                    lat_dest=dest_locality[which.min(dists)]$lat,
                    long_dest=dest_locality[which.min(dists)]$long)
      } else {
        o <- data.table(sampleId=droseu[i]$sampleId,
                    locality_dest=NA,
                    continent_dest=NA,
                    city_dest=NA,
                    country_dest= NA,
                    locality_dist= NA, lat_dest=NA, long_dest=NA)

      }
      return(o)
    }

    droseu <- merge(droseu, locs, by="sampleId")

  ### some double checks:
    ggplot(data=droseu, aes(locality_dist)) + geom_histogram()

    ggplot(data=droseu, aes(y=lat, x=long)) + geom_point()

    ### does the country name from the Microgen and the Metadata sheets align? YES
    table(droseu$Country.x==droseu$Country.y); droseu[droseu$Country.x!=droseu$Country.y] ### these are just some simple discrepancies beteween names; Believe Country.x

    ### does country name from Metadata with the nearest inferred locality; yes, except for newly defined localities
    table(droseu$Country.x==droseu$country_dest);
    droseu[droseu$Country.x!=droseu$country_dest, c("Country.x", "country_dest", "Location_name", "locality_dest", "city_dest", "locality_dist"), with=F]

### Locality name generation
  droseu[,countryCode:=countrycode(Country.x, "country.name", "iso2c")]
  droseu[city_dest=="Charlottesville", countryCode:="VA"]
  droseu[continent_dest!="NorthAmerica",
        new_locality:=paste(countryCode, substr(Location_name, 0, 3), sep="_")]
  droseu[continent_dest=="NorthAmerica",
        new_locality:=paste(countryCode, tolower(substr(Location_name, 0, 2)), sep="_")]

  ggplot(droseu[Country.x!="Australia"], aes(y=as.numeric(new_locality==locality_dest), x=locality_dist)) +
      geom_point(size = 2) +
      geom_smooth(method = "glm", , se = F,
          method.args = list(family = "binomial"))

  droseu[droseu$new_locality==droseu$locality_dest][locality_dist>5][,
          c("Location_name", "Country.x","new_locality", "lat", "lat_dest", "long", "long_dest", "locality_dist"), with=F][order(locality_dist)]

  droseu[,sampleId_use:=paste(new_locality, Date_Comp, sep="_")]
  droseu[,sampleId_use:=gsub("-", "_", sampleId_use)]

### date stuff
  droseu[is.na(Month), date:=paste(Year, NA, NA, se="-")]
  droseu[,date:=ymd(paste(Year, Month, Day_apprx, sep="-"))]
  droseu[,yday:=yday(date)]
  droseu[,year:=Year]

## other stuff
  droseu[,type:="pooled"]
  droseu[,set:="DrosEU_3"]

### rename new droseu data
  setnames(droseu,
           c("sampleId", "continent_dest", "new_locality", "Sampling strategy", "Wild/F1", "sampleId_use", "Country.x"),
           c("sampleId_internal", "continent", "locality", "Sampling_strategy", "Wild_F1", "sampleId", "country"))

  setnames(droseu, "Location_name", "city")

### Cville data
  cville <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/TableS1.Sample Metadata.xlsx"))
  cville <- cville[source_data=="This Study"][Seq_strategy=="pooled"]
  cville[,yday:=yday(ymd(collectionDate))]
  cville[,year:=year(ymd(collectionDate))]
  cville[,continent:="North_America"]
  cville[, lat:=37.9790]
  cville[, long:=-78.4897]
  cville[,type:="pooled"]
  cville[,set:="cville"]
  cville[,country:="US"]

### merge together
  dest_v2 <- rbindlist(list(
    dest_v1[set%in%c("DrosEU", "DrosRTEC", "dgn")][,c("sampleId", "year", "yday", "locality", "lat", "long", "continent", "country", "city", "type", "set"), with=F],
    droseu[,c("sampleId", "year", "yday", "locality", "lat", "long", "continent", "country", "city", "type", "set"), with=F],
    cville[,c("sampleId", "year", "yday", "locality", "lat", "long", "continent", "country", "city", "type", "set"), with=F]), fill=T)

  dest_v2[continent=="NorthAmerica", continent:="North_America"]

### fix continents
  dest_v2[lat<20 & continent=="Europe", continent:="North_America"]
  dest_v2[country=="Morocco", continent:="Africa"]
  dest_v2[lat< -30, continent:="Australia"]

### output
  write.csv(dest_v2, "DESTv2/populationInfo/dest_v2_samples.csv", quote=F)

### basic plot
  sampPlot <- ggplot(dest_v2[type=="pooled"]) +
              geom_line(aes(x=yday, y=lat, group=locality), linetype="dashed") +
              geom_point(aes(x=yday, y=lat, color=continent), size=.75) +
              facet_grid(~year) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  sampPlot



  nperLocale <- dest_v2[set!="dgn",list(.N,
                              nyears=length(unique(year)),
                              ave_per_year=median(prop.table(table(year)))*length(year),
                              max_per_year=max(prop.table(table(year)))*length(year)),
                        list(locality)]

  ave_countPlot <- ggplot(data=nperLocale, aes(x=nyears, y=ave_per_year)) +
  geom_jitter() +
  geom_text(data=nperLocale[nyears>7 & ave_per_year>2], aes(x=nyears, y=ave_per_year, label=locality)) +
  ylab("Median number of\nsamples per year") + xlab("Total number of years")

  max_countPlot <- ggplot(data=nperLocale, aes(x=nyears, y=max_per_year)) +
  geom_jitter() +
  geom_text(data=nperLocale[max_per_year>5], aes(x=nyears, y=max_per_year, label=locality)) +
  ylab("Max number of\nsamples per year") + xlab("Total number of years")

  countPlot <- ggplot(data=nperLocale, aes(x=nyears, y=N)) +
  geom_jitter() +
  geom_text(data=nperLocale[N>10], aes(x=nyears, y=N, label=locality)) +
  ylab("Total number of \nsamples per year") + xlab("Total number of years")


### plot
  layout <- "
  AAA
  BCD"

  mega <-
  sampPlot +
  ave_countPlot + max_countPlot + countPlot +
  plot_layout(design=layout) + plot_annotation(tag_level="A")

  ggsave(mega, file="~/mega.png", height=8, w=16)
