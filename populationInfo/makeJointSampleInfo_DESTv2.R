### libraries
  library(data.table)
  library(readxl)
  library(foreach)
  library(lubridate)
  library(ggplot2)
  library(countrycode)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(sf)
  library(stringi)
  #library(rworldmap)

### set wd
  setwd("/Users/alanbergland/Documents/GitHub/")

### load DESTv1
  dest_v1 <- fread("DESTv2/populationInfo/OriginalMetadata/DEST_v1.samps_10Nov2020.csv")

#  ### add in "exact date" info from Machado et al
#    drosRTEC <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/elife-67577-supp1-v2.xlsx"))
#    drosRTEC[,day:=as.numeric(as.character(Day))]
#    drosRTEC[,exactDate:=!is.na(day)]
#    drosRTEC[,c("Day", "Month", "Year", "exactDate"), with=F]
#    setnames(drosRTEC, "Sample", "sampleId")
#    dim(dest_v1)
#    dest_v1 <- merge(dest_v1, drosRTEC[,c("sampleId", "exactDate"), with=F], all.x=T)
#    dim(dest_v1) ### should have 272 rows

#  ### get year month date info
#    ### for DrosRTEC/DrosEU
#      dest_v1[set%in%c("DrosRTEC", "DrosEU"),year:= as.numeric(paste("20", tstrsplit(collectionDate, "/")[[3]], sep=""))]
#      dest_v1[set%in%c("DrosRTEC", "DrosEU"),month:=as.numeric(tstrsplit(collectionDate, "/")[[1]])]
#      dest_v1[set%in%c("DrosRTEC", "DrosEU"),day:=  as.numeric(tstrsplit(collectionDate, "/")[[2]])]
#      dest_v1[is.na(exactDate) & is.na(day), exactDate:=F]
#      dest_v1[exactDate==F & is.na(day), day:=15]
#      dest_v1[sampleId=="NC_ra_03_n", year:=2003]
#      dest_v1[sampleId=="NC_ra_03_n", month:=NA]
#
#    ### for DGN
#      dest_v1[set%in%c("dgn"),year:= as.numeric(tstrsplit(collectionDate, "/")[[1]])]
#      dest_v1[set%in%c("dgn"),month:=as.numeric(tstrsplit(collectionDate, "/")[[2]])]
#
#    ### padd the missing info
#      dest_v1[is.na(exactDate) & is.na(month), exactDate:=F]
#      dest_v1[exactDate==F & is.na(month), month:=6]
#
#      dest_v1[is.na(exactDate) & is.na(year), exactDate:=F]
#
#      dest_v1[is.na(exactDate), exactDate:=T]
#
#    dest_v1[,c("exactDate", "year", "month", "day"), with=F]
#
#    dest_v1[,yday:=yday(ymd(paste(year, month, day, sep="-")))]
#
#    dest_v1[is.na(yday)]
#    table(dest_v1$exactDate)



  ### some naming conventions
    setnames(dest_v1, "type", "library_type")
    setnames(dest_v1, "sampleType", "fly_type")

    setnames(dest_v1, "Model", "seq_platform")
    setnames(dest_v1, "SRA_accession", "SRA_Accession")
    setnames(dest_v1, "collectionDate", "date_orig")
    dest_v1[,reference:="https://doi.org/10.1093/molbev/msab259"]
    dest_v1[,sampling_strategy:=NA]

### load new DrosEU (v3)
  ### this is the collection metadata
    #droseu_sample <- as.data.table(read_excel("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/OriginalMetadata/DrosEUExtraction_metadata_12 Dec 2022.xlsx"))
    droseu_sample <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/AOB.JCBN.DrosEUExtraction_metadata_12 Jan 2023.xlsx"))
    setnames(droseu_sample,
            c("SampleID", "Plate number", "Plate position"),
            c("sampleId", "Plate_number", "Plate_position"))
    droseu_sample <- droseu_sample[!is.na(sampleId)]
    droseu_sample[,lat:=as.numeric(lat)]
    droseu_sample[,long:=as.numeric(iconv(long, 'utf-8', 'ascii', sub=''))]
    droseu_sample[Altitude=="60 m", Altitude:=60]
    droseu_sample[,altitude:=as.numeric(iconv(Altitude, 'utf-8', 'ascii', sub=''))]

    #droseu_sample[is.na(as.numeric(as.character(droseu_sample$lat)))]
    #droseu_sample[is.na(as.numeric(as.character(droseu_sample$long)))]

    droseu_sample[city=="Charlottesville", lat:=37.9790]
    droseu_sample[city=="Charlottesville", long:=-78.4897]

  ### this is the DNA_library metadata
    #droseu_lib <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/Final Extraction data  _Microgen_DrosEU_2017-2021.xlsx"))
    droseu_lib <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/AOB.JCBN.Final Extraction data _Microgen_DrosEU_2017-2021.xlsx"))

    setnames(droseu_lib,
             c("SampleID", "Sample ID**\r\n(<15 letters)"),
             c("sampleId", "SequencingId"))

    setnames(droseu_lib, "Extraction", "nFlies")

  ### First batch of sequencing includes all samples except these:
    low_qual <- c(32, 50, 76, 77, 83, 87, 92, 96, 107, 120, 121, 127,136, 139, 144, 158, 185, 205, 210, 215, 233, 235, 236, 241, 243, 244, 246, 247, 249, 250, 252, 253)
    low_qual <- paste("DrosEu-", low_qual, sep="")
    #samps[SequencingId%in%low_qual]

  ### merge droseu
    droseu.full <- merge(droseu_sample, droseu_lib, by="sampleId", all=T)

  ### clean up a bit
    setnames(droseu.full, "Year", "year")
    droseu <- droseu.full[,
                          c("sampleId", "SequencingId",
                            "city_orig", "city", "Country.x", "Country.y",
                            "min_day", "max_day", "min_month", "max_month", "year",
                            "lat", "long", "altitude",
                            "nFlies",   "bio_rep", "tech_rep", "exp_rep", "loc_rep",
                            "Sampling strategy", "Wild/F1", "Fruit type", "Collectors name.x"), with=F]

  ### some double checks:
    ggplot(data=droseu, aes(locality_dist)) + geom_histogram()

    ggplot(data=droseu, aes(y=lat, x=long)) + geom_point()

    ### does the country name from the Microgen and the Metadata sheets align? YES
    table(droseu$Country.x==droseu$Country.y); droseu[droseu$Country.x!=droseu$Country.y] ### these are just some simple discrepancies beteween names; Believe Country.x

## other stuff
  droseu[,type:="pooled"]
  droseu[,set:="DrosEU_3"]

### rename new droseu data
  setnames(droseu,
           c("sampleId", "Sampling strategy", "Wild/F1", "Country.x"),
           c("sampleId_internal", "sampling_strategy", "Wild_F1", "country"))

  setnames(droseu, "type", "library_type")
  setnames(droseu, "Wild_F1", "fly_type")
  setnames(droseu, "Collectors name.x", "collector")
  setnames(droseu, "Fruit type", "fruit_type")

### Cville data
  cville <- as.data.table(read_excel("DESTv2/populationInfo/OriginalMetadata/TableS1.Sample Metadata.xlsx"))
  cville <- cville[source_data=="This Study"][Seq_strategy=="pooled"]
  cville[,yday:=yday(ymd(collectionDate))]
  cville[,year:=year(ymd(collectionDate))]
  cville[,continent:="North_America"]
  cville[, lat:=37.9790]
  cville[, long:=-78.4897]
  cville[,library_type:="pooled"]
  cville[,set:="cville"]

  setnames(cville, "sample_type", "fly_type")
  setnames(cville, "Seq_strategy", "library_type")
  cville[,sampling_strategy:="Sweep Netting + Aspiration"]
  cville[,fruit_type:=loc_rep]

### Other published data
  pub <- as.data.table(read_excel("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/OriginalMetadata/DESTplus_Metadata.xlsx"))
  pub[,Date.orig:=paste(year, Month, Day, sep="-")]

  pub[,year:=as.numeric(as.character(year))]

  pub[,set:="dest_plus"]
  setnames(pub, "sampleType", "fly_type")
  setnames(pub, "Sequencing Platform", "seq_platform")
  setnames(pub, "type", "library_type")
  setnames(pub, "Collector", "collector")
  setnames(pub, "SRA_Accesion", "SRA_Accession")
  setnames(pub, "REFERENCE", "reference")

### merge together

  commonCols <- c("sampleId", "set", "SequencingId",
                  "year", "min_day", "max_day", "min_month", "max_month",
                  "locality", "lat", "long",
                  "country", "city",
                  "bio_rep", "tech_rep", "exp_rep", "loc_rep", "fruit_type",
                  "nFlies", "fly_type", "library_type", "sampling_strategy", "seq_platform",
                  "collector", "SRA_Accession", "reference")

  dest_v2 <- rbindlist(list(dest_v1[,commonCols[commonCols%in%names(dest_v1)], with=F],
                            droseu[,commonCols[commonCols%in%names(droseu)], with=F],
                            cville[,commonCols[commonCols%in%names(cville)], with=F],
                            pub[,commonCols[commonCols%in%names(pub)], with=F]),
                          fill=T)
  dest_v2[,tmpId:=c(1:dim(dest_v2)[1])]

### fix a few problematic longitudes
  dest_v2[city=="Muthill", long:=ifelse(long>0, -1*long, long)]
  dest_v2[city=="Edinburgh", long:=ifelse(long>0, -1*long, long)]

### get political identifiers
  states <- ne_states(returnclass="sf")
  countries <- ne_countries(returnclass="sf")

  DT <- dest_v2[!is.na(lat) & !is.na(long),c("tmpId", "long", "lat"), with=F]
  setnames(DT, c("long", "lat"), c("longitude", "latitude"))
  DT_sf = st_as_sf(DT, coords = c("longitude", "latitude"),
                crs = 4326, agr = "constant")

  sf_use_s2(FALSE)
  dest_states <- cbind(DT, as.data.table(states)[st_nearest_feature(DT_sf, states)])
  dest_countries <- cbind(DT, as.data.table(countries)[st_nearest_feature(DT_sf, countries)])


  dest_v2 <- merge(dest_v2, dest_states[,c("tmpId", "iso_a2", "name"), with=F], by="tmpId", all.x=T)
  dest_v2 <- merge(dest_v2, dest_countries[,c("tmpId", "continent"), with=F], by="tmpId", all.x=T)

  setnames(dest_v2, "name", "province")
  dest_v2[,province:=stri_trans_general(str = province,
                                      id = "Latin-ASCII")]


  dest_v2[,continent:=gsub(" ", "_", continent)]

### construct locality ID
  setnames(dest_v2, "locality", "locality_old")
  dest_v2[city=="Ithaca NY", city:="Ithaca"]
  dest_v2[city=="Karensminde Orchard", city:="Karensminde"]
  dest_v2[city=="Bowdoin" & is.na(loc_rep), loc_rep:="Rocky Ridge Orchard"]
  dest_v2[city=="Iraklio", city:="Irakilo"]
  dest_v2[city=="Mariupol-2", loc_rep:="Mariupol-2"]
  dest_v2[city=="Mariupol-1", loc_rep:="Mariupol-1"]
  dest_v2[city=="Mariupol-2", city:="Mariupol"]
  dest_v2[city=="Mariupol-1", city:="Mariupol"]

  dest_v2[city=="Raleigh" & set=="dgn", exp_rep:="insilico pool"]
  dest_v2[city=="Raleigh" & set=="DrosRTEC", exp_rep:="pool-seq"]


  dest_v2[,locality:=paste(iso_a2, substr(province, 0, 3), substr(city, 0, 3), sep="_")]

  tmp <- dest_v2[,list(nCities=length(unique(city)), vec=paste(unique(city), collapse=";")), list(locality)]
  tmp[nCities>1]

### make date string

  dest_v2[,min_jday:=yday(ymd(paste(year, min_month, min_day, sep="-")))]
  dest_v2[,max_jday:=yday(ymd(paste(year, max_month, max_day, sep="-")))]

  dest_v2[is.na(max_jday)]

  dest_v2[,jday:=round(min_jday/2 + max_jday/2)]
  dest_v2[,exactDate:=jday==min_jday]

  dest_v2[,min_month:=month(make_date(year) + min_jday - 1)]
  dest_v2[,max_month:=month(make_date(year) + max_jday - 1)]


  table(dest_v2$exactDate, dest_v2$set)

  dest_v2[,date_string:=as.character(make_date(year) + jday - 1)]
  dest_v2[is.na(jday)]
  dest_v2[,pre:=paste(locality, date_string, sep="_")]

  tmp <- dest_v2[,.N, tmpId]
  tmp[N>1]

### make subsample Id
  ssid <- function(br, tr, er, lr) {
    samps <- paste(br, tr, er, lr, sep=";")
    as.numeric(as.factor(samps))
  }

  dest.ag <- dest_v2[,list(subsample=ssid(br=bio_rep, tr=tech_rep, er=exp_rep, lr=loc_rep), tmpId), list(pre)]
  setkey(dest.ag, pre, tmpId)
  setkey(dest_v2, pre, tmpId)

  dest_v2 <- merge(dest_v2, dest.ag)

### make final sampleId
  setnames(dest_v2, "sampleId", "sampleId_orig")
  dest_v2[,sampleId:=paste(locality, subsample, date_string, sep="_")]
    table(dest_v2[,.N,sampleId]$N)
    dest_v2[,.N,sampleId][N>1]

### get elevation
  library(elevatr)


### output
  dest_v2 <- dest_v2[,-c("pre", "tmpId", "date_string"), with=F]

  colOrder <- c("sampleId", "sampleId_orig", "locality",
                "lat", "long",
                "continent", "country", "city", "province",
                "min_day", "max_day", "min_month", "max_month", "year", "jday", "exactDate",
                "bio_rep", "tech_rep", "exp_rep", "loc_rep", "fruit_type", "subsample",
                "nFlies", "fly_type", "library_type", "sampling_strategy", "seq_platform", "set",
                "collector", "SRA_Accession", "reference")

  dest_v2 <- dest_v2[,colOrder, with=F]
  setkey(dest_v2, sampleId)



  write.csv(dest_v2, quote=F, row.names=F, file="DESTv2/populationInfo/dest_v2.samps_12Jan2023.csv")





















### fix continents
  dest_v2[continent=="NorthAmerica", continent:="North_America"]
  dest_v2[lat<20 & continent=="Europe", continent:="North_America"]
  dest_v2[country=="Morocco", continent:="Africa"]
  dest_v2[lat< -30, continent:="Oceania"]

  table(dest_v2$continent)

### fix set
  table(dest_v2$set, dest_v2$continent)
  dest_v2[set=="extra"]

### fix type
  dest_v2[,type:=tolower(type)]
  table(dest_v2$type)

### fix city names
  dest_v2[city=="Chalet \xe0 Gobet", city:="Chalet Gobet"]
  dest_v2[city=="Purullena. El Bejar\xedn (Guadix)", city:="Purullena"]
  dest_v2[city=="Dond\xe9", city:="Donde"]

### get state/province
  # devtools::install_github("ropensci/rnaturalearth")
  # devtools::install_github("ropensci/rnaturalearthdata")
  # devtools::install_github("ropensci/rnaturalearthhires")

  library(rnaturalearth)
  library(sf)
  library(data.table)
  states <- ne_states(returnclass="sf")

  DT <- dest_v2[!is.na(lat) & !is.na(long),c("sampleId", "long", "lat"), with=F]
  setnames(DT, c("long", "lat"), c("longitude", "latitude"))
  DT_sf = st_as_sf(DT, coords = c("longitude", "latitude"),
                crs = 4326, agr = "constant")

  sf_use_s2(FALSE)
  dest_states <- cbind(DT, as.data.table(states)[st_nearest_feature(DT_sf, states)])
  dest_states[,c("sampleId", "iso_a2", "name", "name_alt", "postal"), with=F][is.na(postal)]



  dest_v2a <- merge(dest_v2, , by="sampleId", all.x=T)
  dest_v2a[,cc:=tstrsplit(sampleId, "_")[[1]]]
  table(dest_v2a[set=="DrosEU"]$cc==dest_v2a[set=="DrosEU"]$iso_a2)
  dest_v2a[cc!=iso_a2]

### output
  write.csv(dest_v2, "DESTv2/populationInfo/dest_v2_samples.csv", quote=F, row.names=F)

### basic plot
  sampPlot <- ggplot(dest_v2) +
              geom_line(aes(x=yday, y=lat, group=locality), linetype="dashed") +
              geom_point(aes(x=yday, y=lat, color=continent, shape=exactDate), size=.75) +
              facet_grid(set~year) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  sampPlot


### who to resequence?
  low_qual <- c(32, 50, 76, 77, 83, 87, 92, 96, 107, 120, 121, 127,136, 139, 144, 158, 185, 205, 210, 215, 233, 235, 236, 241, 243, 244, 246, 247, 249, 250, 252, 253)

  dest_v2 <- merge(dest_v2, droseu[,c("sampleId", "sampleId_internal"), with=F], all.x=T)
  dest_v2[sampleId_internal%in%low_qual, low_qual:=T]
  dest_v2[is.na(low_qual), low_qual:=F]

  tab <- dest_v2[,list(N_low_qual=sum(low_qual), N_total=.N), list(locality)]
  sum(tab[N_low_qual>0][N_low_qual!=N_total]$N_low_qual)





  ggplot(data=dest_v2) +  geom_point(aes(x=yday, y=lat), size=.75)



                 +
                facet_grid(~year) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




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
