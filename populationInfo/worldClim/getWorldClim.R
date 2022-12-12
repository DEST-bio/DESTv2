	library(raster)
  library(data.table)

# first load WC bio variables at the resolution of 2.5 deg
  biod <- getData("worldclim", var="bio", res=2.5)
  tmind <- getData("worldclim", var="tmin", res=2.5)
  tmaxd <- getData("worldclim", var="tmax", res=2.5)
  precd <- getData("worldclim", var="prec", res=2.5)
  ​
​
# extact for each coordinate bio clim variables
  bio<-extract(biod, samps[,c("long", "lat"), with=F])
  tmin<-extract(tmind, samps[,c("long", "lat"), with=F])
  tmax<-extract(tmaxd, samps[,c("long", "lat"), with=F])
  precd<-extract(precd, samps[,c("long", "lat"), with=F])
  ​
# create a full dataset
  bio.data <- as.data.table(cbind(samps[,c("sampleId"), with=F],bio,tmin,tmax,precd))
​
# save into external file
  write.csv(bio.data, file="./DEST_freeze1/populationInfo/worldClim/dest.worldclim.csv", row.names=FALSE)
​
​
