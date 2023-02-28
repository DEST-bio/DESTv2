# module load  htslib/1.10.2 bcftools/1.9 intel/18.0 intelmpi/18.0 parallel/20200322 R/3.6.3; R

### libraries
  library(SeqArray)
  library(data.table)

### open GDS
  setwd("/project/berglandlab/DEST/gds")
  genofile <- seqOpen("dest.all.PoolSNP.001.50.25Feb2023.norep.ann.gds")
  length(seqGetData(genofile, "sample.id"))
  length(seqGetData(genofile, "variant.id"))


  genofile2 <- seqOpen("dest.all.PoolSNP.001.50.10Nov2020.ann.gds")
  length(seqGetData(genofile2, "sample.id"))
  length(seqGetData(genofile2, "variant.id"))

### open metadata
  setwd("/scratch/aob2x/DESTv2/")
  samps <- fread("populationInfo/dest_v2.samps_25Feb2023.csv")

### function
  w2l <- function(d, sampleId, variantId, var) {
    # d <- rd; sampleId=seqGetData(genofile, "sample.id"); variantId=seqGetData(genofile, "variant.id")
    dt <- data.table(x=expand.grid(d$data)$Var1, sampleId=rep(sampleId, length(variantId)), variantId=rep(variantId, each=length(sampleId)))
  }

### depth
  snps.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqGetData(genofile, "$num_allele"),
                        chr=seqGetData(genofile, "chromosome"))

  seqSetFilter(genofile, variant.id=sample(snps.dt[nAlleles==2]$variant.id, 1e5))

  rd <- seqGetData(genofile, "annotation/format/DP")

  rd.dt <-

  dt.ag <- dt[,list(nMissing=sum(is.na(x)), mrd=(mean(x, na.rm=T))), sampleId]

  dt.ag <- merge(dt.ag, samps, by="sampleId")
  dt.ag[,nEff:=(mrd*2*nFlies)/(mrd+2*nFlies)]
  summary(lm(mrd~set, dt.ag))

  dt.ag[,list(mu=median(nEff, na.rm=T)), list(set)]
  dt.ag[,list(mu=median(nEff, na.rm=T)), list(continent)]


### who is missing
  setkey(samps, sampleId)
  samps2 <- merge(samps, data.table(sampleId=seqGetData(genofile, "sample.id"), mapped=T, key="sampleId"), all.x=T, all.y=T)
  samps2[is.na(mapped), mapped:=F]

  table(samps2$mapped, samps2$set)
