### libraries
  library(data.table)
  library(SeqArray)
  library(binom)

### load meta-data file
  samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

### open GDS for common SNPs (PoolSNP)
  genofile <- seqOpen("~/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds", allow.duplicate=T)

### all of the samples there?
  table(samps$sampleId%in%seqGetData(genofile, "sample.id"))
  samps$sampleId[!samps$sampleId%in%seqGetData(genofile, "sample.id")]
  seqGetData(genofile, "sample.id")[!seqGetData(genofile, "sample.id")%in%samps$sampleId]

### common SNP.dt
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        nAlleles=seqGetData(genofile, "$num_allele"),
                        id=seqGetData(genofile, "variant.id"))

  snp.dt <- snp.dt[nAlleles==2]
  seqSetFilter(genofile, snp.dt$id)

  snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

### samples
  newSamps <- seqGetData(genofile, "sample.id")
  dim(samps)
  newSamps[grepl("AU_Que_Inn_-1_2014-02-15", newSamps)]

### function
  getData <- function(chr="2L", start=14617051, end=14617051, conf.int="exact") {
    # chr="2L"; start=14617051; end=14617051

    ### parameters
    ## chr, start and end are the positions of the SNP
    ## conf.int is a charactger vector that is passed to the "methods" parameter of the binomial binom.confint function in the binomial package

    ### filter to target
      snp.tmp <- data.table(chr=chr, pos=start:end)
      setkey(snp.tmp, chr, pos)
      setkey(snp.dt, chr, pos)
      seqSetFilter(genofile, variant.id=snp.dt[J(snp.tmp), nomatch=0]$id)

    ### get annotations
      message("Annotations")
      tmp <- seqGetData(genofile, "annotation/info/ANN")
      len1 <- tmp$length
      len2 <- tmp$data

      snp.dt1 <- data.table(len=rep(len1, times=len1),
                            ann=len2,
                            id=rep(snp.dt[J(snp.tmp), nomatch=0]$id, times=len1))

    # Extract data between the 2nd and third | symbol
      snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
      snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]

    # Collapse additional annotations to original SNP vector length
      snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste(gene, collapse=",")),
                            list(variant.id=id)]

      snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
      snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]

    ### get frequencies
      message("Allele Freqs")

      ad <- seqGetData(genofile, "annotation/format/AD")
      dp <- seqGetData(genofile, "annotation/format/DP")

      if(class(dp)[1]!="SeqVarDataList") {
        dp.orig <- dp
        dp <- list()
        dp$length<-1
        dp$data <- dp.orig

      }
      af <- data.table(ad=expand.grid(ad$data)[,1],
                       dp=expand.grid(dp$data)[,1],
                       sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                       variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))

    ### tack them together
      message("merge")
      afi <- merge(af, snp.dt1.an, by="variant.id")
      afi <- merge(afi, snp.dt, by.x="variant.id", by.y="id")

      afi[,af:=ad/dp]

    ### calculate effective read-depth
      afis <- merge(afi, samps[,c("sampleId", "nFlies")], by="sampleId")

      afis[chr=="X", nEff:=round((dp*nFlies   )/(dp+nFlies-1))]
      afis[chr!="X", nEff:=round((dp*2*nFlies )/(dp+2*nFlies-1))]
      afis[,af_nEff:=round(af*nEff)/nEff]

    ### calculate confidence intervals
      if(!is.null(conf.int)) {
        # conf.int <- "exact"

        tmp <- binom.confint(x=afis[!is.na(af_nEff)]$af_nEff*afis[!is.na(af_nEff)]$nEff, n=afis[!is.na(af_nEff)]$nEff, methods=conf.int)
        afis[!is.na(af_nEff),lci:=tmp$lower]
        afis[!is.na(af_nEff),uci:=tmp$upper]

      }
    ### return
      afis[,-c("n"), with=F]
  }

### test
  data <- getData(start=14617051, end=14617051, chr="2L", conf.int="exact")
  data[sampleId=="AU_Que_Inn_-1_2014-02-15"]
