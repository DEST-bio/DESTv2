# ijob -A berglandlab_standard -c1 -p dev --mem=4G

# module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
  library(data.table)
  library(SeqArray)

### load meta-data file
  samps <- fread("/scratch/aob2x/DESTv2/populationInfo/dest_v2.samps_13Jan2023.csv")

### open GDS for common SNPs (PoolSNP)
  genofile <- seqOpen("/scratch/aob2x/DESTv2_output/dest.all.PoolSNP.001.5.test.norep.ann.gds", allow.duplicate=T)

### common SNP.dt
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        nAlleles=seqGetData(genofile, "$num_allele"),
                        id=seqGetData(genofile, "variant.id"))

  snp.dt <- snp.dt[nAlleles==2]
  seqSetFilter(genofile, snp.dt$id)

  snp.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]


### function
  getData <- function(chr="2L", start=14617051, end=14617051) {
    # chr="2L"; start=14617051; end=14617051

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
      afis <- merge(afi, samps, by="sampleId")

      afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
      afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
      afis[,af_nEff:=round(af*nEff)/nEff]

    ### return
      afis[,-c("n"), with=F]
  }

### test
  data <- getData(start=13749337, end=13749337, chr="2L")
