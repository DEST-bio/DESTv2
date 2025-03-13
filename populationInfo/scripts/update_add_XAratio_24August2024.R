### libraries
  library(data.table)
  library(ggplot2)

### load in XA data
  xa_1 <- readRDS("/Users/alanbergland/coverage.rds")
  xa_1$SampleID<-rownames(xa_1)
  xa_2 <- fread("/Users/alanbergland/Xwide_Autosomes_PiEstimates_SexRatio")

  m <- as.data.table(merge(xa_1, xa_2, by="SampleID", all=T))
  m[,sr:=2*(1-(covX/covA))]
  ggplot(data=m, aes(x=ratioAX, y=SexRatio, color=set)) + geom_point() + facet_grid(~set) + geom_hline(yintercept=1)
  ggplot(data=m, aes(x=ratioAX, y=sr, color=set)) + geom_point() + facet_grid(~set) + geom_hline(yintercept=1)





# ijob -A berglandlab -c20 -p standard --mem=100G
### module load gcc/11.4.0  openmpi/4.1.4 R/4.3.1; echo "R_LIBS_USER=~/R/goolf/4.3" > ~/.Renviron; R

### double check
  library(data.table)
  library(SeqArray)
  library(foreach)
  library(doMC)
  registerDoMC(20)

  samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv")

### open GDS
  genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.24Aug2024.ann.gds", allow.duplicate=TRUE)

### make SNP table
  seqResetFilter(genofile)
  snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                       nAlleles=seqGetData(genofile, "$num_allele"),
                       id=seqGetData(genofile, "variant.id"))

  snp.dt <- snp.dt[nAlleles==2] ### subset to sites with only two alleles
  seqSetFilter(genofile, snp.dt$id)

### get coverage for test sample
  cov.mat <- foreach(samp.i=samps$sampleId)%dopar%{
    #samp.i <- "US_Vir_Cha_1_2018-07-05"
    message(samp.i)
    seqSetFilter(genofile, variant.id=snp.dt$id, sample.id=samp.i)

    dp.mat <- seqGetData(genofile, "annotation/format/DP"  )
    snp.dt[,dp:=dp.mat[1,]]

    tmp.ag <- snp.dt[,list(X.mean=mean(dp[chr=="X"], na.rm=T), X.median=median(dp[chr=="X"], na.rm=T),
                           A.mean=mean(dp[chr%in%c("2L", "2R", "3L", "3R")], na.rm=T), A.median=median(dp[chr%in%c("2L", "2R", "3L", "3R")], na.rm=T),
                           Y.mean=mean(dp[chr=="X"], na.rm=T), Y.median=median(dp[chr=="X"], na.rm=T),
                          sampleId=samp.i)]
    tmp.ag
  }
  cov.mat <- rbindlist(cov.mat)
  cov.mat
