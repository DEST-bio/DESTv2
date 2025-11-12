### Generate coverage stats on X and autosomes from VCF file

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
  cov.mat <- merge(cov.mat, samps[,c("sampleId", "nFlies", "Recommendation")])
  cov.mat[,sr.median:=2*(1 - (X.median/A.median))]
  cov.mat[,sr.mean:=2*(1 - (X.mean/A.mean))]

  save(cov.mat, file="~/dest2_XAcoverage.Rdata")

### merge with sample metadata
  ### libraries
    library(data.table)
    library(ggplot2)
    library(patchwork)

    samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv")

    load("/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/OriginalMetadata/dest2_XAcoverage.Rdata")

    samps <- merge(samps, cov.mat[,c("sampleId", "sr")])

    write.csv(samps, quote=F, row.names=F, file="/Users/alanbergland/Documents/GitHub/DESTv2/populationInfo/dest_v2.samps_24Aug2024.xa.csv")


### Basic plot
  ggplot(data=samps[Recommendation=="Pass"][set!="dgn"], aes(x=set, y=sr, color=set, group=set)) + geom_boxplot() + ylab("Estimated male proportion") + ylim(-.29, 1.29)






samps[sampling_strategy=="Fly trap",list(nfruits=length(unique(na.omit(fruit_type_curated))), nSamps=.N,
            nSamps_withfruit=length(sampleId[!is.na(fruit_type_curated)])), list(locality, continent)][order(nSamps_withfruit)]

samps[locality%in%c("TR_Ank_Yes", "UA_Kie_Vys")][,c("sampleId", "fruit_type_curated")]

fruits <- samps[sampling_strategy=="Fly trap" & continent=="Europe", list(.N), fruit_type_curated]
fruits[,apple:=grepl("Apple", fruit_type_curated)]
fruits[,list(n=sum(N)), list(apple)]

table(samps$sampling_strategy)


samps[,nEff:=(Cov*2*nFlies)/(Cov + 2*nFlies - 1)]
cor.test(samps[Recommendation=="Pass"]$pNpS, samps[Recommendation=="Pass"]$nFlies)
cor.test(samps[Recommendation=="Pass"]$pNpS, samps[Recommendation=="Pass"]$Cov)
cor.test(samps[Recommendation=="Pass"]$pNpS, samps[Recommendation=="Pass"]$nEff)

nf <- ggplot(data=samps[Recommendation=="Pass"][!is.na(pNpS)], aes(y=pNpS, x=nFlies)) + geom_point() + ggtitle("Number of Flies")
nc <- ggplot(data=samps[Recommendation=="Pass"][!is.na(pNpS)], aes(y=pNpS, x=Cov)) + geom_point() + ggtitle("Nominal Coverage")
ec <- ggplot(data=samps[Recommendation=="Pass"][!is.na(pNpS)], aes(y=pNpS, x=nEff)) + geom_point() + ggtitle("Effective Coverage")

library(patchwork)
layout <- "
A
B
C"

nf + nc + ec + plot_layout(design=layout) + plot_annotation(tag_level="A")

summary(lm(pNpS~nFlies+Cov, samps[Recommendation=="Pass"]))


  ### load in data
    xa_2 <- fread("/Users/alanbergland/Xwide_Autosomes_PiEstimates_SexRatio")
    load("~/dest2_XAcoverage.Rdata")
    cov.mat[,sr.median:=2*(1 - (X.median/A.median))]
    m <- merge(xa_2, cov.mat, by.x="SampleID", by.y="sampleId", all=T)
    samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv")

    m <- merge(m, samps, by.x="SampleID", by.y="sampleId", all=T)
    corPlot <- ggplot(data=m[set!="dgn"][Recommendation.y=="Pass"][!is.na(SexRatio)], aes(x=sr.mean, y=SexRatio, color=set)) + geom_point() + facet_grid(~set) + geom_hline(yintercept=1) + geom_vline(xintercept=1) + geom_abline(aes(slope=1, intercept=0))

    t1 <- lm(SexRatio~sr.mean, data=m[set!="dgn"][Recommendation.y=="Pass"][!is.na(SexRatio)])
    m[set!="dgn" & Recommendation.y=="Pass" & !is.na(SexRatio),resid:=resid(t1)]

    residPlot <- ggplot(data=m[set!="dgn"][Recommendation.y=="Pass"][!is.na(SexRatio)], aes(x=resid, y=pcrdup, color=set)) + geom_point() + facet_grid(~set)

    corPlot / residPlot








  ### load in XA data
    xa_1 <- readRDS("/Users/alanbergland/coverage.rds")
    xa_1$SampleID<-rownames(xa_1)
    xa_2 <- fread("/Users/alanbergland/Xwide_Autosomes_PiEstimates_SexRatio")
    load("~/dest2_XAcoverage.Rdata")
    cov.mat[,sr.median:=2*(1 - (X.median/A.median))]

    m <- as.data.table(merge(xa_1, xa_2, by="SampleID", all=T))

    m[,sr:=2*(1-(covX/covA))]

    m <- merge(m, cov.mat, by.x="SampleID", by.y="sampleId", all=T)
    ggplot(data=m, aes(x=ratioAX, y=SexRatio, color=set)) + geom_point() + facet_grid(~set) + geom_hline(yintercept=1)
    ggplot(data=m, aes(x=ratioAX, y=sr, color=set)) + geom_point() + facet_grid(~set) + geom_hline(yintercept=1)

    ggplot(data=m[!is.na(set)][Recommendation.y=="Pass"], aes(x=sr.median, y=SexRatio, color=set)) + geom_point() + facet_grid(~set) + geom_hline(yintercept=1) + geom_abline(aes(slope=1, intercept=0)) + geom_vline(xintercept=1)
    ggplot(data=m[!is.na(set)][Recommendation.y=="Pass"], aes(x=sr.median, y=sr.mean, color=set)) + geom_point() + facet_grid(~set) + geom_hline(yintercept=1) + geom_abline(aes(slope=1, intercept=0))
    ggplot(data=m[!is.na(set)][Recommendation.y=="Pass"], aes(x=sr.mean, y=SexRatio, color=set)) + geom_point() + facet_grid(~set) + geom_hline(yintercept=1) + geom_vline(xintercept=1) + geom_abline(aes(slope=1, intercept=0))

    m[!is.na(set)][Recommendation.y=="Pass"][set=="DrosEU_3"][sr.mean<.8][SexRatio>1]
