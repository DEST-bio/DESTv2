### taken from `data-paper/Figure2/MakeFigure.r`


### libraries
library(tidyverse)

## at first read SNAPE pnps data and filter genomewide values
DATA.snape=read.table("./DEST_freeze1/populationInfo/sampleFilter/SNAPE_full.pnps",
                      header=T,
                      stringsAsFactors = F)
summary(DATA.snape)

DATA.snape.group <- DATA.snape %>%
  filter(Chrom =="genomewide") %>%
  group_by(MAF,Chrom) %>%
  summarise(
    n=n(),
    pnps.m = mean(pNpS),
    pnps.sd=sd(pNpS),
    pnps.se =sd(pNpS)/sqrt(sum(n()))
  )
#DATA.snape.group

## keep pnps without MAF filtering only
DATA.snape.group.MAF0 <- DATA.snape %>%
  filter(MAF == 0 & Chrom =="genomewide")
#DATA.snape.group.MAF0

# read CSV of private SNPs
DATA.ps=read.csv("./DEST_freeze1/populationInfo/sampleFilter/SNAPE_full.ps",
                 header=T,
                 stringsAsFactors = F,
                 sep = "\t")

DATA=merge(DATA.ps,DATA.snape.group.MAF0, by.x="POP",by.y="POP")

DATA$private=log10(DATA$N)

## calculated Mean/SD and threshold based on Mean+2SD for pNpS data
Mean.pNpS=mean(DATA$pNpS)
SD.pNpS=sd(DATA$pNpS)
th.pNpS=Mean.pNpS+1.96*SD.pNpS

##classify
DATA$TH.pNpS <-DATA$pNpS
DATA$TH.pNpS[DATA$TH.pNpS<th.pNpS]<-NA
DATA$TH.pNpS[DATA$TH.pNpS>=th.pNpS]<-"Exclude"
DATA$TH.pNpS[is.na(DATA$TH.pNpS)]<-"Keep"

## calculated Mean/SD and threshold based on Mean+1.96SD for private SNP data
Mean.private=mean(DATA$private)
SD.private=sd(DATA$private)
th.private=Mean.private+1.96*SD.private

#classify
DATA$TH.private <-DATA$private
DATA$TH.private[DATA$TH.private<th.private]<-NA
DATA$TH.private[DATA$TH.private>=th.private]<-"Exclude"
DATA$TH.private[is.na(DATA$TH.private)]<-"Keep"

# make new final classification, where any pop will be excluded that is excluded based on either measurement (pNpS and/or private SNPs)
DATA$Status <- rep("Keep",length(DATA$TH.private))
DATA$Status[DATA$TH.pNpS=="Exclude"]<-"Exclude"
DATA$Status[DATA$TH.private=="Exclude"]<-"Exclude"

write.table(DATA,"./DEST_freeze1/populationInfo/sampleFilter/classify_pops.txt",quote = F,row.names = F)
