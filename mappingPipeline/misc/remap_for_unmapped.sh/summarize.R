# ijob -A biol4559-aob2x -c20 -p standard --mem=40G
# module load gcc/11.4.0 openmpi/4.1.4 R/4.3.1 samtools/1.17; R

### libraries
  .libPaths(c("/project/berglandlab/Rlibs_4.3.1/")); .libPaths()

  library(data.table)
  library(foreach)
  library(foreach)
  library(doMC)
  registerDoMC(20)

###
  setwd("/scratch/aob2x/dest/bam/")
  files <- system("ls -d *.nonDros.bam", intern=T)

  nonDros <- foreach(samp.i=files, .errorhandling="remove")%dopar%{
    #samp.i=files[1]
    message(samp.i)

    if(length(system(paste("ls -d ", samp.i, ".bai", sep=""), intern=T))==0) {
      system(paste("samtools index ", samp.i, sep=""))
    }

    dat <- as.data.table(system(paste("samtools idxstats ", samp.i, sep=""), intern=T))
    dat[,samp:=samp.i]
    dat[,chr:=tstrsplit(V1, "\t")[[1]]]
    dat[,chrLen:=as.numeric(tstrsplit(V1, "\t")[[2]])]
    dat[,nReads:=as.numeric(tstrsplit(V1, "\t")[[3]])]

    dat
  }
  nonDros <- rbindlist(nonDros)


  nonDros[,sampleId:=tstrsplit(samp, "\\.")[[1]]]
  #nonDros <- simContam

  save(nonDros, file="~/nonDros_DESTv2.Rdata")
  setkey(nonDros, chr)

### merge and summarize
  nonDros.ag <- nonDros[,list(nReads=sum(nReads)), list(sampleId)]
  samps <- fread("/scratch/aob2x/DESTv2/populationInfo/dest_v2.samps_8Jun2023.csv")
  nonDros.ag <- merge(nonDros.ag, samps, by="sampleId", all.y=T)
  nonDros.ag[is.na(nReads), state:="noInfo"]
  nonDros.ag[nReads==0, state:="noReads"]
  nonDros.ag[nReads>0, state:="pass"]

  table(nonDros.ag$state, nonDros.ag$set)
  table(nonDros.ag[set=="cville"]$state, nonDros.ag[set=="cville"]$year)
  table(nonDros.ag[set=="dest_plus"]$state, nonDros.ag[set=="dest_plus"]$collector)
  table(nonDros.ag[set=="DrosEU_3"]$state, nonDros.ag[set=="DrosEU_3"]$locality)


  write.table(nonDros.ag[state!="pass"][set%in%c("cville","DrosEU","DrosEU_3","DrosEU_3_sa","DrosRTEC"), c("state", "sampleId"), with=F], file="/scratch/aob2x/dest/missingSamples.delim", quote=F, row.names=F)
  write.table(nonDros.ag[state!="pass"][set%in%c("cville","DrosEU","DrosEU_3","DrosEU_3_sa","DrosRTEC"), -c("nReads", "state"), with=F][!is.na(SRA_Accession)], file="/scratch/aob2x/dest/missingSamples.sra.delim", quote=F, row.names=F, sep=",")

  guide <- fread("/scratch/aob2x/dros3.all.mapping.guide.txt")
  setkey(guide, sampleId)
  guideRedo <- merge(guide, nonDros.ag[state!="pass"][set%in%c("cville","DrosEU","DrosEU_3","DrosEU_3_sa","DrosRTEC"), c("state", "sampleId"), with=F])
  write.table(guideRedo, "/scratch/aob2x/dest/missingSamples.DrosEU_3.delim", quote=F, row.names=F, sep="\t")

###
head /scratch/aob2x/dest/missingSamples.delim
