### libraries
    library(data.table)
    library(foreach)
    library(doMC)
    registerDoMC(20)
    library(ggplot2)
    library(cowplot)
    library(epiR)

### load data
    dat <- fread("/mnt/pricey_4/sim_mel_contamination_simulation/summary/output.delim")

    dat[grepl("mel", V2),readOrigin:="mel"]
    dat[grepl("sim", V2),readOrigin:="sim"]

    dat[grepl("mel", V3),genome:="mel"]
    dat[grepl("sim", V3),genome:="sim"]

### some basic statistics
    dat.ag <- dat[,list(mm=sum(readOrigin=="mel" & genome=="mel"),
                        ss=sum(readOrigin=="sim" & genome=="sim"),
                        ms=sum(readOrigin=="mel" & genome=="sim"),
                        sm=sum(readOrigin=="sim" & genome=="mel")),
                   list(V1)]
   dat.ag[,estContamination := (ms+ss)/(ss+mm+sm+ms)]

   extractNum <- function(x) {
      # x <- "/mnt/pricey_4/sim_mel_contamination_simulation/mapped_reads/mix.simDown_0.melDown_1000000.bam"
        temp <- tstrsplit(rev(tstrsplit(x, "/"))[[1]][1], "\\.")
        simN <- as.numeric(gsub("simDown_", "", temp[[2]][1]))
        melN <- as.numeric(gsub("melDown_", "", temp[[3]][1]))
        simN/(simN+melN)
   }
   dat.ag[,trueContamination := sapply(dat.ag$V1, extractNum)]
   dat.ag[,estContamination := (sm+ss)/(sm+ss+ms+mm)]
   dat.ag[,fracSim_mappingTo_mel:=sm/(mm+ss+sm+ms)]

   epi.ccc(dat.ag$trueContamination, dat.ag$estContamination)


   est.plot <- ggplot(data=dat.ag[trueContamination>=1e-4], aes(x=qlogis(trueContamination), y=qlogis(estContamination))) +
   geom_point() + geom_abline(aes(slope=1, intercept=0)) +
   xlab("Simulated contamination level \n(logit scale)") + ylab("Estimated contamination level \n(logit scale)")


   resid.plot <- ggplot(data=dat.ag[trueContamination>=1e-4], aes(x=qlogis(trueContamination), y=fracSim_mappingTo_mel/estContamination)) + geom_point() +
   ylab("Residual cross mapping rate") + xlab("Simulated contamination level \n(logit scale)")

   combo.plot <- plot_grid(est.plot, resid.plot, labels=c("A", "B"))
   ggsave(combo.plot, file="~/contamPlot.png", heigh=4, w=8)


### sliding window

    setkey(dat, V3)
    stepSize <- 100000
    windowSize <- 100000
    wins <- foreach(chr.i=c("dmel_2L", "dmel_2R", "dmel_3L", "dmel_3R"), .combine="rbind")%do%{
            temp <- dat[J(chr.i)]

            data.table(V3=chr.i,
                       start=seq(from=0, to=max(temp$V4)-windowSize, by=stepSize),
                       stop=seq(from=0, to=max(temp$V4)-stepSize, by=stepSize) + windowSize)
    }

    setkey(dat, V3, V4)
    dat.win <- foreach(win.i=1:dim(wins)[1], .combine="rbind")%dopar%{
        print(paste(win.i, dim(wins)[1]), sep=" / ")

        dat.ag <- dat[J(data.table(V3=wins[win.i]$V3, V4=wins[win.i]$start:wins[win.i]$stop, by="V3,V4")),
                    list(mm=sum(readOrigin=="mel" & genome=="mel"),
                            ss=sum(readOrigin=="sim" & genome=="sim"),
                            ms=sum(readOrigin=="mel" & genome=="sim"),
                            sm=sum(readOrigin=="sim" & genome=="mel"),
                            chr=wins[win.i]$V3,
                            mid=(wins[win.i]$start + wins[win.i]$stop)/2,
                            start=wins[win.i]$start,
                            stop=wins[win.i]$stop),
                       list(V1),
                       nomatch=0]

       dat.ag[,estContamination := (ms+ss)/(ss+mm+sm+ms)]

       dat.ag
    }

    extractNum <- function(x) {
       # x <- "/mnt/pricey_4/sim_mel_contamination_simulation/mapped_reads/mix.simDown_0.melDown_1000000.bam"
         temp <- tstrsplit(rev(tstrsplit(x, "/"))[[1]][1], "\\.")
         simN <- as.numeric(gsub("simDown_", "", temp[[2]][1]))
         melN <- as.numeric(gsub("melDown_", "", temp[[3]][1]))
         simN/(simN+melN)
    }

    dat.win[,trueContamination := sapply(dat.win$V1, extractNum)]
    dat.win[,fracSim_mappingTo_mel:=sm/(mm+sm)]
    dat.win[,rr:=fracSim_mappingTo_mel/trueContamination]


    p <- c(3e-5, .002, 5e-4, .01, 2e-5, .003, .008, .006, .16, .007, .005, 6e-4, .009, .01, .02, .69, .26, .38, .06, .81, .14)
    pa <- p.adjust(p)


### is there a correlation between regions of the genome with high cross-species mappability and the number of seasonal SNPs?
    load("/mnt/pricey_1/dropPop/slidindWindow_switch_orig.100.5.Rdata")
    setkey(slidindWindow, chr.x, start, stop)
    o <- foreach(i=1:dim(dat.win)[1])%dopar%{
        print(paste(i, dim(dat.win)[1], sep=" / "))
        temp <- slidindWindow[J(data.table(chr.x=gsub("dmel_", "", dat.win$chr[i]))), nomatch=0][minPos>=dat.win$start[i]][maxPos>=dat.win$stop[i]]
        cbind(dat.win[i], data.table(frac.mu=mean(temp$frac)))
    }
    o <- rbindlist(o)

    o.ag <- o[,list(frac.mu=mean(frac.mu), fracSim_mappingTo_mel.mu=mean(fracSim_mappingTo_mel)), list(chr, mid)]


    ggplot(data=o, aes(x=fracSim_mappingTo_mel, y=frac.mu, color=chr)) + geom_point() +
    facet_wrap(~trueContamination)

    ggplot(data=dat.win, aes(x=mid, y=or, color=chr)) +
    geom_point() + facet_grid(trueContamination~chr)

















    dat.ag <- dat[,list(mm=sum(readOrigin=="mel" & genome=="mel"),
                        ss=sum(readOrigin=="sim" & genome=="sim"),
                        ms=sum(readOrigin=="mel" & genome=="sim"),
                        sm=sum(readOrigin=="sim" & genome=="mel")),
                   list(V1)]
   dat.ag[,estContamination := (ms+ss)/(ss+mm+sm+ms)]

   extractNum <- function(x) {
      # x <- "/mnt/pricey_4/sim_mel_contamination_simulation/mapped_reads/mix.simDown_0.melDown_1000000.bam"
        temp <- tstrsplit(rev(tstrsplit(x, "/"))[[1]][1], "\\.")
        simN <- as.numeric(gsub("simDown_", "", temp[[2]][1]))
        melN <- as.numeric(gsub("melDown_", "", temp[[3]][1]))
        simN/(simN+melN)
   }
   dat.ag[,trueContamination := sapply(dat.ag$V1, extractNum)]
   dat.ag[,estContamination := (sm+ss)/(sm+ss+ms+mm)]
   dat.ag[,fracSim_mappingTo_mel:=sm/(mm+ss+sm+ms)]

   plot(estContamination~trueContamination, dat.ag)
   abline(0,1)
