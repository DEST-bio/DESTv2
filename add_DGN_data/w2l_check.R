 ### module load gcc/7.1.0  openmpi/3.1.4 R/3.6.1; R


### libraries
  library(data.table)
  library(foreach)

### loop
  o <- foreach(i=c("2L",  "2R", "3L", "3R", "X"))%do%{
    fl <- system(paste("ls /scratch/aob2x/dest/dgn/longData/*",
                        i,
                        "* | rev | cut -f1 -d'/' | rev | cut -f1 -d'_' | awk '{A[$1]++}END{for(i in A)print i,A[i]}' | sort -k2,2n",
                        sep=""), intern=T)

    data.table(pop=tstrsplit(fl, " ")[[1]],
               N=tstrsplit(fl, " ")[[2]],
               chr=paste("chr", i,  sep=""))

  }
  o <- rbindlist(o)
  o.ag <- o[,list(chr2L=(N[chr=="chr2L"]),
                  chr2R=(N[chr=="chr2R"]),
                  chr3L=(N[chr=="chr3L"]),
                  chr3R=(N[chr=="chr3R"]),
                  chrX=(N[chr=="chrX"])), list(pop)]
