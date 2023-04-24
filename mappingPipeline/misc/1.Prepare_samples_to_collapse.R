#### REMAP "COLLAPSED" SAMPLES .... 
#### 

library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(doMC)
registerDoMC(4)

### Iterator --->
args = commandArgs(trailingOnly=TRUE)
i=as.numeric(args[1])
#model=as.character(args[2]) # all_seas ; NoCore20_seas
#nPerm = as.numeric(args[3])

### obtain the metadata -->
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_25Feb2023.qc_merge.csv"

## Fread can load directly from the web
meta <- fread(meta_git)
setDT(meta)
meta %>% class
### extract samples to collapse 
meta %>%
  filter(!is.na(collapsedSamples)) ->
  meta.collapsed

####
samp.ith = meta.collapsed[i]
####
####
sample.name = samp.ith$sampleId

samps.sra = strsplit(samp.ith$SRA_Accession, ";")
samps.iths = length(samps.sra[[1]])

##### create folder for download
system(paste("mkdir", sample.name))

##### start the download process
foreach(k = 1:samps.iths)%dopar%{
  message(k)
  download.command <-
    paste("module load sratoolkit/2.10.5; fastq-dump --outdir", paste(sample.name, "/",samps.sra[[1]][k], sep = ""), "--split-e", samps.sra[[1]][k], sep = " ")
  system(download.command)
}


