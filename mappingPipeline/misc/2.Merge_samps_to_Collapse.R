### Merge reads from the same locality replicate, given thay pass filter.
### 

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

###

samps_to_join1 <- paste(sample.name, "/",samps.sra[[1]], "/", paste(samps.sra[[1]], "_1.fastq", sep = "") , sep = "")
samps_to_join2 <- paste(sample.name, "/",samps.sra[[1]], "/", paste(samps.sra[[1]], "_2.fastq", sep = "") , sep = "")

command_for_1 =
paste("cat ", paste(samps_to_join1, collapse = " ") , " > ", paste(sample.name, "/",sample.name,".joint_1.fastq", sep = "")  )

command_for_2 =
paste("cat ", paste(samps_to_join2, collapse = " ") , " > ", paste(sample.name, "/",sample.name,".joint_2.fastq", sep = "")  )

#### Run commands!
system(command_for_1)
system(command_for_2)

#### --> done
##
##
