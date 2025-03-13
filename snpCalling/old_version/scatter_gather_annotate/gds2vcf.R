#module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3; R

library(SeqArray)
library(Rsamtools)

#seqVCF2GDS("/scratch/aob2x/dest/dest.June14_2020.ann.vcf", "/scratch/aob2x/dest.June14_2020.ann.gds", storage.option="ZIP_RA")

args = commandArgs(trailingOnly=TRUE)
# args <- "dest.all.PoolSNP.001.50.26April2023.norep.ann.gds"
gds.fn=args[[1]]
vcf.fn=gsub(".gds", ".new.vcf.gz", gds.fn)

#vcf.fn=paste(vcf.fn, ".gz", sep="")
#vcf.fn="dest.all.PoolSNP.001.5.test.ann.vcf"
seqParallelSetup(cluster=20, verbose=TRUE)

genofile <- seqOpen(gds.fn)

seqGDS2VCF(genofile, vcf.fn, use_Rsamtools=T, verbose=T)
