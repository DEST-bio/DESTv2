library(adegenet)
library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(reshape2)
library(FactoMineR)
library(factoextra)
library(vcfR)
library(patchwork)
library(reshape2)
library(zoo)
library(matrixStats)
library(data.table)
library(SeqArray)
library(LEA)
library(gdsfmt)
library(SNPRelate)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggExtra)

source("./ThinLDinR_SNPtable.R")

# Part 1 -- Generate SNP file

### open GDS file
  genofile <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.ann.gds")

### get target populations
  samps <- fread("/scratch/yey2sn/DEST/samps.csv")

  samps <- rbind(samps[set=="DrosRTEC"],
                samps[set=="DrosEU"],
                samps[set=="dgn"]
                )
 
### get subsample of data to work on
  seqResetFilter(genofile)
  seqSetFilter(genofile, sample.id=samps$sampleId)

  snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile),
                        missing=seqMissing(genofile, .progress=T))

## choose number of alleles
 snps.dt <- snps.dt[nAlleles==2]

seqSetFilter(genofile, sample.id=samps$sampleId, variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

### select sites
  seqSetFilter(genofile, sample.id=samps$sampleId,
              snps.dt[chr%in%c("2L", "2R", "3L", "3R")][missing<.05][af>.2]$variant.id)

### get allele frequency data
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")

  dat <- ad$data/dp
  dim(dat)  
  
## Add metadata
    colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")
  
  rownames(dat) <- seqGetData(genofile, "sample.id")

# Remove unwated samples -- Because 

samples_to_remove = c(
                      "AT_gr_12_fall",
                      "AT_gr_12_spring",
                      "ES_ba_12_fall",
                      "ES_ba_12_spring",
                      "SIM",
                      "B",
                      "T"
                      )

dat_filt = dat[-which(rownames(dat) %in% samples_to_remove),]

left_join(data.frame(sampleId=rownames(dat_filt)), as.data.frame(samps)) -> DEST_DGN_metadata

# Remove Physical Linkage

snp_info = colnames(dat_filt) 

data.frame(t(dat_filt)) -> dat_filt_t

dat_filt_t %<>% mutate(SNPid=rownames(.)) %>% separate(SNPid, into = c("chr","pos","id"))

rownames(dat_filt_t) = snp_info

dat_filt_t$chr = as.character(dat_filt_t$chr) 
dat_filt_t$pos = as.numeric(dat_filt_t$pos)

picksnps_500<- pickSNPs(dat_filt_t,dist=500)

dat_filt_t_LD500 = dat_filt_t[picksnps_500,] 

dat_filt_t_LD500[,-which(names(dat_filt_t_LD500) %in% c("chr","pos","id") )] %>% t() -> dat_filt_LD500

#save(dat_filt_LD500,DEST_DGN_metadata, file="./DEST_DGN_AllSNPs_Metadata.Rdata")
#load("./DEST_DGN_AllSNPs_Metadata.Rdata")

#Convert NA cells into loci means
dat_filt_maf_LD500_naimp = na.aggregate(dat_filt_LD500)

dat_filt_maf_LD500_naimp %>% PCA(scale.unit = F, graph = F) -> WorldWide_PCA_object

WorldWide_PCA_object$ind$coord %>% as.data.frame() %>% mutate(sampleId = rownames(.)) %>% left_join(., DEST_DGN_metadata) -> PCA_coords_metadata

save(PCA_coords_metadata,WorldWide_PCA_object, file = "./Worlwide_PCA_objects.Rdata")

PCA_coords_metadata$country[which(PCA_coords_metadata$country == "United States")] = "USA"

# run cluster analysis
# 1. Loading and preparing data
df <- WorldWide_PCA_object$ind$coord

cluster_discovery = fviz_nbclust(df, kmeans, method = "gap_stat")
ggsave("cluster_discovery.pdf",cluster_discovery,  width =4, height = 4)

# 2. Compute k-means at K=4
set.seed(123)
km.res <- kmeans(df, 4, nstart = 25)

data.frame(Continental_clusters = km.res$cluster) %>% mutate(sampleId = rownames(.)) -> cluster_data

full_join(as.data.frame(samps[,c("sampleId","lat","long")]), cluster_data) -> Pop_Clusters

#Graph maps
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

ggplot(data = world) +
  geom_sf(fill= "antiquewhite") +
  coord_sf(xlim = c(-125.15, 45.00), ylim = c(-43.00, 65.00), expand = FALSE) + theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.5), panel.background = element_rect(fill = "aliceblue")) + geom_point(data =Pop_Clusters , aes(x=as.numeric(long), y=as.numeric(lat), fill = as.factor(Continental_clusters)),size = 3, shape = 21,  alpha = 0.7)   + xlab("Lon") + ylab("Lat") + ggtitle("Clustering") + theme(legend.position = "bottom") -> World_clusters

ggsave("Pop_clusters.pdf",World_clusters,  width =5, height = 4)

Pop_Clusters$Continental_clusters = gsub("1","1.Europe_W", Pop_Clusters$Continental_clusters)
Pop_Clusters$Continental_clusters = gsub("2","2.North_America", Pop_Clusters$Continental_clusters)
Pop_Clusters$Continental_clusters = gsub("3","3.Europe_E", Pop_Clusters$Continental_clusters)
Pop_Clusters$Continental_clusters = gsub("4","4.Africa", Pop_Clusters$Continental_clusters)

write.table( left_join(samps[,c("sampleId","country","city")], Pop_Clusters), 
             file = "./DEST_Sample_clusters.txt", 
             sep = "\t",quote = F ,row.names = F, col.names = T, append = F)
