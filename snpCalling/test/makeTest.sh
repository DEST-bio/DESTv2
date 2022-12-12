module load tabix


tabix -b 2 -s 1 -e 2 /project/berglandlab/DEST/dest_mapped/RAL/RAL.masked.sync.gz 2L:1-1000000 |
bgzip -c > /scratch/aob2x/DEST_freeze1/snpCalling/test/dgn/RAL/RAL.masked.test.sync.gz
tabix -b 2 -s 1 -e 2 -f /scratch/aob2x/DEST_freeze1/snpCalling/test/dgn/RAL/RAL.masked.test.sync.gz


tabix -b 2 -s 1 -e 2 /project/berglandlab/DEST/dest_mapped/pipeline_output/WA_se_14_spring/WA_se_14_spring.SNAPE.monomorphic.masked.sync.gz 2L:1-1000000 |
bgzip -c > /scratch/aob2x/DEST_freeze1/snpCalling/test/PoolSeq/WA_se_14_spring/WA_se_14_spring.SNAPE.monomorphic.masked.test.sync.gz
tabix -b 2 -s 1 -e 2 -f /scratch/aob2x/DEST_freeze1/snpCalling/test/PoolSeq/WA_se_14_spring/WA_se_14_spring.SNAPE.monomorphic.masked.test.sync.gz


tabix -b 2 -s 1 -e 2 /project/berglandlab/DEST/dest_mapped/pipeline_output/WA_se_14_spring/WA_se_14_spring.masked.sync.gz 2L:1-1000000 |
bgzip -c > /scratch/aob2x/DEST_freeze1/snpCalling/test/PoolSeq/WA_se_14_spring/WA_se_14_spring.masked.test.sync.gz
tabix -b 2 -s 1 -e 2 -f /scratch/aob2x/DEST_freeze1/snpCalling/test/PoolSeq/WA_se_14_spring/WA_se_14_spring.masked.test.sync.gz
