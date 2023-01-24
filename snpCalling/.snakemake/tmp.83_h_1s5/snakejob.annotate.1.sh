#!/bin/sh
# properties = {"type": "single", "rule": "annotate", "local": false, "input": ["/scratch/aob2x/DESTv2_output/sub_bcf/dest.2L.all.PoolSNP.001.5.test.bcf"], "output": ["/scratch/aob2x/DESTv2_output/dest.all.PoolSNP.001.5.test.ann.vcf.gz"], "wildcards": {"popset": "all", "method": "PoolSNP", "maf": "001", "mac": "5", "version": "test"}, "params": {}, "log": [], "threads": 1, "resources": {"ntasks_per_node": 10, "memory_limit": 20, "time_limit": 1440, "log_dir": "/scratch/aob2x/DESTv2_output/logs"}, "jobid": 1, "cluster": {}}
 cd /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling && \
PATH='/apps/software/standard/core/snakemake/6.0.5/bin':$PATH /apps/software/standard/core/snakemake/6.0.5/bin/python3.8 \
-m snakemake /scratch/aob2x/DESTv2_output/dest.all.PoolSNP.001.5.test.ann.vcf.gz --snakefile /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/.snakemake/tmp.83_h_1s5 /scratch/aob2x/DESTv2_output/sub_bcf/dest.2L.all.PoolSNP.001.5.test.bcf --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules annotate --nocolor --notemp --no-hooks --nolock \
--mode 2  --default-resources "ntasks_per_node=1" "memory_limit=10" "time_limit=120"  && touch /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/.snakemake/tmp.83_h_1s5/1.jobfinished || (touch /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/.snakemake/tmp.83_h_1s5/1.jobfailed; exit 1)

