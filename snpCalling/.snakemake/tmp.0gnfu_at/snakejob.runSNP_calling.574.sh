#!/bin/sh
# properties = {"type": "single", "rule": "runSNP_calling", "local": false, "input": ["/scratch/aob2x/DESTv2_output/jobs.csv", "/project/berglandlab/DEST/dest_mapped"], "output": ["/scratch/aob2x/DESTv2_output/sub_vcfs/3R_1101441_1239120.all.PoolSNP.001.5.test.norep.vcf.gz"], "wildcards": {"jobid": "3R_1101441_1239120", "popset": "all", "method": "PoolSNP", "maf": "001", "mac": "5", "version": "test"}, "params": {}, "log": [], "threads": 1, "resources": {"ntasks_per_node": 4, "memory_limit": 36, "time_limit": 120, "log_dir": "/scratch/aob2x/DESTv2_output/logs"}, "jobid": 574, "cluster": {}}
 cd /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling && \
PATH='/apps/software/standard/core/snakemake/6.0.5/bin':$PATH /apps/software/standard/core/snakemake/6.0.5/bin/python3.8 \
-m snakemake /scratch/aob2x/DESTv2_output/sub_vcfs/3R_1101441_1239120.all.PoolSNP.001.5.test.norep.vcf.gz --snakefile /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/.snakemake/tmp.0gnfu_at /scratch/aob2x/DESTv2_output/jobs.csv /project/berglandlab/DEST/dest_mapped --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules runSNP_calling --nocolor --notemp --no-hooks --nolock \
--mode 2  --default-resources "ntasks_per_node=4" "memory_limit=36" "time_limit=120"  && touch /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/.snakemake/tmp.0gnfu_at/574.jobfinished || (touch /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/.snakemake/tmp.0gnfu_at/574.jobfailed; exit 1)

