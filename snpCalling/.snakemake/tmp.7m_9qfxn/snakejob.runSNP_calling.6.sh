#!/bin/sh
# properties = {"type": "single", "rule": "runSNP_calling", "local": false, "input": ["/scratch/aob2x/DESTv2_output/jobs.csv", "/project/berglandlab/DEST/dest_mapped"], "output": ["/scratch/aob2x/DESTv2_output/sub_vcfs/2L_412525_550032.all.PoolSNP.001.5.test.vcf.gz"], "wildcards": {"jobid": "2L_412525_550032", "popset": "all", "method": "PoolSNP", "maf": "001", "mac": "5", "version": "test"}, "params": {}, "log": [], "threads": 1, "resources": {"ntasks_per_node": 1, "memory_limit": 10, "time_limit": 120, "log_dir": "/scratch/aob2x/DESTv2_output/logs"}, "jobid": 6, "cluster": {}}
 cd /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling && \
PATH='/apps/software/standard/core/snakemake/6.0.5/bin':$PATH /apps/software/standard/core/snakemake/6.0.5/bin/python3.8 \
-m snakemake /scratch/aob2x/DESTv2_output/sub_vcfs/2L_412525_550032.all.PoolSNP.001.5.test.vcf.gz --snakefile /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/Snakefile \
--force -j --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/.snakemake/tmp.7m_9qfxn /scratch/aob2x/DESTv2_output/jobs.csv /project/berglandlab/DEST/dest_mapped --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules runSNP_calling --nocolor --notemp --no-hooks --nolock \
--mode 2  --default-resources "ntasks_per_node=1" "memory_limit=10" "time_limit=120"  && touch /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/.snakemake/tmp.7m_9qfxn/6.jobfinished || (touch /gpfs/gpfs0/scratch/aob2x/DESTv2/snpCalling/.snakemake/tmp.7m_9qfxn/6.jobfailed; exit 1)

