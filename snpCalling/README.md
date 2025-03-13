## Running the Pipeline with Snakemake

### Description
This is a short walkthrough to generate annotated VCF files from previously generated masked SYNC files with the scripts at [DEST_freeze1/mappingPipeline](https://github.com/DEST-bio/DEST_freeze1/tree/main/mappingPipeline).

This pipeline has three basic steps:
  1. Slice each gSYNC file into non-overlapping regions and call SNPs. SNP-calling can be performed using PoolSNP or SNAPE. Note that the VCF files produced for SNAPE are processed through PoolSNP using a special flag. (`DEST_freeze1/snpCalling/scatter_gather_annotate/scatter.sh`)
  2. Combine sub-sections of the VCF (`DEST_freeze1/snpCalling/scatter_gather_annotate/gather.sh`)
  3. Annotate the VCF and convert to various formats (`DEST_freeze1/snpCalling/scatter_gather_annotate/annotate.sh`)

### Dependencies
 * SLURM based cluster
 * Snakemake (tested with version 6.1.1, pip installable `pip install snakemake==6.1.1`)
 * snpEff (tested with version 4.3t)
 * Modules (loaded on the cluster)
   * htslib, bcftools, parallel, intel/18.0, intelmpi/18.0, mvapich2/2.3.1, R/3.6.3, python/3.6.6, vcftools/0.1.1, gcc/7.1.0 , openmpi/3.1.4

### Setup
The config file `slurm/config.yaml` defines the cluster specific snakemake profile. It tells snakemake how to interact with SLURM to schedule jobs and correctly allocate resources. The `cluster` field defines the default command to submit a job and should be changed to fit your available allocations/partitions and any other resource limits or preferences.

The config file `workflow.yaml` holds other pipeline parameters which should be changed to fit your needs:
 * `script_directory`: Where the `snpCalling` scripts are located. Should be the path to this directory (`DEST_freeze1/snpCalling`) wherever you have cloned this repo.
 * `working_directory`: Directory where all the data will be processed and where output will be written
 * `pipeline_output_directory`: Directory holding the masked SYNC files from the main pipeline output. They should all have the common suffix `*masked.sync.gz`
 * `popSet`: Population to use: `all` (to include DGN) or `PoolSeq`
 * `method`: Method to use for variant calling: `PoolSNP` or `SNAPE`
 * `maf`: Minimum allele frequency (only used with `PoolSNP` method; use NA for `SNAPE`)
 * `mac`: Minmum allele count (only used with `PoolSNP` method; use NA for `SNAPE`)
 * `version`: Name for the run (e.g. the date).
 * `poolsnp_jobs`: umber of jobs to break the `run_poolsnp.sh` step into
 * `snpEff_path`: Where to find the .jar for snpEff

### Running

add to path

First, do a dry run with snakemake. This outputs the jobs which will be submitted, checks that everything snakemake needs for initialization is present, checks for syntax issues, etc. From `DEST_freeze1/snpCalling`, run

```bash
module load gcc/11.4.0  openmpi/4.1.4 python/3.11.4 snakemake/7.24.2
cd ~/CompEvoBio_modules/utils/snpCalling
snakemake -f --profile ~/CompEvoBio_modules/utils/snpCalling/slurm -n

snakemake --profile ~/CompEvoBio_modules/utils/snpCalling/slurm --unlock
snakemake -Sf
```

Then, if everything looks OK, run:


```bash
mkdir /scratch/aob2x/compBio_SNP_25Sept2023
mkdir /scratch/aob2x/compBio_SNP_25Sept2023/logs

module load gcc/11.4.0  openmpi/4.1.4 python/3.11.4 snakemake/7.24.2
cd ~/CompEvoBio_modules/utils/snpCalling

#cp jobs_genome.csv /scratch/aob2x/DESTv2_output/jobs.csv

snakemake --profile ~/CompEvoBio_modules/utils/snpCalling/slurm --unlock

sbatch ~/CompEvoBio_modules/utils/snpCalling/runSnakemake.sh
sacct -j 64522247
sacct -u aob2x
cd /scratch/aob2x/compBio_SNP_28Sept2024
cat /scratch/aob2x/compBio_SNP_28Sept2024/logs/runSnakemake.64522247
```













cd /scratch/aob2x/DESTv2_output_26April2023/

ls -lh /scratch/aob2x/DESTv2_output | head
ls -lS /scratch/aob2x/DESTv2_output_PoolSNP_50/sub_vcfs/*SNAPE* | wc -l
ls -lh /scratch/aob2x/DESTv2_output_PoolSNP_50/sub_bcf/*2L*

ls -lth /scratch/aob2x/DESTv2_output_PoolSNP_50/sub_vcfs/ | less

less -S /scratch/aob2x/DESTv2_output_PoolSNP_50/sub_vcfs/3R_30151921_30289600.all.PoolSNP.001.50.25Feb2023.vcf.gz

less -S /scratch/aob2x/DESTv2_output/jobs.csv
ls -lh /scratch/aob2x/DESTv2_output/logs/

find /scratch/aob2x/DESTv2_output_PoolSNP_50/sub_vcfs/ -type f -name "*.vcf.gz" -size -10k | \
sed 's/norep.vcf.gz/\*/g' | sed 's/vcf.gz/\*/g' | sed 's/^/rm /g'



zcat /scratch/aob2x/DESTv2_output_PoolSNP_50/sub_vcfs/3L_6170581_6307704.all.PoolSNP.001.50.25Feb2023.vcf.gz | less
ls -lh /scratch/aob2x/DESTv2_output_PoolSNP_50/sub_vcfs/3L_6170581_6307704.all.PoolSNP.001.50.25Feb2023.vcf.gz



 rm /scratch/aob2x/DESTv2_output_26April2023/snpEff*
 rm /scratch/aob2x/DESTv2_output_26April2023/dest*
 rm /scratch/aob2x/DESTv2_output_26April2023/sub_vcfs/*
 rm /scratch/aob2x/DESTv2_output_26April2023/sub_bcf/*
 rm /scratch/aob2x/DESTv2_output_26April2023/logs/*


rm /scratch/aob2x/DESTv2_output/jobs.csv

cat /scratch/aob2x/DESTv2_output_PoolSNP_50/logs/runSNP_calling.46464670.err
cat /scratch/aob2x/DESTv2_output/logs/*.47142169*err | less -S

47122813
cd /scratch/aob2x/DESTv2_output/sub_vcfs/

### Output files
VCF, BCF, and GDS files are output to `<working_directory>/dest.<popSet>.<method>.<maf>.<mac>.<version>.ann.*`.
