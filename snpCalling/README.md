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
 * `popSet`: Population to use: only option is `all`
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
module load gcc/9.2.0 openmpi/3.1.6 python/3.7.7 snakemake/6.0.5
cd /scratch/aob2x/DESTv2/snpCalling
snakemake --profile /scratch/aob2x/DESTv2/snpCalling/slurm -n
```

Then, if everything looks OK, run:
```bash
mkdir /scratch/aob2x/DESTv2_output/
mkdir /scratch/aob2x/DESTv2_output/logs/

module load gcc/9.2.0 openmpi/3.1.6 python/3.7.7 snakemake/6.0.5
cd /scratch/aob2x/DESTv2/snpCalling
snakemake --profile /scratch/aob2x/DESTv2/snpCalling/slurm
```

cd /scratch/aob2x/DESTv2_output/logs

ls -lh /scratch/aob2x/DESTv2_output
ls -lh /scratch/aob2x/DESTv2_output/sub_vcfs/
less -S /scratch/aob2x/DESTv2_output/dest.all.PoolSNP.001.5.test.ann.vcf
ls -lh /scratch/aob2x/DESTv2_output/sub_bcf/


rm /scratch/aob2x/DESTv2_output/snpEff*
rm /scratch/aob2x/DESTv2_output/dest*
rm /scratch/aob2x/DESTv2_output/sub_vcfs/*
rm /scratch/aob2x/DESTv2_output/sub_bcf/*

cat /scratch/aob2x/DESTv2_output/logs/runSNP_calling.46009067.err
cat /scratch/aob2x/DESTv2_output/logs/*.46009067.err

cd /scratch/aob2x/DESTv2_output/sub_vcfs/

### Output files
VCF, BCF, and GDS files are output to `<working_directory>/dest.<popSet>.<method>.<maf>.<mac>.<version>.ann.*`.
