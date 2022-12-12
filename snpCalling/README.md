## Running the Pipeline with Snakemake

### Description
This is a short walkthrough to generate annotated VCF files from previously generated masked SYNC files with the scripts at [DEST_freeze1/mappingPipeline](https://github.com/DEST-bio/DEST_freeze1/tree/main/mappingPipeline).

This pipeline has three basic steps:
  1. Slice the each gSYNC file into non-overlapping regions and call SNPs. SNP-calling can be performed using PoolSNP or SNAPE. Note that the VCF files produced for SNAPE are processed through PoolSNP using a special flag. (`DEST_freeze1/snpCalling/scatter_gather_annotate/scatter.sh`)
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
 * `poolseq_sync_directory`: Directory holding the masked SYNC files from the main pipeline output. They should all have the common suffix `*masked.sync.gz`
 * `other_sync_directory`: Directory holding the other masked SYNC files. Again, should have the common suffix `*masked.sync.gz`. Currently, this directory is where the DGN data are.
 * `popSet`: Population to use: either `all` or `PoolSeq`. `all` uses both paths listed above to find masked SYNC files. `PoolSeq` only uses `poolseq_sync_directory`
 * `method`: Method to use for variant calling: `PoolSNP` or `SNAPE`
 * `maf`: Minimum allele frequency (only used with `PoolSNP` method; use NA for `SNAPE`)
 * `mac`: Minmum allele count (only used with `PoolSNP` method; use NA for `SNAPE`)
 * `version`: Name for the run (e.g. the date).
 * `poolsnp_jobs`: umber of jobs to break the `run_poolsnp.sh` step into
 * `snpEff_path`: Where to find the .jar for snpEff

### Running

First, do a dry run with snakemake. This outputs the jobs which will be submitted, checks that everything snakemake needs for initialization is present, checks for syntax issues, etc. From `DEST_freeze1/snpCalling`, run
```bash
module load gcc/9.2.0  openmpi/3.1.6 python/3.7.7
snakemake --profile slurm -n
```

Then, if everything looks OK, run:
```bash
snakemake --profile slurm
```

### Output files
VCF, BCF, and GDS files are output to `<working_directory>/dest.<popSet>.<method>.<maf>.<mac>.<version>.ann.*`.
