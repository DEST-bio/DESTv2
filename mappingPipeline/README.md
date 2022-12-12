# Scripts to download, map, call polymorphism in pooled sequencing data-sets for Drosophila

## Description

This set of scripts provides a pipeline to build wholeGenomeSync files for each population sample from raw FASTQ data and defines a Dockerfile to build a docker image which can act as a standalone tool to run the pipeline.

Please read this repo to learn about the various options of our pipeline. We recommend running through our [TUTORIAL](https://github.com/DEST-bio/DEST_freeze1/tree/main/mappingPipeline/Tutorial) before using this script on the entire dataset. The tutorial uses a small toy dataset and only takes ~15-20 minutes to run (it may take a little longer if you have to build the docker image for the first time).

This repo is divided into various sections and users may start it at different point depending on their starting data. For example, users seeking to replicate our results from the DEST paper are advised to execute all steps. On the other hand, those using the data set on new data may start the pipeline at a different point depending if the data is to be downloaded from the SRA archive or if it exist locally in the user's cluster.

Please be advised that this repo assumes that the user will run the program on a cluster computer as it takes advantage of array jobs. Nevertheless, the script can be modified to run on a different configuration.

## Before we start: Download the DEST pipeline
### Define working directory
Prior to running this script, or any other in this pipeline for that matter, we advise the user to define a working directory link (e.g., stored in a variable). We assume that the git repo will be cloned here

```bash
#example
wd=./DEST_example
```

### Clone  git repo

```bash
cd ${wd}
git clone https://github.com/DEST-bio/DEST_freeze1.git
```

### Create a folder to store the SLURM output (e.g. ".out" and ".err")
This is optional. However, all of our scripts assume that all SLURM output will be dumped into a common folder created in the working folder, i.e., ${wd}.
```bash
mkdir ${wd}/slurmOutput
```
Now we can proceed to download some data

## Obtaining and preparing the data

### If the data is to be downloaded from the SRA
We have included in our pipeline a script designed to download SRA data from NCBI to the user's cluster.  This script can be found [here](https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/scripts/download_SRA.sh).</br>

#### Download data from SRA (specify 72 hour time limit)
This script uses the sratoolkit/2.10.5 to download desired samples from NCBI. For this script to run properly, the user will need to have installed the sratoolkit.  The user will also need a metadata file. This is looks like  [this](https://github.com/DEST-bio/DEST_freeze1/blob/main/populationInfo/samps.csv)</br> and our code depends on it to extract the information needed to complete the pipeline. In our particular case, our metadata file is  a csv file calles "samps.csv" and can be found in the repo in /DEST/populationInfo/

**about the script:** our script is launched using an array job of *1-N*, where *N* is the numbre of lines in the metadata file (thus we assume that there is one SRA file per metadata line). Because the metadata file has a header, we remove the first line. **Also, remember to modify the SLURM header of the script so it can be run in the user's cluster.** The script has requires **two** user inputs: First, the metadata file (e.g.,  ./downloadSRA.sh "metadata_file") the exact format shown [here](https://github.com/DEST-bio/DEST_freeze1/blob/main/populationInfo/samps.csv)</br> . Second, the location of the output file where the reads will be stored.

```bash
sbatch --array=1-$( sed '1d' ${wd}/DEST_freeze1/populationInfo/samps.csv | wc -l  ) \
${wd}/DEST_freeze1/mappingPipeline/scripts/download_SRA.sh  \
${wd}/DEST_freeze1/populationInfo/samps.csv \
${wd}/fastq
```
**Expected Output:** This script will produce fastq files in the form of SRR_1.fastq.gz. Where "_1" or "_2" are the forward and reverse reads of each run.

### If the data exist locally
If the data exist locally, i.e. pair-end reads, you must ensure that the file names have a structure that is identical to the expected outcome described above, e.g., FILEID_1.fastq.gz.
It is a good idea to rename your files to rename your files using this uniform convention
```bash
#For example
mv .../myfile.F.some_name.fastq.gz ./ID_1.fastq.gz
mv .../myfile.R.some_name.fastq.gz ./ID_2.fastq.gz
#Where "ID" is some identifier which makes sense to the user...
```
## Check that data are in the correct FASTQ format
Double check that all downloaded data are in Fastq33. Uses script from [here](https://github.com/brentp/bio-playground/blob/master/reads-utils/guess-encoding.py). </br> This script has 2 inputs: the location of the repo and the place where the reads are stored

```bash
sbatch ${wd}/DEST_freeze1/mappingPipeline/scripts/check_fastq_encoding.sh \
${wd}/DEST_freeze1 \
${wd}/fastq

#Finally run to show what files which have a different encoding.
grep -v "1.8" ${wd}/fastq/qualEncodings.delim
```

## Build singularity container from docker image
Download (and create the SIF image of) the docker image of the DEST mapping pipeline. This process may take a few minutes.
```bash
module load singularity # <- Remember you may have to load "singularity" differently in your cluster.
singularity pull docker://destbio/dest_freeze1:latest
```

## Pipeline options
These are the details of the array job script [runDocker.sh](https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/scripts/runDocker.sh) which executes the mapping pipeline. This script can run the pipeline in its entirety (from Fastq to final file), or partially. As such, it requires user attention on a couple of important options.

### User input for the script
The script takes 4 inputs:
1. Current folder address (where the SIF image is located)
2. The address to the folder containing the reads
3. The address to the output folder
4. The address to the metadata file

```bash
### This is an example
runDocker.sh \
<arg $1: working dir of the SIF file> \
<arg $2: reads> \
<arg $3: output folder> \
<arg $4: metadata>
```
### Options of the script
Our code has 2 parts. The first part "Get Sample information" leverages the array job task id "SLURM_ARRAY_TASK_ID" to iterate over the metadata file and extract all important information about the files. Accordingly, the script removes the header of the meatadata file and saves the population name "pop", SRX id "srx", and the number of flies pooled "numFlies". The script assumes that this information corresponds to the columns 1, 14, and 12, respectively in the metadata file (here stored in the $4 variable; i.e., argument 4) --> [make sure the file looks like this](https://github.com/DEST-bio/DEST_freeze1/blob/main/populationInfo/samps_10Nov2020.csv)</br>.


```bash
# This is an example. Do not Run
###################################
# Part  1. Get Sample information #
###################################

  pop=$( cat $4  | sed '1d' | cut -f1,14 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f1 -d',' )
  srx=$( cat $4 | sed '1d' | cut -f1,14 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )
  numFlies=$( cat $4  | sed '1d' | cut -f1,12 -d',' | grep -v "NA" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f2 -d',' )

```
The second part of the script executes the Docker image.  This script uses the reads in the format  "${srx}_1.fastq.gz".

The paramters for the docker pipeline are:

* **read1** Read name and path to PE1 fastq
* **read2** Read name and path to PE2 fastq
* **sample** Output sample name
* **output="."** Keep as is
* **threads** Number of threads for mapping. Needs to be consistent with SLURM request
* **max_cov=0.95** Upper threshold (quantile) for read depth filter
* **min_cov=10** Lower threshold (absolute number) for read depth filter
* **theta=0.005** Genetic diversity prior for SNAPE-pooled
* **D=0.01** Genetic differentiation prior between reference genome and population for SNAPE-pooled
* **priortype="informative"** Flag passed to SNAPE-pooled
* **fold="unfolded"** Flag passed to SNAPE-pooled
* **maxsnape=0.9** The posterior probability threshold for classifying site as polymorphic
* **nflies=40** Number of flies in sample
* **base_quality_threshold=25** Base quality filtering threshold
* **illumina_quality_coding=1.8** Illumin encoding
* **minIndel=5** Filtering distance to indels
* **do_prep** Runs steps to map raw reads and produces output BAM files. 1=run prep stage; 0=skip prep stage
* **do_poolsnp** Runs PoolSNP on BAM files. 1=run PoolSNP; 0=skip PoolSNP
* **do_snape** Runs SNAPE-pooled on BAM files. PoolSNP must be run before running snape. 1=run SNAPE-pooled; 0=skip SNAPE-pooled

In this example we are loading the program "singularity" using the option "module load". This may vary in your cluster. Make sure you are loading singularity in your environment before proceeding.

```bash
# This is an example. Do not Run
###################################
# Part  2. Run Docker             #
###################################

  module load singularity #<- Remember that loading singularity in your cluster may be different!

  singularity run \
  $1/dest_freeze1_latest.sif \
  $2/${srx}_1.fastq.gz \
  $2/${srx}_2.fastq.gz \
  ${pop} \
  $3 \
  --cores $SLURM_CPUS_PER_TASK \
  --max-cov 0.95 \
  --min-cov 4 \
  --base-quality-threshold 25 \
  --num-flies ${numFlies} \
  --do_poolsnp \
  --do_snape

```

### The SLURM header must be modified to your cluster or machine
Notice that our example scripts all use an example SLURM header that is unique to our supercomputer. **You will need to modify this for it to run in your machine!!**

```{sh}
#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 11
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00 #<= this may depend on your resources
#SBATCH --mem 90G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/RunDest.%A_%a.err # Standard error
```
Depending on your cluster you may have to add additional options such as:
```{sh}
#SBATCH -p <your partition, if applicable>
#SBATCH --account <your account name, if applicable>
```

## Running the singularity container across list of populations
Finally, the user can run the pipeline on their data using the following array job:
```bash
sbatch --array=1-$( sed '1d' ${wd}/DEST_freeze1/populationInfo/samps.csv | wc -l  ) \
${wd}/DEST_freeze1/mappingPipeline/scripts/runDocker.sh \ # The script
${wd} \ # Argument 1: Where the SIF file is located
${wd}/fastq \ # Argument 2: Where the reads are located
${wd}/pipeline_output \ # Argument 3: Output folder
${wd}/DEST_freeze1/populationInfo/samps_10Nov2020.csv
```

After this step has been completed, you are ready to move to the second step of the pipeline. "SNP calling, VCF and GDS generation"
