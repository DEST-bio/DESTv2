# DEST mapping tutorial

The following script will guide the user through a example of the mapping pipeline described in our paper using a toy dataset.

### Step 1. Define your working directory and  Download the DEST pipeline
Throughout our example we will use ${wd} as the base directory.
```{sh}
wd=.
git clone https://github.com/DEST-bio/DEST_freeze1.git
```

### Step 2. Check that you have the toy dataset and metadata
The toy dataset is provided in this github. Make sure you find the reads in the git folder.
Start by declaring this function:

```{sh}
Test_toy_Files () {
if [[ -e "$1" ]]
then
echo "Good News: The file ${1} exists."
else
echo "Bad News: Somthing is wrong! Did you downloaded the DEST repo?"
fi
}
```
Ok, now run:
```{sh}
Test_toy_Files ${wd}/DEST_freeze1/mappingPipeline/Tutorial/ToyReads_1.fastq.gz
Test_toy_Files ${wd}/DEST_freeze1/mappingPipeline/Tutorial/ToyReads_2.fastq.gz
Test_toy_Files ${wd}/DEST_freeze1/mappingPipeline/Tutorial/ToyExample_samps.csv
```
**If everything looks good proceed to step 3**. Keep on reading if you would like to learn more about the toy dataset.

#### Extra notes on the toy dataset:
These toy reads were generated using  [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/). In our case, we generated a toy dataset of 20k reads (10k paired end reads to be exact) by sampling a real FASTQ file.
```{sh}
#### NOT PART OF THE TUTORIAL DO NOT RUN #####
#Load bbmap --
module load gcc/9.2.0
module load bbmap

#Replace these with your own D. melanogaste reads!
forward=real_reads_1.fastq.gz
reverse=real_reads_2.fastq.gz

#Run the bbmap sampler
reformat.sh \
in1=$forward \
in2=$reverse \
out1=./ToyReads_1.fastq.gz \
out2=./ToyReads_2.fastq.gz \
samplereadstarget=10000 \
-Xmx 30G \
usejni=t

gzip ToyReads_1.fastq
gzip ToyReads_2.fastq
#### NOT PART OF THE TUTORIAL DO NOT RUN #####
```

### Step 3. Create a folder for the slurm output
Now  create a dump folder for all outputs. This is most relevant if running dest on a cluster.
```{sh}
mkdir ${wd}/slurmOutput
```

### Step 4. Create the singularity image
Now lets build the docker image. **Skip this step if you have already built the image!**
In this example we are loading the program "singularity" using the option "module load". **This may vary in your cluster**.
Make sure you are loading singularity in your environment before proceeding.

*If this is your first time doing this, this step may take a while. But you only have to do this once.*

```{sh}
module load singularity # <- Remember you may have to load "singularity" differently in your cluster.
<<<<<<< HEAD
singularity pull docker://destbio/dest_freeze1:latest

=======
singularity pull docker://destbio/dest_freeze1:latest
>>>>>>> 2a7e04846d756706b5ae35ee45b8e1993c03220c
```

### Step 5. Personalize your pipeline options
Remember to update your SLURM header for the file [runDocker.sh](https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/scripts/runDocker.sh). The default header (shown below) is provided as an example, and will not work on your cluster.
If you are unfamiliar with the SLURM header. Read more [here](https://slurm.schedmd.com/documentation.html).
```{sh}
#### NOT PART OF THE TUTORIAL DO NOT RUN #####
#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 11
#SBATCH -N 1 # on one node
#SBATCH -t 72:00:00 #<= this may depend on your resources
#SBATCH --mem 90G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/RunDest.%A_%a.err # Standard error
#### NOT PART OF THE TUTORIAL DO NOT RUN #####
```
Depending on your cluster you may have to add additional options such as:
```{sh}
#SBATCH -p <your partition, if applicable>
#SBATCH --account <your account name, if applicable>
```

Remember to update your options in the file [runDocker.sh](https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/scripts/runDocker.sh). **If you are running this code for tutorial purposes, then the default options (shown below) is what you want!**

```{sh}
#### NOT PART OF THE TUTORIAL DO NOT RUN #####
# The default options
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
<<<<<<< HEAD
  --do-snape

=======
  --do_snape

>>>>>>> 2a7e04846d756706b5ae35ee45b8e1993c03220c
  #### NOT PART OF THE TUTORIAL DO NOT RUN #####
```

### Step 6. Run the mapping pipeline
This step will run the tutorial pipeline. Notice that we are using tutorial-specific files. these tutorial only files are:

 * [Toy reads](https://github.com/DEST-bio/DEST_freeze1/tree/main/mappingPipeline/Tutorial)
 * [Toy metadata](https://github.com/DEST-bio/DEST_freeze1/blob/main/mappingPipeline/Tutorial/ToyExample_samps.csv)

The mapping pipeline works by extracting information from the metadata file. for example, declare this function.
```{sh}
What_will_i_run () {

pop=$( cat $1  | sed '1d' | cut -f1,14 -d',' | grep -v "NA" | sed "1q;d" | cut -f1 -d',' )
srx=$( cat $1 | sed '1d' | cut -f1,14 -d',' | grep -v "NA" | sed "1q;d" | cut -f2 -d',' )
numFlies=$( cat $1  | sed '1d' | cut -f1,12 -d',' | grep -v "NA" | sed "1q;d" | cut -f2 -d',' )

echo "The sample name will be --> ${pop}"
echo "The read name is --> ${srx}"
echo "The number of flies pooled is --> ${numFlies}"

}
```

Now run:
```{sh}
What_will_i_run ${wd}/DEST_freeze1/mappingPipeline/Tutorial/ToyExample_samps.csv
```

**Finally**, we are ready to run the tutorial pipeline:

```{sh}
sbatch --array=1-$( sed '1d' ${wd}/DEST_freeze1/mappingPipeline/Tutorial/ToyExample_samps.csv | wc -l  ) \
${wd}/DEST_freeze1/mappingPipeline/scripts/runDocker.sh \
${wd} \
${wd}/DEST_freeze1/mappingPipeline/Tutorial/ \
${wd}/example_output_freeze1 \
${wd}/DEST_freeze1/mappingPipeline/Tutorial/ToyExample_samps.csv
```
<<<<<<< HEAD
If the toy example looks fine, then you are ready to run the pipeline on the whole dataset!
=======
If the toy example looks fine (i.e., no errors were produced in the pipeline), then you are ready to run the pipeline on the whole dataset!

A succesful run should produce various files, for example:

* FastQC reports for the reads
* BAM and MPILEUP file for the population being processed
* Reports on read mapping (e.g., duplicates)
* A sync.gz files (including masked sites) and a bed.gz file.
* Various sync files specifically for SNAPE (e.g., masked sites, monomorphic sites, etc.)
>>>>>>> 2a7e04846d756706b5ae35ee45b8e1993c03220c
