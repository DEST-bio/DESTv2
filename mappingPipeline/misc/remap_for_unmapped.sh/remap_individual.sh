
## SLURM_ARRAY_TASK_ID=200

### combine keep lists:
### four columns: "sranum", "sampleId", "type" (iso or wild), "set" (2016 or OW)

  # cat /standard/vol186/bergland-lab/alyssa/2016filters/keepList_v4.txt | awk '{print $1","$1",CartersMountainIsofemale,2016"}' > /scratch/aob2x/dest/individual_data/samples2use.csv
  # grep "CM" /project/berglandlab/Dmel_Single_Individuals/Overwintering_2018_2019/Individuals_Metadata_OW_n5_n7_corrected.csv | \
  # cut -f20,28,22 -d',' | sed 's/ //g' | sed 's/\"//g' | awk -F',' '{print $1","$3","$2",OW"}' >> /scratch/aob2x/dest/individual_data/samples2use.csv

### 2016 flies: 4 sets of PE files per sample
  sranum=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dest/individual_data/samples2use.csv | cut -f1 -d',')
  sample=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dest/individual_data/samples2use.csv | cut -f2 -d',')
  set=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/dest/individual_data/samples2use.csv | cut -f4 -d',')

  echo $sample " / " $sranum " / " $set

### standardize input files
  if [[ "$set" == "2016" ]]; then
    echo "2016 sample"
    # ls -lh /standard/vol186/bergland-lab/alyssa/2016fastq/*_1_${sranum}*.gz
    # ls -lh /standard/vol186/bergland-lab/alyssa/2016fastq/*_2_${sranum}*.gz

    cat /standard/vol186/bergland-lab/alyssa/2016fastq/*_1_${sample}*.gz > /scratch/aob2x/dest/individual_data/fastq/${sranum}_1.fastq.gz
    cat /standard/vol186/bergland-lab/alyssa/2016fastq/*_2_${sample}*.gz > /scratch/aob2x/dest/individual_data/fastq/${sranum}_2.fastq.gz

    r1_filename=/scratch/aob2x/dest/individual_data/fastq/${sranum}_1.fastq.gz
    r2_filename=/scratch/aob2x/dest/individual_data/fastq/${sranum}_2.fastq.gz

  else
    echo "OW sample"
    r1_filename=/project/berglandlab/Dmel_Single_Individuals/Overwintering_2018_2019/raw_data/usftp21.novogene.com/raw_data/${sranum}/*1.fq.gz
    r2_filename=/project/berglandlab/Dmel_Single_Individuals/Overwintering_2018_2019/raw_data/usftp21.novogene.com/raw_data/${sranum}/*2.fq.gz
  fi



  
