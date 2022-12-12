# Add DGN data

## Description
> This set of scripts downloads DGN data and formats the data into a gSYNC format in dm6. One file per chromosome per population


## -1. specify working directory
```bash
wd="/scratch/aob2x/dest"
```

## 0. Download all DGN data
> Needs a tab delimited file with jobID, prefix, path to DGN bz2 file: `DEST/add_DGN_data/dgn.list` <br/>
> Note that job 4 will fail. Why? Because 4 is the fourth line on DGN website for the DSPR. I don't think that we need to include that one.<br/>
```bash
sbatch --array=1-8 ${wd}/DEST/add_DGN_data/downloadDGN.sh
```
> OUT: ${wd}/dest/dgn/rawData<br/>


## 1. Unpack
> Each tarball is a bit different so the unpack scripts are different for each 1-8 (minus 4), from above. <br/>
 ```bash
sbatch --array=1-8 ${wd}/DEST/add_DGN_data/unpack.sh
```

## 2. Wide to long
> should be 4725 jobs <br/>
```bash
cd ${wd}/dgn/wideData/; ls * | tr '\t' '\n' | awk '{print >NR"\t"$0}' > ${wd}/dgn/dgn_wideFiles.delim
sbatch --array=1-$( tail -n1 ${wd}/dgn/dgn_wideFiles.delim | >cut -f1 ) ${wd}/DEST/add_DGN_data/wide2long.sh
```
> A quick check to make sure things look good:
> `w2l_check.R`

## 3. Format reference genome
```bash
 sbatch ${wd}/DEST/add_DGN_data/getRefGenome.sh
```

> Check that length of reference geomes are the same as the DGP long files: <br/>
> `wc -l ${wd}/referenceGenome/r5/2L.long` <br/>
> `wc -l ${wd}/dgn/longData/B_B04_Chr2L.seq.long` <br/>
> `wc -l ${wd}/referenceGenome/r5/3L.long` <br/>
> `wc -l ${wd}/dgn/longData/B_B04_Chr3L.seq.long` <br/>
> the reference geome is one line longer than the DGN data but that is due to a trailing empty line. No prob. <br/>

## 4. Generate gSYNC files for each population
> `pop_chr_maker.sh` generates the population/chromosome job id file

```bash
 ${wd}/DEST/add_DGN_data/pop_chr_maker.sh
 nJobs=$( tail -n1 ${wd}/dgn/pops.delim | cut -f1 )
 sbatch --array=1-${nJobs} ${wd}/DEST/add_DGN_data/makePopGenomeSync.sh
```

## 5. Liftover to dm6 and generate bgzipped gSYNC file
```bash
nJobs=$( tail -n1 ${wd}/dgn/pops.delim | cut -f1 )
sbatch --array=1-${nJobs} ${wd}/DEST/add_DGN_data/liftover_r5_to_r6.sh
```

## 6. Concatenate into one file per population, same as pooled; fix missing data fields
```bash
nJobs=$( cat ${wd}/dgn/pops.delim | cut -f3 | sort | uniq | awk '{print NR}'| tail -n1 )
sbatch --array=1-${nJobs} ${wd}/DEST/add_DGN_data/concatenate.sh
```

## 7. Move to output directory
```bash
nJobs=$( cat ${wd}/dgn/pops.delim | cut -f3 | sort | uniq | awk '{print NR}'| tail -n1 )
sbatch --array=1-${nJobs} ${wd}/DEST/add_DGN_data/move.sh
```
sacct -j 12449704
