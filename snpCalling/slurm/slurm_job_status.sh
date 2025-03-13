#!/bin/bash
# Check status of Slurm job. 
# Source: https://www.embl.org/groups/bioinformatics-rome/blog/2022/05/snakemake-profile-5-handling-memory-and-timeout-errors/.  
# More info: https://snakemake.readthedocs.io/en/v7.3.8/tutorial/additional_features.html#using-cluster-status

jobid="$1"

output=`sacct -j "$jobid" --format State --noheader | head -n 1 | awk '{print $1}'`

if [[ $output =~ ^(COMPLETED).* ]]
then
  echo success
elif [[ $output =~ ^(RUNNING|PENDING|COMPLETING|CONFIGURING|SUSPENDED).* ]]
then
  echo running
else
  echo failed
fi