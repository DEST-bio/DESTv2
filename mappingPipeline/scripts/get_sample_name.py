#!/usr/bin/env python

import sys

path_to_samps = "/scratch/jho5ze/DEST_docker/samps.csv"
with open(path_to_samps, encoding = "ISO-8859-1") as src:
    sra_name = sys.argv[1].strip().split("/")[-1].split("_")[0]
    for ix, line in enumerate(src.readlines()):
        line = line.split(",")
        if line[13] == sra_name:
            print(line[0])
            break
