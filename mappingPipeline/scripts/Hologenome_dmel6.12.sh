#!/bin/bash

################################################
# Script to construct reference "hologenome"   #
# for D. melanogaster and associated microbes  #
# modified from Casey M. Bergman (Manchester)  #
# 	Martin Kapun 17/10/2016		       #
################################################


# download Drosophila melanogaster reference genome from FlyBase.org data (v.6.12)
# unzip file (the fasta-file contains all chromosomes)
# clean-up header a bit to improve readability

mkdir /Volumes/HD2/Hologenome/raw
cd /Volumes/HD2/Hologenome/raw

## D. melanogaster (only contains ONE copy of mitochondrion)
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.12_FB2016_04/fasta/dmel-all-chromosome-r6.12.fasta.gz
gunzip -c dmel-all-chromosome-r6.12.fasta.gz | sed 's/ type=.*//g'> D_melanogaster_r6.12.fasta

## S. cerevisiae
curl -O http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-2-1_20150113.tgz
tar -xvzf S288C_reference_genome_R64-2-1_20150113.tgz
cat S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa | sed 's/^>.*\[chromosome=/>Saccharomyces_cerevisiae_/g' | sed 's/]//g' > S_cerevisiae.fasta
rm -r S288C_reference_genome_R64-2-1_20150113


# download genomes for microbes known to be associated with D. melanogaster from NCBI
# if necessary, combine multi-contig draft assemblies into single fasta file
# rename fasta headers to include species name

# Wolbachia pipientis
curl -O ftp.ensemblgenomes.org/pub/release-32/bacteria/fasta/bacteria_0_collection/wolbachia_endosymbiont_of_drosophila_melanogaster/dna/Wolbachia_endosymbiont_of_drosophila_melanogaster.ASM802v1.dna.chromosome.Chromosome.fa.gz
gunzip -c Wolbachia_endosymbiont_of_drosophila_melanogaster.ASM802v1.dna.chromosome.Chromosome.fa.gz | sed '/^>/s/.*/>W_pipientis/g'  > W_pipientis.fasta

# Pseudomonas entomophila
curl -O ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria/fasta/bacteria_2_collection/pseudomonas_entomophila_l48/dna/Pseudomonas_entomophila_l48.ASM2610v1.dna.toplevel.fa.gz
gunzip -c Pseudomonas_entomophila_l48.ASM2610v1.dna.toplevel.fa.gz | sed '/^>/s/.*/>P_entomophila/g'  > P_entomophila.fasta

# Commensalibacter intestini
curl -O ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria/fasta/bacteria_14_collection/commensalibacter_intestini_a911/dna/Commensalibacter_intestini_a911.ASM23144v1.dna.toplevel.fa.gz
gunzip -c Commensalibacter_intestini_a911.ASM23144v1.dna.toplevel.fa.gz | sed 's/>* .*/_C_intestini/g'  > C_intestini.fasta

# Acetobacter pomorum
curl -O ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria/fasta/bacteria_5_collection/acetobacter_pomorum_dm001/dna/Acetobacter_pomorum_dm001.ASM19324v1.dna.toplevel.fa.gz
gunzip -c Acetobacter_pomorum_dm001.ASM19324v1.dna.toplevel.fa.gz | sed 's/>* .*/_A_pomorum/g'  > A_pomorum.fasta

# Gluconobacter morbifer
curl -O ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria/fasta/bacteria_24_collection/gluconobacter_morbifer_g707/dna/Gluconobacter_morbifer_g707.ASM23435v1.dna.nonchromosomal.fa.gz
gunzip -c Gluconobacter_morbifer_g707.ASM23435v1.dna.nonchromosomal.fa.gz | sed 's/>* .*/_C_morbifer/g'  > C_morbifer.fasta

# Providencia burhodogranariea
curl -O ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria/fasta/bacteria_3_collection/providencia_burhodogranariea_dsm_19968/dna/Providencia_burhodogranariea_dsm_19968.ASM31485v2.dna.toplevel.fa.gz
gunzip -c Providencia_burhodogranariea_dsm_19968.ASM31485v2.dna.toplevel.fa.gz | sed 's/>* .*/_P_burhodogranariea/g'  > P_burhodogranariea.fasta

# Providencia alcalifaciens
curl -O ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria/fasta/bacteria_18_collection/providencia_alcalifaciens_dmel2/dna/Providencia_alcalifaciens_dmel2.ASM31487v2.dna.toplevel.fa.gz
gunzip -c Providencia_alcalifaciens_dmel2.ASM31487v2.dna.toplevel.fa.gz | sed 's/>* .*/_P_alcalifaciens/g'  > P_alcalifaciens.fasta

# Providencia rettgeri
curl -O ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria/fasta/bacteria_23_collection/providencia_rettgeri_dmel1/dna/Providencia_rettgeri_dmel1.ASM31483v2.dna.toplevel.fa.gz
gunzip -c Providencia_rettgeri_dmel1.ASM31483v2.dna.toplevel.fa.gz | sed 's/>* .*/_P_rettgeri/g'  > P_rettgeri.fasta

# Enterococcus faecalis
curl -O ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria/fasta/bacteria_72_collection/enterococcus_faecalis/dna/Enterococcus_faecalis.ASM69626v1.dna.toplevel.fa.gz
gunzip -c Enterococcus_faecalis.ASM69626v1.dna.toplevel.fa.gz | sed 's/>* .*//g'  > E_faecalis.fasta


## remove raw gzipped files
rm *.gz

gzip *


##### create D. melanogaster "hologenome" from individual species fasta files
cat /Volumes/HD2/Hologenome/raw/D_melanogaster_r6.12.fasta /Volumes/HD2/Hologenome/raw/S_cerevisiae.fasta /Volumes/HD2/Hologenome/raw/A_pomorum.fasta /Volumes/HD2/Hologenome/raw/C_intestini.fasta /Volumes/HD2/Hologenome/raw/C_morbifer.fasta /Volumes/HD2/Hologenome/raw/E_faecalis.fasta /Volumes/HD2/Hologenome/raw/P_alcalifaciens.fasta /Volumes/HD2/Hologenome/raw/P_burhodogranariea.fasta /Volumes/HD2/Hologenome/raw/P_entomophila.fasta /Volumes/HD2/Hologenome/raw/P_rettgeri.fasta /Volumes/HD2/Hologenome/raw/W_pipientis.fasta | gzip > ../holo_dmel_6.12.fa.gz
