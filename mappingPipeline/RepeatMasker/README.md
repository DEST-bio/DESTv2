## README

Repeatmasker needs to be installed, e.g. using homebrew

```bash
brew install repeatmasker
```

First create and change into directory

```bash
mkdir /Users/mkapun/Documents/GitHub/DEST/RepeatMasker/ref
cd /Users/mkapun/Documents/GitHub/DEST/RepeatMasker/ref
```

get transposon libary and ref genome
```bash
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster//dmel_r6.12_FB2016_04/fasta/dmel-all-transposon-r6.12.fasta.gz
curl -O ftp://ftp.flybase.net/genomes/Drosophila_melanogaster//dmel_r6.12_FB2016_04/fasta/dmel-all-chromosome-r6.12.fasta.gz
```

only keep contig name in headers (no spaces)
```bash
python /Users/mkapun/Documents/GitHub/DEST/RepeatMasker/scripts/adjust-id.py \
/Users/mkapun/Documents/GitHub/DEST/RepeatMasker/ref/dmel-all-transposon-r6.12.fasta > /Users/mkapun/Documents/GitHub/DEST/RepeatMasker/ref/dmel-all-transposon-r6.12_fixed-id.fasta
```

repeat mask reference genome
```bash
RepeatMasker \
-pa 20 \
--lib /Users/mkapun/Documents/GitHub/DEST/RepeatMasker/ref/dmel-all-transposon-r6.12_fixed-id.fasta \
--gff \
--q \
/Users/mkapun/Documents/GitHub/DEST/RepeatMasker/ref/dmel-all-chromosome-r6.12.fasta.gz
```
