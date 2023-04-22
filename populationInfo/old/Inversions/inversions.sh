## first convert VCF to sync

## convert VCF to SYNC
gunzip -c dest.all.PoolSNP.001.50.10Nov2020.ann.vcf.gz \
| parallel \
--jobs 20 \
--pipe \
-k \
--cat python3 vcf2sync.py \
--input {} \
| gzip > dest.all.PoolSNP.001.50.10Nov2020.sync.gz

mkdir inversion

## get count at inversion specific marker SNPs
gunzip -c dest.all.PoolSNP.001.50.10Nov2020.sync.gz \
| parallel \
--pipe \
--jobs 20 \
-k \
--cat python3 overlap_in_SNPS.py \
--source inversion_markers_v6.txt_pos \
--target  {} \
> inversion_markers.sync

NAMES=$(gunzip -c dest.all.PoolSNP.001.50.10Nov2020.ann.vcf.gz | head -50 | awk '/^#C/' | cut -f10- | tr '\t' ',')

## calculate average frequencies for marker SNPs
python3 inversion-freqs.py \
inversion_markers_v6.txt \
inversion_markers.sync \
$NAMES \
> inversion.af
