# ijob -c20 --mem=20G -p standard -A berglandlab

module load samtools

cd /project/berglandlab/DEST/dest_mapped/pipeline_output/AT_gr_12_fall


samtools view -h -@ 20 -b -f 4 -F 264  AT_gr_12_fall.original.bam  > /scratch/aob2x/tmps1.bam
samtools view -h -@ 20 -b -f 8 -F 260  AT_gr_12_fall.original.bam  > /scratch/aob2x/tmps2.bam
samtools view -h -@ 20 -b -f 12 -F 256 AT_gr_12_fall.original.bam  > /scratch/aob2x/tmps3.bam


ls -lh /scratch/aob2x/tmps1.bam
ls -lh /scratch/aob2x/tmps2.bam
ls -lh /scratch/aob2x/tmps3.bam

tar czvf ~/AT_gr_12_fall.unmapped.tar.gz /scratch/aob2x/tmps1.bam /scratch/aob2x/tmps2.bam /scratch/aob2x/tmps3.bam
