ls ./bam/*ted.bam | while read id
do
name=${id##*/}
neme=${name%%.bam}"_count.txt"

featureCounts -T 10 -p -t exon -g gene_name -a /lustre/user/taowlab/wanghc/work/rna-seq/ng/hg38/hg38.gtf -o ./featurecounts/$neme $id 
#-T thread -p pair-end -t type exon -a annotation

done
