ls ./hisat/ | while read id
do
name=${id##*/}
neme=${name%%.sam}".bam"
nume=${name%%.sam}"_sorted.bam"
samtools view -@ 20 -S ./hisat/$id -b > ./bam/$neme #-@ thread
samtools sort -@ 20 ./bam/$neme -o ./bam/$nume
samtools index ./bam/$nume
done
