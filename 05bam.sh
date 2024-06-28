ls ./hisat/ | while read id
do
name=${id##*/}
neme=${name%%.sam}".bam"
nume=${name%%.sam}"_sorted.bam"
samtools view -@ 20 -S ./hisat/$id -b > ./bam2/$neme #-@ thread
samtools sort -@ 20 ./bam2/$neme -o ./bam2/$nume
samtools index ./bam2/$nume
done
