ls ./trim/*_1_val_1.fq.gz | while read id1
do
id2=${id1%%_1_val_1.fq.gz}"_2_val_2.fq.gz"
name=${id1##*/}
neme=${name%%_1_val_1.fq.gz}

hisat2 -p 20 -t -x /lustre/user/taowlab/wanghc/study/index/hg38/genome/genome -1 $id1 -2 $id2 -S ./hisat/${neme}.sam #-p thread -t output time 
done
