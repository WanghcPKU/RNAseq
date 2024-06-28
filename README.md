# RNASEQ DATA ANALYSIS
## STEP 0 
确认工作目录（Confirm your working directory）  
所在目录下面有raw文件夹储存原始文件（'raw' file store rawdata）  
## STEP 1 01fastq.sh  
质量控制（Quality Control）  
1.fastqc 主要查看GC% ; 接头 ; 碱基质量 （GC% ; Adaptor ; Base quality）  
2.multiqc  
## STEP 2 02trim.sh  
去接头（The raw RNAseq data is first trimmed by Trim Galore to remove adaptors）  
## STEP 3 03fastq.sh  
## STEP 4 04hisat.sh  
回帖基因组（The high quality RNAseq reads were mapped to the hg38 human / mm10 mouse genome by using Hisat2 ）  
nohup sh 04hisat.sh > mapped_ratio.out &  
index download : https://daehwankimlab.github.io/hisat2/download/  
## STEP 5 05bam.sh  
## STEP 6 06featurecounts.sh  
 'The number of reads mapped to each gene was counted using featureCounts software'  
 featureCounts software was used to map the number of reads to each gene.  
gtf download : https://www.gencodegenes.org/  
## STEP 7 Deseq2.r    
差异基因分析（Differentially expressed genes (DEGs) analysis）    
