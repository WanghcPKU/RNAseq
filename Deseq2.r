###处理 counts 矩阵
s01 <- read.table("./featurecounts/S0-1-LGE12080_L2_sorted_count.txt",sep="\t",header=T)
s02 <- read.table("./featurecounts/S0-2-LGE12081_L2_sorted_count.txt",sep="\t",header=T)
s51 <- read.table("./featurecounts/5-1-LGE12082_L2_sorted_count.txt",sep="\t",header=T)
s52 <- read.table("./featurecounts/5-2-LGE12083_L2_sorted_count.txt",sep="\t",header=T)

head(s01)

s01 <- s01[,c(1,ncol(s01))]
s02 <- s02[,c(1,ncol(s02))]
s51 <- s51[,c(1,ncol(s51))]
s52 <- s52[,c(1,ncol(s52))]

colnames(s01) <- c("gene_id","s01")
colnames(s02) <- c("gene_id","s02")
colnames(s51) <- c("gene_id","s51")
colnames(s52) <- c("gene_id","s52")

li <- list(s01,s02,s51,s52)

all <- Reduce(merge,li)
head(all)

write.csv(all,"raw_count.csv")

####FPKM转换
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}
####TPM转换
Counts2TPM <- function(counts, effLen){
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))}

###featurecounts 的结果里有gene length
len <- select(s01,Geneid,Length)

row.names(len) <- len$Geneid
row.names(all) <- all$gene_id
eff_len <- len[row.names(all),]

all <- all[,-1]
fpkm <- data.frame(apply(all,2,countToFpkm, effLen=eff_len$Length))
write.csv(fpkm,"fpkm.csv")

###correlation analysis

cc <- cor(all, method = "spearman")
head(cc)

library(pheatmap)

pdf("cor.pdf")
pheatmap(cc)
dev.off()

###DEG analysis
library(DESeq2)

condition <- factor(c(rep("s0",2),rep("s5",2)),
levels= c("s0","s5"))
colData <- data.frame(row.names=colnames(all), condition)
dds <- DESeqDataSetFromMatrix(all, colData, design= ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 3   #### 至少在3个样本中counts > 10 的gene
dds <- dds[keep, ]

vsdata <- vst(dds, blind=FALSE)  #归一化
rlogdata <- rlog(dds,blind = FALSE)
dds <- DESeq(dds)

res <- results(dds, contrast=c("condition","s5","s0"),independentFiltering = FALSE) ## independentFiltering 是否独立过滤 TRUE会在计算p值之前过滤掉噪声高的基因 但会有许多padj = 0 的情况

write.csv(res,"0_10.csv")

### GO term 
library(clusterProfiler)
library(org.Hs.eg.db)

res_sign <- subset(res, abs(log2FoldChange) > 0.5 & padj < 0.05)
write.csv(res_sign, "./res/res_sign.csv")


gene_up <- subset(res_sign, log2FoldChange > 0)$X

geneid <- bitr(gene_up, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
gene.GO <- enrichGO(gene = geneid$ENTREZID,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "ALL",
            pAdjustMethod ="BH",
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            readable = T)

gene.GO <- as.data.frame(gene.GO)
head(gene.GO)
write.csv(gene.GO, "./res/upgo.csv")

####可以使用simplify函数 或 simplifyEnrichment 简化GO结果
go_sim <- simplify(gene.GO,cutoff=0.7,by="p.adjust",select_fun=min) 

library(simplifyEnrichment) # 只能对 BP MF CC的一种进行简化 要先subset结果
go_id_bp <- gene.GO$ID

mat <- GO_similarity(go_id_bp, ont = "MF", db="org.Hs.eg.db")
df <- simplifyGO(mat, plot = T)



gene_down <- subset(res_sign, log2FoldChange < 0)$X

geneid <- bitr(gene_down, fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Hs.eg.db)
gene.GO <- enrichGO(gene = geneid$ENTREZID,
            OrgDb = org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = "ALL",
            pAdjustMethod ="BH",
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            readable = T)

gene.GO <- as.data.frame(gene.GO)
head(gene.GO)
write.csv(gene.GO, "./res/downgo.csv")

### KEGG
### 只演示down up同理
ekegg <- enrichKEGG(gene = geneid$ENTREZID,
organism='hsa',keyType = 'kegg',pAdjustMethod = 'BH',
pvalueCutoff = 1,qvalueCutoff = 1,use_internal_data=T) ### use_internal_data 使用本地库 联网可能不返回结果

###KEGG 只能返回 geneid geneid 与 symbol 手动转换

replace_ids_with_genes <- function(id_string, mapping_df) {
  # 分割字符串为 ID 向量
  id_vector <- unlist(strsplit(id_string, "/"))
  # 将 ID 转换为数字类型（排除空字符串）
  id_vector <- as.numeric(id_vector[id_vector != ""])
  # 查找每个 ID 的对应基因名称
  gene_vector <- mapping_df$SYMBOL[match(id_vector, mapping_df$ENTREZID)]
  # 重新组合成所需格式，去掉 NA 的部分
  result <- paste(gene_vector[!is.na(gene_vector)], collapse = "/")
  return(result)
}

df_result <- ekegg %>%
  mutate(Gene_Column = sapply(geneID, replace_ids_with_genes, mapping_df = geneid))

write.csv(df_result, "./res/kegg.csv")
