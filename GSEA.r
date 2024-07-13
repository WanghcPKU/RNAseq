library(clusterProfiler)
library(org.Mm.eg.db)
res <- read.csv("kidney_select.csv")
rownames(res) <- res$X
kideny_sorted <- res[order(-res$log2FoldChange),] %>% as.data.frame()

geneid <- bitr(rownames(kideny_sorted), fromType = "SYMBOL",
            toType = c("ENTREZID"),
            OrgDb = org.Mm.eg.db)

kideny_sorted$SYMBOL <- rownames(kideny_sorted)
head(geneid)

need_DEG <- merge(kideny_sorted, geneid, by='SYMBOL')
geneList <- need_DEG$log2FoldChange
names(geneList) <- need_DEG$ENTREZID
geneList <- sort(geneList, decreasing = T)


gsea.GO <- gseGO(geneList     = geneList,
               ont          = "BP",  # "BP"、"MF"和"CC"或"ALL"
               OrgDb        = org.Mm.eg.db,#人类org.Hs.eg.db 鼠org.Mm.eg.db
               keyType      = "ENTREZID",
               pvalueCutoff = 1)   #实际为padj阈值可调整

GO_kk <- DOSE::setReadable(gsea.GO, 
                           OrgDb=org.Mm.eg.db,
                           keyType='ENTREZID')
head(GO_kk)
