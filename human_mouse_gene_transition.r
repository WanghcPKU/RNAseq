
library(biomaRt)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl", 
host = "https://dec2021.archive.ensembl.org/")
mouse <- useMart("ensembl",dataset="mmusculus_gene_ensembl", 
host = "https://dec2021.archive.ensembl.org/")

oxidative_trans <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
       values = oxidative$X,mart = human,
       attributesL = c("mgi_symbol","chromosome_name","start_position"),
       martL = mouse,uniqueRows = T)
