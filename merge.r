folder_path <- "/lustre/user/taowlab/wanghc/work/lvwc/rnaseq6/kidney/featurecounts/" # 替换为你的文件夹路径,该文件夹下面存放有featurecounts的结果

txt_files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE) # 匹配读取.txt 文件
#批量读取
txt_contents <- lapply(txt_files, read.table, header = TRUE, sep = "\t")
#修改列名
extract_column_name <- function(full_name) {
  match <- sub(".*\\.bam\\.([^.]+\\.[^.]+)\\..*", "\\1", full_name)
  return(match)
} #featurecounts 的最后一列counts的名字 形如 ..bam.kidney.EP3.LGF12286_L2_sorted.bam 仅保留kidney.EP3

processed_contents <- lapply(txt_contents, function(df) {
  df <- df[, c(1, ncol(df))]
  rownames(df) <- df[, 1]
  colnames(df)[1] <- "geneid"
  colnames(df)[2] <- extract_column_name(colnames(df)[2])
  return(df)
})
all <- Reduce(function(x, y) merge(x, y, by = "geneid"), processed_contents)
head(all)
