file_path <- "/chenfeilab/Coconut_srv1/PNUTS/Analysis/merge_RNA_TT"
setwd(file_path)

file_list <- list.files()

for (file in file_list) {
  df_temp <- read.table(file)
  df_temp_1 <- df_temp[, "log2FoldChange", drop = FALSE]
  df_temp_1$gene <- rownames(df_temp_1)
  assign(gsub("_vs_.*", "", file), df_temp_1)
}


df_merged <- merge(merge(merge(dTAG.1h_PNUTS.rep1,
                               dTAG.3h_PNUTS.rep1,
                               by = "gene"),
                         dTAG.3h_PNUTS.fkbp.rep1,
                         by = "gene"),
                  dTAG.12h_PNUTS.fkbp.rep1,
                  by = "gene")

names(df_merged) <- c("gene", "TT_1hvsDMSO", "TT_3hvsDMSO", "RNA_3hvsDMSO", "RNA_24hvsDMSO")
rownames(df_merged) <- df_merged$gene

df_merged_1 <- df_merged[, -1] %>% arrange(desc(TT_1hvsDMSO))


pdf("FC_compare.pdf", height = 10, width = 5)
df_merged_2 <- df_merged_1
df_merged_2[df_merged_2 > 2] <- 2
df_merged_2[df_merged_2 < -2] <- -2
ComplexHeatmap::Heatmap(df_merged_2,
                         cluster_columns = FALSE,
                         cluster_rows = FALSE,
                         col = colorRampPalette(c('blue', 'white', 'red'))(200),
                         column_split = split_col,
                         show_row_names = FALSE)
dev.off()

