# export_species_matrix.R - 种级聚合丰度矩阵（分辨率对照用）
# 策略：Species 有明确种名（非 uncultured/sp./unclassified）-> 特征 = "Genus Species"
#       否则退化聚合为 "Genus unclassified"（属级残余）
# 输出：<DATA_DIR>/feat/<CO>_species_rel.csv（行名 = 样本 ID，格式同属级矩阵）
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({
  library(phyloseq)
})

out_dir <- "<DATA_DIR>/feat"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
cohorts <- if (length(args) >= 1) args else c("PRJNA302078", "PRJNA592384", "PRJNA691964")

for (co in cohorts) {
  ps <- readRDS(file.path("<OUT_DIR>", co, "phyloseq.rds"))
  otu <- as.data.frame(as.matrix(otu_table(ps)))   # samples x ASVs
  tax <- as.data.frame(tax_table(ps))

  gen <- as.character(tax$Genus)
  gen[is.na(gen) | gen %in% c("", "NA")] <- "Unassigned"
  sp <- as.character(tax$Species)
  sp[is.na(sp) | sp %in% c("", "NA")] <- "unclassified"
  # 明确种名：非 uncultured / 非 "sp." / 非 unclassified
  is_sp <- !grepl("uncultured|unclassified| sp\\.$", sp, ignore.case = TRUE)
  feat <- ifelse(is_sp, paste(gen, sp, sep = " "), paste(gen, "unclassified", sep = " "))
  # 种名里的多余空格/括号清理
  feat <- gsub("\\s+", " ", feat)

  smat <- t(apply(otu, 1, function(row) tapply(row, feat, sum)))
  smat <- as.data.frame(smat)
  smat <- smat[, order(colSums(smat), decreasing = TRUE)]
  smat_rel <- smat / rowSums(smat)

  write.csv(smat_rel, file.path(out_dir, paste0(co, "_species_rel.csv")), row.names = TRUE)
  cat(co, ": samples", nrow(smat_rel), "species-features", ncol(smat_rel), "-> saved\n")
}
cat("DONE\n")
