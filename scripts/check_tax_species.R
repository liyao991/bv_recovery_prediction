# check_tax_species.R - 检查各队列 SILVA 种级注释比例
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages(library(phyloseq))
for (co in c("PRJNA302078", "PRJNA592384", "PRJNA691964")) {
  cat("===", co, "===\n")
  ps <- readRDS(file.path("<OUT_DIR>", co, "phyloseq.rds"))
  cat("  readRDS OK, class:", class(ps), "\n")
  tb <- tax_table(ps)
  cat("  tax_table dim:", nrow(tb), "x", ncol(tb), "\n")
  cat("  tax列:", paste(colnames(tb), collapse="/"), "\n")
  sp <- as.character(tb[, "Species"])
  sp[sp %in% c("", "NA", "uncultured")] <- NA
  cat("  Species非NA:", round(mean(!is.na(sp)), 3), "(", sum(!is.na(sp)), "/", length(sp), ")\n")
  ex <- head(sp[!is.na(sp)], 4)
  if (length(ex) > 0) cat("  示例:", paste(ex, collapse=" | "), "\n")
  cat("\n")
}
cat("DONE\n")
