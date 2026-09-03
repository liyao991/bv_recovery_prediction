# export_genus_matrix.R - export genus-level abundance + alpha diversity
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
})

out_dir <- "<DATA_DIR>/feat"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
cohorts <- if (length(args) >= 1) args else c("PRJNA302078", "PRJNA592384", "PRJNA691964")
for (co in cohorts) {
  ps <- readRDS(file.path("<OUT_DIR>", co, "phyloseq.rds"))
  otu <- as.data.frame(as.matrix(otu_table(ps)))   # samples x ASVs
  tax <- as.data.frame(tax_table(ps))
  
  # genus-level: aggregate ASV abundances by genus
  gen <- as.character(tax$Genus)
  gen[is.na(gen)] <- "Unassigned"
  # genus x sample matrix
  gmat <- t(apply(otu, 1, function(row) tapply(row, gen, sum)))
  gmat <- as.data.frame(gmat)
  gmat <- gmat[, order(colSums(gmat), decreasing = TRUE)]
  
  # relative abundance
  gmat_rel <- gmat / rowSums(gmat)
  
  # alpha diversity (on ASV level)
  shannon <- diversity(otu, index = "shannon")
  chao1 <- estimateR(round(otu))["S.chao1", ]
  n_asv_obs <- specnumber(otu)
  
  meta_df <- data.frame(
    sample = rownames(gmat),
    shannon = shannon,
    chao1 = chao1,
    n_asv = n_asv_obs,
    stringsAsFactors = FALSE
  )
  
  write.csv(gmat_rel, file.path(out_dir, paste0(co, "_genus_rel.csv")), row.names = TRUE)
  write.csv(meta_df, file.path(out_dir, paste0(co, "_alpha.csv")), row.names = FALSE)
  cat(co, ": samples", nrow(gmat_rel), "genera", ncol(gmat_rel), "-> saved\n")
}
