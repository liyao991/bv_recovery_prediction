# check_cst_feasibility.R - check species-level resolution for CST typing
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({
  library(dada2)
  library(phyloseq)
})

cohorts <- c("PRJNA302078", "PRJNA592384", "PRJNA398590")
for (co in cohorts) {
  cat("=====", co, "=====\n")
  ps <- readRDS(file.path("<OUT_DIR>", co, "phyloseq.rds"))
  otu <- as.data.frame(as.matrix(otu_table(ps)))
  tax <- as.data.frame(tax_table(ps))
  cat("n samples:", nrow(otu), " n ASVs:", ncol(otu), "\n")
  # genus-level abundance
  gen <- tax$Genus
  gen[is.na(gen)] <- "NA"
  abund <- colSums(otu)
  gsum <- tapply(abund, gen, sum)
  gsum <- sort(gsum, decreasing = TRUE)
  top <- head(gsum, 12)
  cat("top genera (% of total reads):\n")
  tot <- sum(abund)
  for (i in seq_along(top)) cat(sprintf("  %-20s %6.2f%%\n", names(top)[i], 100 * top[i] / tot))
  # Lactobacillus species-level resolution
  lb <- which(gen == "Lactobacillus")
  cat("Lactobacillus ASVs:", length(lb), "\n")
  sp <- tax$Species[lb]
  sp[is.na(sp)] <- "(unassigned)"
  spsum <- tapply(abund[lb], sp, sum)
  spsum <- sort(spsum, decreasing = TRUE)
  for (i in seq_along(head(spsum, 10))) cat(sprintf("  %-30s %6.2f%% of all reads\n", names(spsum)[i], 100 * spsum[i] / tot))
  # key BV genera present?
  bv <- c("Gardnerella", "Atopobium", "Prevotella", "Sneathia", "Megasphaera", "Mobiluncus", "Fannyhessea", "Dialister")
  cat("BV genera present:", paste(intersect(bv, names(gsum)), collapse = ", "), "\n")
  cat("\n")
}
