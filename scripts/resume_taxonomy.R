# resume_taxonomy.R - resume pipeline from seqtab_nochim (all-ASCII version)
# Steps: assignTaxonomy -> phyloseq -> write.csv
# Usage: Rscript resume_taxonomy.R <cohort>
# NOTE : multithread=TRUE is safe since SILVA trainset was deduplicated to
#        32,397 seqs (the earlier OOM was caused by the 450k-seq raw trainset,
#        NOT by multithreading - confirmed on PRJNA592384/PRJNA398590).
# NOTE : tryRC=TRUE is REQUIRED - our ASV seqs are reverse-complement of the
#        SILVA reference (verified 2026-08-29: fwd chain hits random genera,
#        rc chain hits Lactobacillus; without tryRC all taxa come back NA).
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({ library(dada2); library(phyloseq) })

args <- commandArgs(trailingOnly = TRUE)
cohort <- if (length(args) >= 1) args[1] else "PRJNA302078"
out_dir <- file.path("<OUT_DIR>", cohort)
st_file <- file.path(out_dir, "seqtab_nochim.rds")
if (!file.exists(st_file)) stop("seqtab_nochim.rds not found: ", st_file)
seqtab.nochim <- readRDS(st_file)
cat("Loaded seqtab_nochim:", nrow(seqtab.nochim), "samples x", ncol(seqtab.nochim), "ASVs\n")

taxa_file <- "<DATA_DIR>/silva/silva_138.1_train_set_bac_arch.fa.gz"
if (file.exists(taxa_file)) {
  cat("assignTaxonomy ...\n")
  taxa <- assignTaxonomy(seqtab.nochim, taxa_file, multithread = TRUE, tryRC = TRUE)
  saveRDS(taxa, file.path(out_dir, "taxa.rds"))
  cat("taxa saved:", nrow(taxa), "ASVs x", ncol(taxa), "levels\n")
} else {
  stop("SILVA file missing: ", taxa_file)
}

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(taxa))
saveRDS(ps, file.path(out_dir, "phyloseq.rds"))
cat("phyloseq saved, ntaxa =", ntaxa(ps), "\n")

write.csv(seqtab.nochim, file.path(out_dir, "asv_table.csv"))
cat("asv_table.csv written\n")

# quick summary of taxa assignment rates
ta <- as.data.frame(taxa)
cat("\nTaxa assignment completeness (% non-NA):\n")
for (lvl in colnames(ta)) {
  cat(sprintf("  %-10s %.1f%%\n", lvl, 100 * mean(!is.na(ta[[lvl]]))))
}
cat("\n== DONE:", cohort, "==\n")
