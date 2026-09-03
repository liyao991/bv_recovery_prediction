# summary_cohorts.R - print per-cohort DADA2 summary table (all-ASCII)
# Reads: filter_stats.rds, seqtab.rds, seqtab_nochim.rds, taxa.rds
# Usage: Rscript summary_cohorts.R
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({ library(dada2); library(phyloseq) })

cohorts <- c("PRJNA302078", "PRJNA592384", "PRJNA398590")
out_root <- "<OUT_DIR>"
cat(sprintf("%-14s %-8s %-10s %-10s %-12s %-10s %-8s %-8s %-8s %-8s\n",
    "COHORT", "SAMPLES", "IN_READS", "FILT_READS", "FILT_PCT",
    "RAW_ASVs", "NOCHIM_ASVs", "KEEP_PCT", "PHYLUM", "GENUS"))

for (c in cohorts) {
  d <- file.path(out_root, c)
  if (!dir.exists(d)) { cat(sprintf("%-14s MISSING\n", c)); next }
  fs <- file.path(d, "filter_stats.rds")
  st <- file.path(d, "seqtab.rds")
  sn <- file.path(d, "seqtab_nochim.rds")
  tx <- file.path(d, "taxa.rds")
  if (!(file.exists(fs) && file.exists(st) && file.exists(sn))) {
    cat(sprintf("%-14s INCOMPLETE\n", c)); next
  }
  out <- readRDS(fs)
  seqtab <- readRDS(st)
  seqtab.nc <- readRDS(sn)
  in_reads <- sum(out[, 1])
  filt_reads <- sum(out[, 2])
  filt_pct <- 100 * filt_reads / in_reads
  raw_asv <- ncol(seqtab)
  nc_asv <- ncol(seqtab.nc)
  keep_pct <- 100 * nc_asv / raw_asv
  phylum <- genus <- NA
  if (file.exists(tx)) {
    taxa <- readRDS(tx)
    phylum <- round(100 * mean(!is.na(taxa[, "Phylum"])), 1)
    genus  <- round(100 * mean(!is.na(taxa[, "Genus"])), 1)
  }
  cat(sprintf("%-14s %-8d %-10.0f %-10.0f %-11.1f %-12d %-10d %-8.1f %-8s %-8s\n",
      c, nrow(seqtab.nc), in_reads, filt_reads, filt_pct,
      raw_asv, nc_asv, keep_pct, phylum, genus))
}
cat("\n== SUMMARY_DONE ==\n")
