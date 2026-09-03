# ============================================================================
# Line B - DADA2 16S pipeline (all-ASCII version)
# Cohorts: Xiao PRJNA302078 (single-end 454) / Plummer PRJNA592384,
#          PRJNA398590 (paired-end MiSeq)
# Input : <DATA_DIR>/sra/<PRJNA>/<run>.fastq.gz
# Output: <OUT_DIR>/<PRJNA>/  (seqtab, taxa, phyloseq objects)
# NOTE  : All paths must be ASCII - R on this Windows box cannot open
#         files on Chinese paths. Data moved to <DATA_DIR> (junction at
#         original location keeps old paths working).
# Params: truncLen/maxEE from 01_peek_qc.py + plummer_trial.R calibration.
# ============================================================================
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({
  library(dada2); library(phyloseq); library(DECIPHER); library(ShortRead)
})

args <- commandArgs(trailingOnly = TRUE)
cohort   <- if (length(args) >= 1) args[1] else "PRJNA302078"
out_dir  <- file.path("<OUT_DIR>", cohort)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- cohort parameter table (calibrated 2026-08-29) ------------------------
params <- list(
  PRJNA302078 = list(paired = FALSE, truncLen = 450, maxEE = 2),        # 454 SE
  PRJNA592384 = list(paired = TRUE,  truncLen = c(240, 220), maxEE = c(2, 2)),  # overlap ~33bp, merge 99.9%
  PRJNA398590 = list(paired = FALSE, truncLen = 280, maxEE = 2),         # R2 too poor for merge -> SE R1
  PRJNA691964 = list(paired = TRUE,  truncLen = c(0, 0), maxEE = c(2, 2))  # Gustin V4 2x150: reads 145-150bp, truncLen=150 DROPS short reads (DADA2 discards < truncLen), so use 0=no truncation; truncQ=2 + maxEE only. Retention ~94%.
)
if (!cohort %in% names(params)) stop("undefined params: ", cohort)
P <- params[[cohort]]
paired <- P$paired

# ---- 1. locate fastq -------------------------------------------------------
fq_dir <- file.path("<DATA_DIR>/sra", cohort)
if (!dir.exists(fq_dir)) stop("dir not found: ", fq_dir)
fqs <- sort(list.files(fq_dir, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
if (length(fqs) == 0) fqs <- sort(list.files(fq_dir, pattern = "fastq\\.gz$", full.names = TRUE))
if (length(fqs) == 0) stop("no fastq.gz in dir")
sample_names <- sub("_1\\.fastq\\.gz$|_2\\.fastq\\.gz$|\\.fastq\\.gz$", "", basename(fqs))

if (paired) {
  fqs_R1 <- fqs
  fqs_R2 <- sort(list.files(fq_dir, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))
  if (length(fqs_R1) != length(fqs_R2)) stop("R1/R2 file count mismatch")
}

# ---- 2. filter & trim ------------------------------------------------------
filt_dir <- file.path(out_dir, "filtered"); dir.create(filt_dir, showWarnings = FALSE)
filt_R1 <- file.path(filt_dir, paste0(sample_names, "_F.fastq.gz"))
if (paired) filt_R2 <- file.path(filt_dir, paste0(sample_names, "_R.fastq.gz"))

if (paired) {
  out <- filterAndTrim(fqs_R1, filt_R1, fqs_R2, filt_R2,
                       truncLen = P$truncLen, maxN = 0, maxEE = P$maxEE,
                       truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
} else {
  out <- filterAndTrim(fqs, filt_R1, truncLen = P$truncLen, maxN = 0, maxEE = P$maxEE,
                       truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
}
print(out); saveRDS(out, file.path(out_dir, "filter_stats.rds"))

# ---- 3. learn errors & denoise ---------------------------------------------
errF <- learnErrors(filt_R1, multithread = TRUE)
if (paired) errR <- learnErrors(filt_R2, multithread = TRUE)
dada_F <- dada(filt_R1, err = errF, multithread = TRUE, pool = FALSE)
if (paired) dada_R <- dada(filt_R2, err = errR, multithread = TRUE, pool = FALSE)

# ---- 4. merge (paired) / sequence table (single) ----------------------------
if (paired) {
  mergers <- mergePairs(dada_F, filt_R1, dada_R, filt_R2, verbose = TRUE)
  seqtab  <- makeSequenceTable(mergers)
} else {
  seqtab <- makeSequenceTable(dada_F)
}
cat("ASV count: ", ncol(seqtab), "\n")
saveRDS(seqtab, file.path(out_dir, "seqtab.rds"))

# ---- 5. chimera removal -----------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
cat("ASV after chimera removal: ", ncol(seqtab.nochim), "\n")
saveRDS(seqtab.nochim, file.path(out_dir, "seqtab_nochim.rds"))

# ---- 6. taxonomy (SILVA v138.1, DADA2-compatible trainset) -----------------
# Raw SILVA NR99 trunc has 2-20 semicolon levels (Eukaryota) which breaks
# assignTaxonomy. build_silva_trainset.py filtered Bacteria/Archaea and
# standardized to 7 levels.
# tryRC=TRUE: our ASV seqs are reverse-complement of SILVA reference (verified
# 2026-08-29); without tryRC all taxa assignments come back NA.
taxa_file <- "<DATA_DIR>/silva/silva_138.1_train_set_bac_arch.fa.gz"
if (file.exists(taxa_file)) {
  taxa <- assignTaxonomy(seqtab.nochim, taxa_file, multithread = TRUE, tryRC = TRUE)
  saveRDS(taxa, file.path(out_dir, "taxa.rds"))
} else {
  warning("SILVA training set missing, skip taxonomy: ", taxa_file)
}

# ---- 7. phyloseq export -----------------------------------------------------
if (exists("taxa")) {
  ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(taxa))
  saveRDS(ps, file.path(out_dir, "phyloseq.rds"))
}
write.csv(seqtab.nochim, file.path(out_dir, "asv_table.csv"))
cat("=== DONE:", cohort, "===\n")
