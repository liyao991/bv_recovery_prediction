#!/usr/bin/env Rscript
# ============================================================
# Supplementary Code S1  --  ILLUSTRATIVE, NOT RE-EXECUTED
#
# DADA2 amplicon processing pipeline for the vaginal 16S data sets used
# in this study. Cohorts: PRJNA302078 (Xiao 2016, 454 SE),
# PRJNA592384 (Plummer 2023, PE300), PRJNA691964 (Gustin 2022, PE250).
#
# !!! READ THIS BEFORE USING THIS FILE !!!
# This script documents the processing logic reported in Methods. It has
# NOT been re-executed against the raw reads since the analyses were run,
# because the raw FASTQ files are not redistributed here (see Data
# availability). It is therefore supplied as an ILLUSTRATIVE script, not
# as a one-command reproduction of the ASV tables.
#
# Known limitations, stated openly:
#   * Cohort-specific truncation lengths and the taxonomy reference above
#     are the values recorded in Supplementary Table S1; the exact mixing
#     of `truncLen = c(0, 0)` for PRJNA691964 relies on full overlap and
#     should be re-validated per sequencing run.
#   * Sample-to-patient pairing, exclusion of partner-intervention cohorts
#     (StepUp-1/StepUp-2), and the derived feature matrices are handled by
#     the companion scripts in `scripts/`, not here.
#   * SILVA v138.1 NR99 numbering and the `tryRC = TRUE` setting matter for
#     this data set; changing the reference build will shift assignments.
#
# Usage:  LINEB_DATA_DIR=/path/to/lineB_data Rscript code_S1_dada2_pipeline.R PRJNA302078
# Depends: dada2 >= 1.30, ShortRead, ggplot2
# ============================================================

library(dada2); packageVersion("dada2")  # >= 1.30
library(ShortRead); library(ggplot2)

# ---- 1. Arguments and path conventions (ASCII paths for Windows) ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript code_S1_dada2_pipeline.R <PRJNA_ID>")
}
PRJNA_ID <- args[[1]]

DATA_DIR      <- Sys.getenv("LINEB_DATA_DIR", "<DATA_DIR>")   # set to your local copy
PATH_FASTQ    <- file.path(DATA_DIR, "raw",      PRJNA_ID)
PATH_FILTERED <- file.path(DATA_DIR, "filtered", PRJNA_ID)
PATH_ASV      <- file.path(DATA_DIR, "asv")
PATH_REF      <- file.path(DATA_DIR, "tax", "silva_v138.1_train_set.fa.gz")

dir.create(PATH_FILTERED, recursive = TRUE, showWarnings = FALSE)
dir.create(PATH_ASV,      recursive = TRUE, showWarnings = FALSE)

# ---- 2. Cohort-specific parameters ----
# PRJNA302078 (Xiao 2016):    454 single-end, V3-V4
# PRJNA592384 (Plummer 2023): MiSeq PE300,    V3-V4
# PRJNA691964 (Gustin 2022):  MiSeq PE250,    V3-V4
cohort_params <- list(
  PRJNA302078 = list(type = "SE", truncLen = 450,       maxEE = 2),
  PRJNA592384 = list(type = "PE", truncLen = c(240, 220), maxEE = c(2, 2)),
  PRJNA691964 = list(type = "PE", truncLen = c(0, 0),     maxEE = c(2, 2))  # full overlap, no truncation
)
if (is.null(cohort_params[[PRJNA_ID]])) {
  stop(sprintf("No parameter set defined for cohort '%s'", PRJNA_ID))
}
p <- cohort_params[[PRJNA_ID]]

# ---- 3. Quality profiles (Supplementary Fig. S1) ----
fnFs <- sort(list.files(PATH_FASTQ, pattern = "_R1_001.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(PATH_FASTQ, pattern = "_R2_001.fastq.gz$", full.names = TRUE))
if (length(fnFs) == 0L) stop("No R1 FASTQ files found in ", PATH_FASTQ)
if (p$type == "PE" && length(fnRs) != length(fnFs)) {
  stop("R1/R2 file counts differ; check the download before proceeding")
}

plotQualityProfile(fnFs[seq_len(min(10L, length(fnFs)))])
if (p$type == "PE") plotQualityProfile(fnRs[seq_len(min(10L, length(fnRs)))])

# ---- 4. Filter + trim ----
filtFs <- file.path(PATH_FILTERED, basename(fnFs))

if (p$type == "SE") {
  filts <- filterAndTrim(fnFs, filtFs,
                         truncLen = p$truncLen,
                         maxEE    = p$maxEE,
                         multithread = FALSE)  # Windows: avoid DADA2 multithread
} else {
  filtRs <- file.path(PATH_FILTERED, basename(fnRs))
  filts  <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                          truncLen = p$truncLen,
                          maxEE    = p$maxEE,
                          multithread = FALSE)
}

# ---- 5. Error learning (reverse reads only for paired-end data) ----
errF <- learnErrors(filtFs, multithread = FALSE, randomize = TRUE, verbose = FALSE)
if (p$type == "PE") {
  errR <- learnErrors(filtRs, multithread = FALSE, randomize = TRUE, verbose = FALSE)
}

# ---- 6. Sample inference + merging (PE only) ----
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
if (p$type == "PE") {
  dadaRs  <- dada(filtRs, err = errR, multithread = FALSE)
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = FALSE)
  seqtab  <- makeSequenceTable(mergers)
} else {
  seqtab  <- makeSequenceTable(dadaFs)
}

# ---- 7. Chimera removal + taxonomy assignment ----
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE)
taxa <- assignTaxonomy(seqtab.nochim, refFasta = PATH_REF,
                       multithread = FALSE, tryRC = TRUE)   # tryRC = TRUE is critical for SILVA NR99

# ---- 8. Output ----
saveRDS(seqtab.nochim, file = file.path(PATH_ASV, paste0(PRJNA_ID, "_seqtab.RDS")))
saveRDS(taxa,          file = file.path(PATH_ASV, paste0(PRJNA_ID, "_taxa.RDS")))

cat(sprintf("%s done: %d samples x %d ASVs after chimera removal\n",
            PRJNA_ID, nrow(seqtab.nochim), ncol(seqtab.nochim)))
