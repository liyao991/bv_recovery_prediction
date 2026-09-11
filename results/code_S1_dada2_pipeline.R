#!/usr/bin/env Rscript
# ============================================================
# Supplementary Code S1
# DADA2 amplicion processing pipeline for vaginal microbiome 16S data
# Cohorts: PRJNA302078 (Xiao 2016, 454 SE), PRJNA592384 (Plummer① 2023, PE300), PRJNA691964 (Gustin 2022, PE250)
# Depends: dada2 >= 1.30, ShortRead, ggplot2
# ============================================================

library(dada2); packageVersion("dada2")  # >= 1.30
library(ShortRead); library(ggplot2)

# ---- 1. Path conventions (ASCII for Windows) ----
PATH_FASTQ <- "D:/lineB_data/raw/<PRJNA_ID>"   # per-cohort fastq dir
PATH_FILTERED <- "D:/lineB_data/filtered/<PRJNA_ID>"
PATH_TAX <- "D:/lineB_data/tax/silva_v138.1_assignTaxonomy.RData"

# ---- 2. Cohort-specific parameters ----
# PRJNA302078 (Xiao 2016): 454 single-end, V3-V4
# PRJNA592384 (Plummer 2023): MiSeq PE300, V3-V4
# PRJNA691964 (Gustin 2022): MiSeq PE250, V3-V4

cohort_params <- list(
  PRJNA302078 = list(type = "SE", truncLen = 450, maxEE = 2),
  PRJNA592384 = list(type = "PE", truncLen = c(240, 220), maxEE = c(2, 2)),
  PRJNA691964 = list(type = "PE", truncLen = c(0, 0), maxEE = c(2, 2))  # full overlap, no truncation
)

# ---- 3. Quality profiles (Supplementary Fig. S1) ----
plotQualityProfile(sample(10, file.path(PATH_FASTQ, list.files(PATH_FASTQ)[1])))

# ---- 4. Filter + trim ----
fnFs <- sort(list.files(PATH_FASTQ, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(PATH_FASTQ, pattern = "_R2_001.fastq.gz", full.names = TRUE))

if (cohort_params[[PRJNA_ID]]$type == "SE") {
  filts <- filterAndTrim(fnFs, file.path(PATH_FILTERED, basename(fnFs)),
                         truncLen = cohort_params[[PRJNA_ID]]$truncLen,
                         maxEE = cohort_params[[PRJNA_ID]]$maxEE,
                         multithread = FALSE)  # Windows: avoid DADA2 multithread
} else {
  filts <- filterAndTrim(fnFs, file.path(PATH_FILTERED, basename(fnFs)),
                         fnRs, file.path(PATH_FILTERED, basename(fnRs)),
                         truncLen = cohort_params[[PRJNA_ID]]$truncLen,
                         maxEE = cohort_params[[PRJNA_ID]]$maxEE,
                         multithread = FALSE)
}

# ---- 5. Error learning ----
errF <- learnErrors(file.path(PATH_FILTERED, basename(fnFs)), multithread = FALSE)
errR <- learnErrors(file.path(PATH_FILTERED, basename(fnRs)), multithread = FALSE)

# ---- 6. Sample inference + merging (PE only) ----
dadaFs <- dada(file.path(PATH_FILTERED, basename(fnFs)), err = errF, multithread = FALSE)
if (cohort_params[[PRJNA_ID]]$type == "PE") {
  dadaRs <- dada(file.path(PATH_FILTERED, basename(fnRs)), err = errR, multithread = FALSE)
  mergers <- mergePairs(dadaFs, filts, dadaRs, filts, verbose = TRUE)
  seqtab <- makeSequenceTable(mergers)
} else {
  seqtab <- makeSequenceTable(dadaFs)
}

# ---- 7. Chimera removal + taxonomy assignment ----
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE)
taxa <- assignTaxonomy(seqtab.nochim, refFasta = "D:/lineB_data/tax/silva_v138.1_train_set.fa.gz",
                       multithread = FALSE, tryRC = TRUE)   # tryRC = TRUE is critical for SILVA NR99

# ---- 8. Output ----
saveRDS(seqtab.nochim, file = paste0("D:/lineB_data/asv/<PRJNA_ID>_seqtab.RDS"))
saveRDS(taxa,          file = paste0("D:/lineB_data/asv/<PRJNA_ID>_taxa.RDS"))
