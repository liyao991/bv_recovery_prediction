# plummer_se_trial.R - single-end (R1 only) trial for poor-quality cohorts
# Purpose: Plummer PRJNA398590 R2 quality too poor for reliable merge;
#          validate single-end R1 processing and pick truncLen.
# Usage: Rscript plummer_se_trial.R <cohort> <truncLen> [maxEE]
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
cohort    <- if (length(args) >= 1) args[1] else "PRJNA398590"
truncLen  <- if (length(args) >= 2) as.integer(args[2]) else 280
maxEE     <- if (length(args) >= 3) as.numeric(args[3]) else 2

fq_dir <- file.path("<DATA_DIR>/sra", cohort)
R1s <- sort(list.files(fq_dir, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
if (length(R1s) == 0) R1s <- sort(list.files(fq_dir, pattern = "fastq\\.gz$", full.names = TRUE))
n_trial <- min(3, length(R1s))
R1s <- R1s[1:n_trial]
cat("== SE trial:", cohort, "n =", n_trial, "truncLen =", truncLen, "maxEE =", maxEE, "==\n")

out_dir <- file.path("<OUT_DIR>", paste0(cohort, "_se_trial_", truncLen))
dir.create(out_dir, showWarnings = FALSE)
filt <- file.path(out_dir, paste0("S", 1:n_trial, "_F.fastq.gz"))

out <- filterAndTrim(R1s, filt, truncLen = truncLen, maxN = 0, maxEE = maxEE,
                     truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = 4)
print(out)
cat(sprintf("retention min/median: %.1f%% / %.1f%%\n",
            100 * min(out[, 2] / out[, 1]), 100 * median(out[, 2] / out[, 1])))

errF <- learnErrors(filt, multithread = 4)
dada_F <- dada(filt, err = errF, multithread = 4, verbose = FALSE)
seqtab <- makeSequenceTable(dada_F)
cat("ASV count:", ncol(seqtab), " length distribution:\n")
print(table(nchar(colnames(seqtab))))
saveRDS(seqtab, file.path(out_dir, "seqtab_se_trial.rds"))
cat("== SE TRIAL_DONE ==\n")
