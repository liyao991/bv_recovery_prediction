# ============================================================================
# Plummer 双端 merge 可行性试跑（3 样本）
# 目的：验证 truncLen c(240,200)/c(250,200) 下 R1/R2 重叠是否足够、merge 率是否可接受
# 用法: Rscript plummer_trial.R <cohort>  (默认 PRJNA592384)
# 输出: stdout 打印 merge 统计 + <OUT_DIR>/<cohort>_trial/ 下中间产物
# ============================================================================
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({
  library(dada2); library(ShortRead)
})

args <- commandArgs(trailingOnly = TRUE)
cohort <- if (length(args) >= 1) args[1] else "PRJNA592384"
# CLI: Rscript plummer_trial.R <cohort> <truncR1> <truncR2> [maxEE1] [maxEE2] [maxMismatch]
if (length(args) >= 3) {
  truncLen <- as.integer(c(args[2], args[3]))
  maxEE    <- if (length(args) >= 5) as.numeric(c(args[4], args[5])) else c(2, 2)
} else {
  params <- list(
    PRJNA592384 = list(truncLen = c(240, 220), maxEE = c(2, 2)),
    PRJNA398590 = list(truncLen = c(250, 250), maxEE = c(2, 2))
  )
  P <- params[[cohort]]
  if (is.null(P)) stop("undefined params: ", cohort)
  truncLen <- P$truncLen; maxEE <- P$maxEE
}
maxMismatch <- if (length(args) >= 6) as.integer(args[6]) else 0
cat("params: truncLen =", truncLen, " maxEE =", maxEE, " maxMismatch =", maxMismatch, "\n")

fq_dir <- file.path("<DATA_DIR>/sra", cohort)
R1s <- sort(list.files(fq_dir, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
R2s <- sort(list.files(fq_dir, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))
if (length(R1s) != length(R2s)) stop("双端文件数不匹配")
n_trial <- min(3, length(R1s))
R1s <- R1s[1:n_trial]; R2s <- R2s[1:n_trial]
cat("== 试跑队列:", cohort, "样本数:", n_trial, " truncLen:", truncLen, "==", "\n")
cat("样本:", basename(R1s), "\n")

out_dir <- file.path("<OUT_DIR>", paste0(cohort, "_trial_", truncLen[1], "_", truncLen[2]))
dir.create(out_dir, showWarnings = FALSE)
filt_R1 <- file.path(out_dir, paste0("S", 1:n_trial, "_F.fastq.gz"))
filt_R2 <- file.path(out_dir, paste0("S", 1:n_trial, "_R.fastq.gz"))

out <- filterAndTrim(R1s, filt_R1, R2s, filt_R2,
                     truncLen = truncLen, maxN = 0, maxEE = maxEE,
                     truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = 4)
print(out)
cat("\n== 过滤保留率（min/中位）:",
    sprintf("%.1f%%", min(out[, 2] / out[, 1] * 100)),
    sprintf("%.1f%%", median(out[, 2] / out[, 1] * 100)), "\n")

errF <- learnErrors(filt_R1, multithread = 4)
errR <- learnErrors(filt_R2, multithread = 4)
dada_F <- dada(filt_R1, err = errF, multithread = 4, verbose = FALSE)
dada_R <- dada(filt_R2, err = errR, multithread = 4, verbose = FALSE)

mergers <- mergePairs(dada_F, filt_R1, dada_R, filt_R2, maxMismatch = maxMismatch, verbose = TRUE)
cat("\n== merge 统计（reads 级）==\n")
for (i in seq_along(mergers)) {
  m <- mergers[[i]]
  n_in  <- out[i, 2]
  n_mrg <- if (nrow(m) > 0) sum(m$abundance) else 0
  if (nrow(m) > 0) {
    cat(sprintf("  %s: merge %d / %d (%.1f%%)  重叠长度 中位 %d bp\n",
                basename(R1s[i]), n_mrg, n_in, n_mrg / n_in * 100,
                as.integer(median(m$nmatch))))
  } else {
    cat(sprintf("  %s: 0 merge\n", basename(R1s[i])))
  }
}
seqtab <- makeSequenceTable(mergers)
cat("\nASV 数:", ncol(seqtab), " 序列长度分布:\n")
print(table(nchar(colnames(seqtab))))
saveRDS(seqtab, file.path(out_dir, "seqtab_trial.rds"))
cat("\n== TRIAL_DONE:", cohort, "==\n")
