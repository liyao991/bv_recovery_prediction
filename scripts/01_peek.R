# 01_peek.R —— 查看队列 fastq 质量概况，用于确定 DADA2 truncLen 参数
# 用法: Rscript 01_peek.R <PRJNA号> [样本数上限]
suppressPackageStartupMessages({ library(dada2); library(ShortRead) })
args <- commandArgs(trailingOnly = TRUE)
cohort <- if (length(args) >= 1) args[1] else "PRJNA302078"
n_max  <- if (length(args) >= 2) as.integer(args[2]) else 3
fq_dir <- file.path("<DATA_DIR>/sra", cohort)
fqs <- sort(list.files(fq_dir, pattern = "fastq\\.gz$", full.names = TRUE))
cat("共", length(fqs), "个 fastq，抽查前", min(n_max, length(fqs)), "个\n")
for (f in head(fqs, min(n_max, length(fqs)))) {
  cat("=====", basename(f), "=====\n")
  sr <- readFastq(f, withIds = FALSE, n = 50000)
  q  <- quality(sr)
  cat("reads(≤5万):", length(sr), " 平均读长:", round(mean(width(sr))), "bp\n")
  # 平均质量分（na.rm：454/旧平台质量字符偶有异常）
  qmat <- as(q, "matrix"); qmat <- as.numeric(qmat)
  cat("平均Q:", round(mean(qmat, na.rm = TRUE), 1),
      "  <Q20比例:", round(mean(qmat < 20, na.rm = TRUE) * 100, 1), "%",
      "  NA数:", sum(is.na(qmat)), "\n")
  # 质量字符 ASCII 范围（判断 Phred 编码：+33 正常约 33-73）
  chars <- as.integer(charToRaw(paste(as.character(q), collapse = "")))
  cat("质量字符ASCII:", min(chars), "-", max(chars), "\n")
}
# 画质量曲线存图
if (length(fqs) > 0) {
  pdf(file.path("<REPO_DIR>/results", paste0(cohort, "_quality.pdf")), width = 10, height = 6)
  plotQualityProfile(head(fqs, min(n_max, length(fqs))))
  dev.off()
  cat("质量图已存 results/", cohort, "_quality.pdf\n", sep = "")
}
