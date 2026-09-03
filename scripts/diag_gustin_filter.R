# Diagnose Gustin filterAndTrim low retention (all ASCII, UTF-8 no)
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages(library(dada2))

fq <- "<DATA_DIR>/sra/PRJNA691964"
R1 <- file.path(fq, "SRR13437347_1.fastq.gz")
R2 <- file.path(fq, "SRR13437347_2.fastq.gz")
od <- "<OUT_DIR>/gustin_diag"
dir.create(od, showWarnings = FALSE)

t1 <- file.path(od, "t1_F.fastq.gz"); t2 <- file.path(od, "t1_R.fastq.gz")
o1 <- filterAndTrim(R1, t1, R2, t2, maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE)
print(o1); cat("--- t1: full params ---\n")

t3 <- file.path(od, "t3_F.fastq.gz"); t4 <- file.path(od, "t3_R.fastq.gz")
o2 <- filterAndTrim(R1, t3, R2, t4, maxN = 0, maxEE = c(2,2), rm.phix = FALSE, compress = TRUE)
print(o2); cat("--- t3: no truncQ, no rm.phix ---\n")

t5 <- file.path(od, "t5_F.fastq.gz"); t6 <- file.path(od, "t5_R.fastq.gz")
o3 <- filterAndTrim(R1, t5, R2, t6, maxN = 0, maxEE = c(2,2), compress = TRUE)
print(o3); cat("--- t5: maxEE only ---\n")

t7 <- file.path(od, "t7_F.fastq.gz"); t8 <- file.path(od, "t7_R.fastq.gz")
o4 <- filterAndTrim(R1, t7, R2, t8, truncLen = c(150,150), maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE)
print(o4); cat("--- t7: + truncLen 150 ---\n")

# quality encoding auto-detect check
suppressPackageStartupMessages(library(ShortRead))
q <- readFastq(R1, n = 1000)
cat("qualityType guess (dada2 auto):\n")
print(dada2:::guessQualityEncoding(q, "fastq"))
