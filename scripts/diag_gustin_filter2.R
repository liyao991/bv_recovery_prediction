# Unambiguous isolation test for Gustin filter retention (ASCII only)
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages(library(dada2))

fq <- "<DATA_DIR>/sra/PRJNA691964"
R1 <- file.path(fq, "SRR13437347_1.fastq.gz")
R2 <- file.path(fq, "SRR13437347_2.fastq.gz")
od <- "<OUT_DIR>/gustin_diag2"
dir.create(od, showWarnings = FALSE)

show <- function(label, out, out2 = NULL) {
  cat(">>>", label, "\n")
  if (is.null(out2)) print(out) else print(cbind(out, out2))
}

# A: R1 only, no truncLen
a1 <- file.path(od, "a1_F.fastq.gz")
o <- filterAndTrim(R1, a1, maxN = 0, maxEE = 2, compress = TRUE)
show("A: R1 only maxEE=2 (no truncLen)", o)

# B: R1 only, truncLen=150
b1 <- file.path(od, "b1_F.fastq.gz")
o <- filterAndTrim(R1, b1, truncLen = 150, maxN = 0, maxEE = 2, compress = TRUE)
show("B: R1 only truncLen=150 + maxEE=2", o)

# C: R1 only, truncLen=140
c1 <- file.path(od, "c1_F.fastq.gz")
o <- filterAndTrim(R1, c1, truncLen = 140, maxN = 0, maxEE = 2, compress = TRUE)
show("C: R1 only truncLen=140 + maxEE=2", o)

# D: paired, truncLen=150, maxEE=2 (original culprit)
d1 <- file.path(od, "d1_F.fastq.gz"); d2 <- file.path(od, "d1_R.fastq.gz")
o <- filterAndTrim(R1, d1, R2, d2, truncLen = c(150,150), maxN = 0, maxEE = c(2,2), compress = TRUE)
show("D: PAIRED truncLen=c(150,150) maxEE=c(2,2)", o)

# E: paired, no truncLen
e1 <- file.path(od, "e1_F.fastq.gz"); e2 <- file.path(od, "e1_R.fastq.gz")
o <- filterAndTrim(R1, e1, R2, e2, maxN = 0, maxEE = c(2,2), compress = TRUE)
show("E: PAIRED no truncLen maxEE=c(2,2)", o)

# F: R1 only maxEE=2 + truncQ=2 + rm.phix=TRUE, no truncLen
f1 <- file.path(od, "f1_F.fastq.gz")
o <- filterAndTrim(R1, f1, maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE, compress = TRUE)
show("F: R1 only maxEE=2 + truncQ + rm.phix (no truncLen)", o)

cat("=== DONE ===\n")
