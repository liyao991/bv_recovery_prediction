# test_assign_mini.R - minimal assignTaxonomy test to isolate NA cause
# Test: 10 real ASVs from Plummer1 + full trainset -> check non-NA rate
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({ library(dada2) })

st <- readRDS("<OUT_DIR>/PRJNA592384/seqtab_nochim.rds")
cat("seqtab dim:", dim(st), "\n")
# take 10 most abundant ASVs
top <- sort(colSums(st), decreasing = TRUE)[1:10]
mini <- st[, names(top), drop = FALSE]
cat("testing", ncol(mini), "ASVs\n")

t1 <- assignTaxonomy(mini, "<DATA_DIR>/silva/silva_138.1_train_set_bac_arch.fa.gz",
                     multithread = TRUE)
cat("non-NA counts per level:\n")
print(colSums(!is.na(t1)))
cat("any assignment:\n")
print(t1)
cat("== MINI_TEST_DONE ==\n")
