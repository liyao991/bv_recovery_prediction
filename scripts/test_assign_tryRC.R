# test_assign_tryRC.R - verify tryRC=TRUE fixes taxonomy (all-ASCII)
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({ library(dada2) })
st <- readRDS("<OUT_DIR>/PRJNA592384/seqtab_nochim.rds")
top <- sort(colSums(st), decreasing = TRUE)[1:10]
mini <- st[, names(top), drop = FALSE]
cat("testing", ncol(mini), "ASVs\n")
t1 <- assignTaxonomy(mini, "<DATA_DIR>/silva/silva_138.1_train_set_bac_arch.fa.gz",
                     multithread = TRUE, tryRC = TRUE)
cat("non-NA per level:\n")
print(colSums(!is.na(t1)))
cat("top ASV assignments:\n")
print(t1[1:3, , drop = FALSE])
cat("== TRYRC_TEST_DONE ==\n")
