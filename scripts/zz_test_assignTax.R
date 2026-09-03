# test_assignTaxonomy.R - minimal reproduction of Memory allocation failed
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages(library(dada2))

# load Xiao seqtab, take first 50 ASVs
st <- readRDS("<OUT_DIR>/PRJNA302078/seqtab_nochim.rds")
st50 <- st[, 1:min(50, ncol(st)), drop = FALSE]
cat("test seqtab:", dim(st50), "\n")

# small trainset (5000 refs)
taxa_small <- "<DATA_DIR>/silva/test_5000.fa.gz"
cat("assignTaxonomy with 5000-ref trainset, multithread=FALSE ...\n")
t0 <- Sys.time()
res <- tryCatch(
  assignTaxonomy(st50, taxa_small, multithread = FALSE),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
cat("elapsed:", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")
if (!is.null(res)) {
  print(dim(res))
  print(res[1:3, , drop = FALSE])
  cat("SMALL OK\n")
}
