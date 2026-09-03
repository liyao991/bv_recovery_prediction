# export_alpha_family.R - 三队列 α 多样性家族导出（ASV/属/种 三级 × 6 指标）
# 输入: <OUT_DIR>/{cohort}/phyloseq.rds
# 输出: <DATA_DIR>/feat/{cohort}_alpha_family.csv
#   sample 行 × 列: {level}_{metric}: asv_observed, asv_chao1, asv_shannon, asv_simpson, asv_invsimpson, asv_pielou, genus_*, species_*
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
})

COHORTS <- c("PRJNA302078", "PRJNA592384", "PRJNA691964")  # Xiao / Plummer / Gustin

agg_by <- function(otu, mapping) {
  # otu: taxa × sample 计数矩阵; mapping: 每个 taxa 的层级标签（含 NA）
  lbl <- as.character(mapping)
  lbl[is.na(lbl) | lbl == ""] <- "Unclassified"
  keep <- rownames(otu) %in% names(mapping)
  otu <- otu[keep, , drop = FALSE]
  lbl <- lbl[keep]
  if (any(duplicated(lbl))) {
    out <- rowsum(otu, group = lbl, reorder = TRUE)
  } else {
    out <- otu
    rownames(out) <- lbl
  }
  out
}

alpha_family <- function(mat) {
  # mat: sample × taxa 整数计数矩阵
  n_s <- nrow(mat)
  out <- matrix(NA_real_, n_s, 6,
                dimnames = list(rownames(mat),
                                c("observed", "chao1", "shannon", "simpson", "invsimpson", "pielou")))
  for (i in seq_len(n_s)) {
    x <- as.numeric(mat[i, ])
    S <- sum(x > 0)
    out[i, "observed"] <- S
    if (S == 0) next
    out[i, "chao1"] <- tryCatch(estimateR(x)["S.chao1"], error = function(e) NA_real_)
    H <- tryCatch(diversity(x, index = "shannon"), error = function(e) NA_real_)
    D <- tryCatch(diversity(x, index = "simpson"), error = function(e) NA_real_)  # 1-D
    invD <- tryCatch(diversity(x, index = "invsimpson"), error = function(e) NA_real_)
    out[i, "shannon"] <- H
    out[i, "simpson"] <- D
    out[i, "invsimpson"] <- invD
    out[i, "pielou"] <- if (S > 1 && is.finite(H)) H / log(S) else if (S == 1) 1 else NA_real_
  }
  as.data.frame(out)
}

for (co in COHORTS) {
  ps <- readRDS(file.path("<OUT_DIR>", co, "phyloseq.rds"))
  otu <- as.matrix(otu_table(ps))          # taxa × sample（计数）
  if (taxa_are_rows(ps)) otu <- otu else otu <- t(otu)
  smp <- sample_names(ps)
  # 转 sample × taxa
  otu_s <- t(otu)
  rownames(otu_s) <- smp

  tax <- as.data.frame(tax_table(ps))
  g_map <- setNames(tax$Genus, rownames(tax))
  s_map <- setNames(tax$Species, rownames(tax))

  a_asv <- alpha_family(otu_s)
  otu_g <- agg_by(otu, g_map)
  a_gen <- alpha_family(t(otu_g))
  otu_spp <- agg_by(otu, s_map)
  a_spp <- alpha_family(t(otu_spp))

  cols <- c("observed", "chao1", "shannon", "simpson", "invsimpson", "pielou")
  res <- data.frame(sample = smp, row.names = NULL, check.names = FALSE)
  for (lv in c("asv", "genus", "species")) {
    a <- switch(lv, asv = a_asv, genus = a_gen, species = a_spp)
    for (m in cols) res[[paste0(lv, "_", m)]] <- a[[m]]
  }
  write.csv(res, file.path("<DATA_DIR>/feat", paste0(co, "_alpha_family.csv")), row.names = FALSE)
  cat("=== ", co, " done: ", nrow(res), " samples x ", ncol(res) - 1, " features ===\n", sep = "")
  cat("  示例: ", paste(res[1, 2:7], collapse = ", "), "\n", sep = "")
}
cat("ALL DONE\n")
