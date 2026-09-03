# export_398590_features.R - PRJNA398590 (StepUp) 特征导出（绕过 phyloseq 段错误）
# 背景：phyloseq.rds 加载时 segfault（Exit 139），故改用原始文件：
#   - taxa.rds     : ASV x 7 分类矩阵（Kingdom/Phylum/Class/Order/Family/Genus/Species）
#   - asv_table.csv: samples x ASV 计数矩阵（首列样本名，表头为 ASV 序列）
# 逻辑完全对齐 export_genus_matrix.R / export_species_matrix.R / export_alpha_family.R
# 输出（与其它队列同名格式）：
#   <DATA_DIR>/feat/PRJNA398590_genus_rel.csv
#   <DATA_DIR>/feat/PRJNA398590_species_rel.csv
#   <DATA_DIR>/feat/PRJNA398590_alpha_family.csv
#   <DATA_DIR>/feat/PRJNA398590_alpha.csv
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages(library(vegan))

co <- "PRJNA398590"
out_dir <- "<DATA_DIR>/feat"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 读取原始文件（不依赖 phyloseq）----
tax <- readRDS(file.path("<OUT_DIR>", co, "taxa.rds"))   # ASV x 7 matrix
otu <- read.csv(file.path("<OUT_DIR>", co, "asv_table.csv"),
                row.names = 1, check.names = FALSE)          # samples x ASV

# 对齐：tax 行名（ASV 序列）必须与 otu 列名一致
stopifnot(all(colnames(otu) %in% rownames(tax)))
tax <- as.data.frame(tax[colnames(otu), , drop = FALSE], stringsAsFactors = FALSE)

cat("cohort", co, ": samples", nrow(otu), "x ASVs", ncol(otu), "\n")

# ---- 属级聚合（对齐 export_genus_matrix.R）----
gen <- as.character(tax$Genus)
gen[is.na(gen)] <- "Unassigned"
gmat <- as.data.frame(t(apply(otu, 1, function(row) tapply(row, gen, sum))))
gmat <- gmat[, order(colSums(gmat), decreasing = TRUE)]
gmat_rel <- gmat / rowSums(gmat)
write.csv(gmat_rel, file.path(out_dir, paste0(co, "_genus_rel.csv")), row.names = TRUE)

# ---- 种级聚合（对齐 export_species_matrix.R）----
sp <- as.character(tax$Species)
sp[is.na(sp) | sp %in% c("", "NA")] <- "unclassified"
is_sp <- !grepl("uncultured|unclassified| sp\\.$", sp, ignore.case = TRUE)
feat <- ifelse(is_sp, paste(gen, sp, sep = " "), paste(gen, "unclassified", sep = " "))
feat <- gsub("\\s+", " ", feat)
smat <- as.data.frame(t(apply(otu, 1, function(row) tapply(row, feat, sum))))
smat <- smat[, order(colSums(smat), decreasing = TRUE)]
smat_rel <- smat / rowSums(smat)
write.csv(smat_rel, file.path(out_dir, paste0(co, "_species_rel.csv")), row.names = TRUE)

# ---- α 多样性家族（对齐 export_alpha_family.R）----
agg_by <- function(mat, mapping) {
  lbl <- as.character(mapping)
  lbl[is.na(lbl) | lbl == ""] <- "Unclassified"
  keep <- rownames(mat) %in% names(mapping)
  mat <- mat[keep, , drop = FALSE]
  lbl <- lbl[keep]
  if (any(duplicated(lbl))) {
    out <- rowsum(mat, group = lbl, reorder = TRUE)
  } else {
    out <- mat
    rownames(out) <- lbl
  }
  out
}

alpha_family <- function(mat) {
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
    D <- tryCatch(diversity(x, index = "simpson"), error = function(e) NA_real_)
    invD <- tryCatch(diversity(x, index = "invsimpson"), error = function(e) NA_real_)
    out[i, "shannon"] <- H
    out[i, "simpson"] <- D
    out[i, "invsimpson"] <- invD
    out[i, "pielou"] <- if (S > 1 && is.finite(H)) H / log(S) else if (S == 1) 1 else NA_real_
  }
  as.data.frame(out)
}

otu_mat <- as.matrix(otu)      # samples x ASV
otu_t <- t(otu_mat)            # ASV x samples
smp <- rownames(otu_mat)

g_map <- setNames(gen, rownames(tax))   # ASV -> genus
s_map <- setNames(feat, rownames(tax))  # ASV -> species-feature

a_asv <- alpha_family(otu_mat)
otu_g <- agg_by(otu_t, g_map)
a_gen <- alpha_family(t(otu_g))
otu_spp <- agg_by(otu_t, s_map)
a_spp <- alpha_family(t(otu_spp))

cols <- c("observed", "chao1", "shannon", "simpson", "invsimpson", "pielou")
res <- data.frame(sample = smp, row.names = NULL, check.names = FALSE)
for (lv in c("asv", "genus", "species")) {
  a <- switch(lv, asv = a_asv, genus = a_gen, species = a_spp)
  for (m in cols) res[[paste0(lv, "_", m)]] <- a[[m]]
}
write.csv(res, file.path(out_dir, paste0(co, "_alpha_family.csv")), row.names = FALSE)

# ---- 简单 alpha.csv（对齐 export_genus_matrix.R 的输出列）----
shannon <- diversity(otu_mat, index = "shannon")
n_asv <- specnumber(otu_mat)
chao1 <- sapply(seq_len(nrow(otu_mat)), function(i)
  tryCatch(estimateR(as.numeric(otu_mat[i, ]))["S.chao1"], error = function(e) NA_real_))
meta_df <- data.frame(sample = smp, shannon = shannon, chao1 = chao1, n_asv = n_asv,
                      stringsAsFactors = FALSE)
write.csv(meta_df, file.path(out_dir, paste0(co, "_alpha.csv")), row.names = FALSE)

cat("DONE", co, ": samples", nrow(otu_mat),
    "genera", ncol(gmat_rel), "species-features", ncol(smat_rel), "\n")
