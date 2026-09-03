# 环境验证：用 dada2 自带示例数据跑通最小完整流程
suppressPackageStartupMessages({
  library(dada2); library(phyloseq); library(DECIPHER)
})

cat("== 版本 ==", "\n")
cat("R:", R.version.string, "\n")
cat("dada2:", as.character(packageVersion("dada2")), "\n")
cat("phyloseq:", as.character(packageVersion("phyloseq")), "\n")
cat("DECIPHER:", as.character(packageVersion("DECIPHER")), "\n")

ext <- system.file("extdata", package = "dada2")
path <- file.path(ext, "sam1F.fastq.gz")
cat("\n== 1. filterAndTrim ==", "\n")
filt <- tempfile(fileext = ".fastq.gz")
out <- filterAndTrim(path, filt, truncLen = 200, maxEE = 2, verbose = TRUE)
print(out)

cat("\n== 2. learnErrors + dada ==", "\n")
err <- learnErrors(filt, nbases = 1e6, multithread = TRUE)
dd <- dada(filt, err = err, verbose = TRUE)
cat("denoised 序列数:", sum(dd$denoised), "\n")

cat("\n== 3. makeSequenceTable ==", "\n")
st <- makeSequenceTable(dd)
cat("样本 x ASV 维度:", dim(st), "\n")

cat("\n== 4. removeBimeraDenovo ==", "\n")
st2 <- removeBimeraDenovo(st, method = "consensus", verbose = TRUE)
cat("去嵌合后维度:", dim(st2), "\n")

cat("\n== 5. assignTaxonomy (SILVA 示例训练集) ==", "\n")
tax <- assignTaxonomy(st2, file.path(ext, "example_train_set.fa.gz"), multithread = TRUE)
tax <- addSpecies(tax, file.path(ext, "example_species_assignment.fa.gz"))
print(tax[1, , drop = FALSE])

cat("\n== 6. phyloseq 构建 ==", "\n")
ps <- phyloseq(otu_table(st2, taxa_are_rows = FALSE), tax_table(tax))
cat("phyloseq OK, taxa 数:", ntaxa(ps), "\n")

cat("\n✅ 环境验证通过：dada2 → phyloseq → DECIPHER 完整链路可用\n")
