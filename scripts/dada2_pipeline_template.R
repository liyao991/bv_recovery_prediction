# ============================================================================
# 线B · DADA2 16S 处理流水线（模板）
# 用途：统一处理开放队列（Xiao PRJNA302078 / Plummer PRJNA592384 / PRJNA398590）
# 输入：<DATA_DIR>/sra/<PRJNA号>/<run>.fastq.gz   （原 <DATA_DIR> 已迁移，原路径保留 junction 兼容）
# 输出：<OUT_DIR>/<队列>/ 下的 seqtab、taxa、phyloseq 对象（R 对中文路径的文件打开不可用，
#       故计算走 ASCII 路径；完成后由 Python 复制交付物到工作空间 line_B/results/）
# 参数：各队列 truncLen/maxEE 见下方 params 表（由 01_peek_qc.py 实测定稿；
#       Plummer 双端 merge 可行性待小样本试跑校准）
# 质量概况工具：优先用 01_peek_qc.py（Python，避开 Windows 下 ShortRead
#       quality() 矩阵解析 NA 的 bug；01_peek.R 仅作 DADA2 plotQualityProfile 参考）
# ============================================================================
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
suppressPackageStartupMessages({
  library(dada2); library(phyloseq); library(DECIPHER); library(ShortRead)
})

args <- commandArgs(trailingOnly = TRUE)
cohort   <- if (length(args) >= 1) args[1] else "PRJNA302078"   # 队列登录号
out_dir  <- file.path("<OUT_DIR>", cohort)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 队列参数表（2026-08-29 01_peek_qc.py 实测质量曲线定稿）----------------
params <- list(
  PRJNA302078 = list(paired = FALSE, truncLen = 450, maxEE = 2),        # 454 单端，median~510bp，575bp 后质量崩
  PRJNA592384 = list(paired = TRUE,  truncLen = c(240, 200), maxEE = c(2, 2)),  # MiSeq 2×300，R2 差，待 merge 试跑校准
  PRJNA398590 = list(paired = TRUE,  truncLen = c(250, 200), maxEE = c(2, 2))   # MiSeq 2×300，质量更差，待 merge 试跑校准
)
if (!cohort %in% names(params)) stop("未定义参数: ", cohort)
P <- params[[cohort]]
paired <- P$paired

# ---- 1. 定位 fastq ---------------------------------------------------------
fq_dir <- file.path("<DATA_DIR>/sra", cohort)
if (!dir.exists(fq_dir)) stop("找不到目录: ", fq_dir)
fqs <- sort(list.files(fq_dir, pattern = "_1\\.fastq\\.gz$", full.names = TRUE))
if (length(fqs) == 0) fqs <- sort(list.files(fq_dir, pattern = "fastq\\.gz$", full.names = TRUE))
if (length(fqs) == 0) stop("目录中没有 fastq.gz")
sample_names <- sub("_1\\.fastq\\.gz$|_2\\.fastq\\.gz$|\\.fastq\\.gz$", "", basename(fqs))

if (paired) {
  fqs_R1 <- fqs; fqs_R2 <- sort(list.files(fq_dir, pattern = "_2\\.fastq\\.gz$", full.names = TRUE))
  if (length(fqs_R1) != length(fqs_R2)) stop("双端文件数不匹配")
}

# ---- 2. 过滤与修剪（truncLen/maxEE 来自队列参数表）-------------------------
filt_dir <- file.path(out_dir, "filtered"); dir.create(filt_dir, showWarnings = FALSE)
filt_R1 <- file.path(filt_dir, paste0(sample_names, "_F.fastq.gz"))
if (paired) filt_R2 <- file.path(filt_dir, paste0(sample_names, "_R.fastq.gz"))

if (paired) {
  out <- filterAndTrim(fqs_R1, filt_R1, fqs_R2, filt_R2,
                       truncLen = P$truncLen, maxN = 0, maxEE = P$maxEE,
                       truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
} else {
  out <- filterAndTrim(fqs, filt_R1, truncLen = P$truncLen, maxN = 0, maxEE = P$maxEE,
                       truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
}
print(out); saveRDS(out, file.path(out_dir, "filter_stats.rds"))

# ---- 3. 学习误差与去噪 ------------------------------------------------------
errF <- learnErrors(filt_R1, multithread = TRUE)
if (paired) errR <- learnErrors(filt_R2, multithread = TRUE)
dada_F <- dada(filt_R1, err = errF, multithread = TRUE, pool = FALSE)
if (paired) dada_R <- dada(filt_R2, err = errR, multithread = TRUE, pool = FALSE)

# ---- 4. 合并（双端） / 构建计数表（单端） -----------------------------------
if (paired) {
  mergers <- mergePairs(dada_F, filt_R1, dada_R, filt_R2, verbose = TRUE)
  seqtab  <- makeSequenceTable(mergers)
} else {
  seqtab <- makeSequenceTable(dada_F)
}
cat("ASV 数: ", ncol(seqtab), "\n")
saveRDS(seqtab, file.path(out_dir, "seqtab.rds"))

# ---- 5. 去嵌合体 ------------------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
cat("去嵌合后 ASV 数: ", ncol(seqtab.nochim), "\n")
saveRDS(seqtab.nochim, file.path(out_dir, "seqtab_nochim.rds"))

# ---- 6. 分类注释（SILVA v138.1） ------------------------------------------
# 来源：SILVA 官方文件服务器（zenodo 不可达时的替代源，2026-08-28 已验证可达）
#   https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/
#   SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz
# 注意：官方库 header 为 ">accession 域;门;纲;目;科;属;种"（7 层分号分类），
#       assignTaxonomy 直接输出到种级，无需 addSpecies。
taxa_file <- "<DATA_DIR>/silva/SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"
if (file.exists(taxa_file)) {
  taxa <- assignTaxonomy(seqtab.nochim, taxa_file, multithread = TRUE)
  # 若后续拿到 DADA2 官方物种文件（zenodo 恢复后），可取消注释以 refine 种级：
  # sp_file <- "<DATA_DIR>/silva/silva_species_assignment_v138.1.fa.gz"
  # if (file.exists(sp_file)) taxa <- addSpecies(taxa, sp_file)
  saveRDS(taxa, file.path(out_dir, "taxa.rds"))
} else {
  warning("SILVA 训练集缺失，跳过分类注释：", taxa_file)
}

# ---- 7. 导出 phyloseq 对象与计数表 -----------------------------------------
if (exists("taxa")) {
  ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), tax_table(taxa))
  saveRDS(ps, file.path(out_dir, "phyloseq.rds"))
}
write.csv(seqtab.nochim, file.path(out_dir, "asv_table.csv"))
cat("=== 完成:", cohort, "===\n")
