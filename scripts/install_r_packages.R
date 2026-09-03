# 线B 16S 处理工具链安装脚本（R + Bioconductor dada2/phyloseq/DECIPHER）
# 使用清华 TUNA 镜像加速。用法：Rscript install_r_packages.R
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 核心：16S 质控/ASV 处理（dada2）、菌群生态分析（phyloseq）、SILVA 分类（DECIPHER）
pkgs <- c("dada2", "phyloseq", "DECIPHER")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("== Installing:", p, "==\n")
    tryCatch(
      BiocManager::install(p, update = FALSE, ask = FALSE),
      error = function(e) cat("FAILED:", p, "->", conditionMessage(e), "\n")
    )
  } else {
    cat("== Already installed:", p, "==\n")
  }
}

# 常用辅助包
cran_pkgs <- c("tidyverse", "data.table", "readxl", "writexl", "openxlsx")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("== Installing CRAN:", p, "==\n")
    tryCatch(
      install.packages(p),
      error = function(e) cat("FAILED:", p, "->", conditionMessage(e), "\n")
    )
  }
}

cat("== 安装完成检查 ==\n")
for (p in c(pkgs, cran_pkgs)) {
  cat(sprintf("%-15s %s\n", p, ifelse(requireNamespace(p, quietly = TRUE), "OK", "MISSING")))
}
