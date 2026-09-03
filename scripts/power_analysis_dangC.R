# power_analysis_dangC.R
# 档 C 方案 A 功效分析：入组后 1 个月内复发（3 例）vs 未复发（65 例）
library(pROC)

cat("=== 档 C 方案 A：1 个月内复发 vs 未复发（Plummer PRJNA592384）===\n")
cat(sprintf("n=68  复发=3  未复发=65  事件率=%.3f\n\n", 3/68))

# 最小可检测 AUC
pt <- power.roc.test(ncases = 3, ncontrols = 65,
                     sig.level = 0.05, power = 0.8, alternative = "two.sided")
cat(sprintf("最小可检测 AUC (alpha=0.05, power=0.8) = %.3f\n", pt$auc))

# 检出 0.65 / 0.70 所需 n（保持 3:65 比例 ≈ 事件率 4.4%）
for (target in c(0.65, 0.70, 0.80)) {
  n_need <- NA
  tryCatch({
    lo <- 10; hi <- 5000
    for (iter in 1:60) {
      mid <- round((lo + hi) / 2)
      nc <- max(1, round(mid * 3 / 68))
      nct <- mid - nc
      if (nct < 1) { lo <- mid + 1; next }
      p2 <- power.roc.test(ncases = nc, ncontrols = nct,
                           sig.level = 0.05, power = 0.8, alternative = "two.sided")
      if (p2$auc <= target) { hi <- mid } else { lo <- mid + 1 }
    }
    n_need <- hi
  }, error = function(e) n_need <<- NA)
  cat(sprintf("检出 AUC=%.2f 需 n ≈ %s（保持 3:65 比例）\n", target, n_need))
}

# 反向：现有 3/65 能检出的 AUC 已经极差，补充：若事件数=16（方案 B 任意复发）会怎样
cat("\n=== 对照：方案 B（任意随访复发 16 例）的功效 ===\n")
pt_b <- power.roc.test(ncases = 16, ncontrols = 52,
                       sig.level = 0.05, power = 0.8, alternative = "two.sided")
cat(sprintf("方案 B n=68 复发=16 未复发=52 最小可检测 AUC = %.3f\n", pt_b$auc))
