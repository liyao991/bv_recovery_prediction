# power_analysis.R
# #C2 任务一：功效分析（power analysis）
# 用 pROC::power.roc.test (Hanley-McNeil) 计算每个 CST4 亚组的最小可检测 AUC，
# 以及检出 AUC=0.65 / 0.70 所需的最小样本量。
library(pROC)

# 三队列 CST4 亚组的 n、复发数、康复数（来自 results_richness_direction.csv）
# 复发=结局1（未康复 LB<50%），康复=结局0（康复 LB>=50%）
cohorts <- data.frame(
  cohort    = c("Xiao",    "Plummer",  "Gustin"),
  n         = c(54,        43,         33),
  n_cases   = c(27,        39,         12),   # 复发（结局=1）
  n_controls= c(27,        4,          21)    # 康复（结局=0）
)

cat("=== 各队列 CST4 亚组 三件套 ===\n")
for (i in seq_len(nrow(cohorts))) {
  cc <- cohorts[i, ]
  cat(sprintf("%-8s n=%d  复发=%d  康复=%d  事件率=%.2f\n",
              cc$cohort, cc$n, cc$n_cases, cc$n_controls, cc$n_cases/cc$n))
}

cat("\n=== 任务一：最小可检测 AUC（alpha=0.05, power=0.8, 双侧 vs 0.5）===\n")
results <- data.frame()
for (i in seq_len(nrow(cohorts))) {
  cc <- cohorts[i, ]
  # power.roc.test 需 ncases/ncontrols
  pt <- power.roc.test(ncases = cc$n_cases, ncontrols = cc$n_controls,
                       sig.level = 0.05, power = 0.8, alternative = "two.sided")
  # pt$auc 是能检出的备择 AUC（min detectable）
  min_auc <- pt$auc

  # 反向：给定 AUC=0.65 / 0.70，各需多大 n（假设 case:control 比例保持）
  ratio <- cc$n_cases / cc$n_controls
  for (target in c(0.65, 0.70)) {
    n_need <- NA
    # 二分搜索总 n（保持比例）
    tryCatch({
      lo <- 10; hi <- 2000
      for (iter in 1:50) {
        mid <- round((lo + hi) / 2)
        nc <- round(mid * cc$n_cases / cc$n)  # cases
        nct <- mid - nc                        # controls
        if (nct < 1) { lo <- mid + 1; next }
        p2 <- power.roc.test(ncases = nc, ncontrols = nct,
                             sig.level = 0.05, power = 0.8,
                             alternative = "two.sided")
        # 若 mid 样本量对应的最小可检测 AUC <= target，说明 mid 够大
        if (p2$auc <= target) { hi <- mid } else { lo <- mid + 1 }
      }
      n_need <- hi
    }, error = function(e) n_need <<- NA)
    results <- rbind(results, data.frame(
      cohort = cc$cohort, n = cc$n, n_cases = cc$n_cases, n_controls = cc$n_controls,
      event_rate = round(cc$n_cases/cc$n, 3),
      min_detect_auc = round(min_auc, 3),
      target_auc = target, n_needed = n_need))
  }
}

cat("\n=== 结果表 ===\n")
print(results, row.names = FALSE)

write.csv(results, "<OUT_DIR>/modeling/results_power_analysis.csv", row.names = FALSE)
cat("\n已写入 <OUT_DIR>/modeling/results_power_analysis.csv\n")
