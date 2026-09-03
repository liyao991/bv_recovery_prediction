# c6_metafor.R
# #C6 task 1: random-effects meta-analysis of richness Spearman rho across 3 cohorts
# (baseline CST4 subgroup, binary outcome LB>=50%)
# Input rho corrected to baseline-CST4 subgroup values:
#   Xiao +0.04533 (n=54), Plummer +0.07126 (n=43), Gustin -0.4537 (n=33)
# Gustin rho corrected to LB>=0.5 outcome (was -0.4977 under >log1p(0.5) typo).
# Fisher-z transform, SE = 1/sqrt(n-3), random-effects (REML), back-transform.

library(metafor)

rho <- c(0.04533, 0.07126, -0.4537)
n   <- c(54, 43, 33)
lab <- c("Xiao", "Plummer", "Gustin")

z  <- atanh(rho)
vi <- 1 / (n - 3)

# random-effects (REML, metafor default)
res <- rma(yi = z, vi = vi, method = "REML")
res_dl <- rma(yi = z, vi = vi, method = "DL")

pooled_z  <- res$b[1, 1]
pooled_rho <- tanh(pooled_z)
ci_z  <- c(res$ci.lb, res$ci.ub)
ci_rho <- tanh(ci_z)
tau2  <- res$tau2
I2    <- res$I2
Q     <- res$QE
Qp    <- res$QEp
se_z  <- res$se
pval  <- res$pval

# per-cohort Fisher-z CI for plotting
rows <- data.frame(
  cohort = lab,
  n = n,
  rho = rho,
  z = z,
  z_se = sqrt(vi),
  stringsAsFactors = FALSE
)
rows$z_lo <- rows$z - 1.96 * rows$z_se
rows$z_hi <- rows$z + 1.96 * rows$z_se
rows$rho_lo <- tanh(rows$z_lo)
rows$rho_hi <- tanh(rows$z_hi)

write.csv(rows, "<OUT_DIR>/modeling/c6_metafor_cohort.csv", row.names = FALSE)

pooled_row <- data.frame(
  cohort = "Pooled (REML)",
  n = sum(n),
  rho = pooled_rho,
  z = pooled_z,
  z_se = se_z,
  z_lo = ci_z[1],
  z_hi = ci_z[2],
  rho_lo = ci_rho[1],
  rho_hi = ci_rho[2],
  stringsAsFactors = FALSE
)
write.csv(pooled_row, "<OUT_DIR>/modeling/c6_metafor_pooled.csv", row.names = FALSE)

cat(sprintf("pooled_rho_REML=%.4f\n", pooled_rho))
cat(sprintf("pooled_rho_ci_lo=%.4f\n", ci_rho[1]))
cat(sprintf("pooled_rho_ci_hi=%.4f\n", ci_rho[2]))
cat(sprintf("tau2=%.5f\n", tau2))
cat(sprintf("I2=%.2f\n", I2))
cat(sprintf("Q=%.3f\n", Q))
cat(sprintf("Q_df=%d\n", res$k - 1))
cat(sprintf("Q_p=%.4g\n", Qp))
cat(sprintf("pooled_p=%.4g\n", pval))
cat(sprintf("pooled_rho_DL=%.4f\n", tanh(res_dl$b[1, 1])))
