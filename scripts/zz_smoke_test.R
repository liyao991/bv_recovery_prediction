# 测试 saveRDS / write.csv 到含中文路径（输出目录）是否可行
suppressWarnings(Sys.setlocale("LC_CTYPE", "Chinese (Simplified)_China.936"))
out_dir <- "<REPO_DIR>/results/PRJNA302078"
cat("输出目录:", out_dir, "\n")
ok <- dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("dir.create 结果:", ok, " 存在:", dir.exists(out_dir), "\n")
f1 <- file.path(out_dir, "zz_test.rds")
saveRDS(list(a = 1), f1)
cat("saveRDS 到中文路径:", file.exists(f1), "\n")
f2 <- file.path(out_dir, "zz_test.csv")
write.csv(data.frame(x = 1), f2)
cat("write.csv 到中文路径:", file.exists(f2), "\n")
# 验证正常化后仍可读回
x <- readRDS(f1); cat("readRDS 回读:", x$a == 1, "\n")
unlink(c(f1, f2))
cat("OK 中文输出路径测试通过\n")
