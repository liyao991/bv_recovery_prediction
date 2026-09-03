# -*- coding: utf-8 -*-
"""
00_bench_hardware.py —— 设备能力基线评测
用途：评估本机是否足以支撑 线A(WGS) / 线B(DADA2+ML) / 线B'(shotgun) 三线计算负载。
用法: py_env python 00_bench_hardware.py
注意：内存分配目标 = 可用内存 40%（上限 10GB）；磁盘测速写 256MB 临时文件后删除。
"""
import os
import time
import tempfile
import platform
import subprocess
import psutil
import numpy as np


def cpu_name():
    try:
        out = subprocess.check_output(["wmic", "cpu", "get", "name"],
                                      text=True, timeout=10)
        lines = [l.strip() for l in out.strip().splitlines() if l.strip()]
        return lines[1] if len(lines) > 1 else "?"
    except Exception:
        try:
            out = subprocess.check_output(
                ["powershell", "-NoProfile", "-Command",
                 "(Get-CimInstance Win32_Processor).Name"],
                text=True, timeout=15)
            return out.strip() or "?"
        except Exception:
            return "unknown"


print("=" * 64)
print("设备能力评估（线A WGS / 线B DADA2+ML / 线B' shotgun）")
print("=" * 64)

print("\n[CPU]")
print(f"  型号: {cpu_name()}")
print(f"  物理核心: {psutil.cpu_count(logical=False)}  逻辑线程: {psutil.cpu_count(logical=True)}")
freq = psutil.cpu_freq()
if freq:
    print(f"  当前频率: {freq.current:.0f} MHz (最大 {freq.max:.0f} MHz)")
print(f"  空闲负载: {psutil.cpu_percent(interval=1.0)}%")

print("\n[内存]")
vm = psutil.virtual_memory()
print(f"  物理总内存: {vm.total/1e9:.1f} GB")
print(f"  当前可用:   {vm.available/1e9:.1f} GB ({vm.percent}% 已用)")
sm = psutil.swap_memory()
print(f"  页面文件:   总 {sm.total/1e9:.1f} GB / 可用 {sm.free/1e9:.1f} GB")

print("\n[磁盘]")
for p in psutil.disk_partitions(all=False):
    if p.fstype == "":
        continue
    try:
        u = psutil.disk_usage(p.mountpoint)
        print(f"  {p.device} 挂载 {p.mountpoint} ({p.fstype}) 总 {u.total/1e9:.0f} GB / 剩余 {u.free/1e9:.0f} GB")
    except Exception:
        pass

# ---- 内存分配压力测试（模拟 DADA2 去噪 / 大矩阵驻留）----
print("\n[压力测试 1] 大数组内存分配（模拟 DADA2 峰值驻留）")
avail = vm.available
target = min(avail * 0.4, 10e9)
n = int(target / 4)
allocated = 0
for factor, label in [(1.0, "目标"), (0.5, "半量")]:
    n_try = int(n * factor)
    try:
        t0 = time.time()
        A = np.ones(n_try, dtype=np.float32)
        A[:: 1 << 20] = 1  # touch 一页防懒分配
        dt = time.time() - t0
        allocated = A.nbytes
        print(f"  ✅ 成功分配 {allocated/1e9:.1f} GB float32 数组（耗时 {dt:.1f}s）")
        # 做一次快速规约验证可计算
        t0 = time.time()
        s = float(A.sum())
        print(f"     全数组求和 {time.time()-t0:.2f}s, 校验和={s:.1f}")
        del A
        break
    except MemoryError:
        print(f"  ❌ 分配 {n_try*4/1e9:.1f} GB 失败 (MemoryError)")
print(f"  → 结论: 峰值可安全驻留约 {max(allocated, 0)/1e9:.1f} GB 数据")

# ---- 浮点算力（模拟 ML/DADA2 迭代运算）----
print("\n[压力测试 2] 浮点算力（1000×1000 矩阵乘 × 5 次）")
A = np.random.rand(1000, 1000).astype(np.float64)
t0 = time.time()
for _ in range(5):
    C = A @ A
dt = time.time() - t0
flops = 5 * 2 * 1000**3
print(f"  {dt:.2f}s → {flops/dt/1e9:.0f} GFLOPS")
del A, C

# ---- ML 训练负载（模拟 Line B 特征工程后的建模）----
print("\n[压力测试 3] ML 训练（模拟 Line B：5万样本 × 100 特征）")
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
import xgboost as xgb

rng = np.random.default_rng(42)
X = rng.random((50000, 100), dtype=np.float32)
y = (X[:, 0] + X[:, 1] * 2 > 1.5).astype(int)

t0 = time.time()
lr = LogisticRegression(max_iter=300, solver="liblinear").fit(X, y)
print(f"  LASSO-LR     {time.time()-t0:.1f}s")

t0 = time.time()
rf = RandomForestClassifier(n_estimators=100, n_jobs=-1).fit(X, y)
print(f"  RF 100 树    {time.time()-t0:.1f}s")

t0 = time.time()
xgm = xgb.XGBClassifier(n_estimators=200, n_jobs=-1, tree_method="hist",
                        max_depth=6, learning_rate=0.1).fit(X, y)
print(f"  XGB 200 轮   {time.time()-t0:.1f}s")

# ---- 磁盘顺序写测速（256MB，写完即删，避免干扰下载）----
print("\n[压力测试 4] 磁盘顺序写测速（256 MB，写完删除）")
test_dir = tempfile.gettempdir()  # 系统临时目录，勿写入 WPS 同步的工作空间
fpath = os.path.join(test_dir, "bench_tmp.bin")
buf = np.random.bytes(1 << 20)
try:
    t0 = time.time()
    with open(fpath, "wb") as f:
        for _ in range(256):
            f.write(buf)
    dt = time.time() - t0
    print(f"  256 MB 写入 {dt:.2f}s → {256/dt:.0f} MB/s")
except Exception as e:
    print(f"  写入失败: {e}")
finally:
    try:
        if os.path.exists(fpath):
            os.remove(fpath)
    except OSError:
        pass

print("\n" + "=" * 64)
print("评测完成")
print("=" * 64)
