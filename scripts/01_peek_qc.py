# -*- coding: utf-8 -*-
"""
01_peek_qc.py —— 队列 fastq 质量概况（Python 版，替代 01_peek.R 更稳）
用法: py_env python 01_peek_qc.py <PRJNA号> [样本数上限] [每位置输出间隔bp]
输出: 每样本 reads/读长分布/平均Q/Q20比例 + 每位置平均Q曲线（文本+PNG）
目的: 确定 DADA2 truncLen 参数
"""
import gzip
import sys
import os

BASE = "<DATA_DIR>/sra"


def peek(fastq_gz, max_reads=50000):
    rlens = []
    qsum = []   # 每位置 Q 累计
    qcnt = []   # 每位置计数
    with gzip.open(fastq_gz, "rb") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip(b"\r\n")
            fh.readline()  # +
            qual = fh.readline().rstrip(b"\r\n")
            rlens.append(len(seq))
            for i, b in enumerate(qual):
                q = b - 33  # Phred+33
                if i >= len(qsum):
                    qsum.append(0)
                    qcnt.append(0)
                qsum[i] += q
                qcnt[i] += 1
            if len(rlens) >= max_reads:
                break
    return rlens, qsum, qcnt


def main():
    args = sys.argv[1:]
    cohort = args[0] if len(args) >= 1 else "PRJNA302078"
    n_max = int(args[1]) if len(args) >= 2 else 3
    step = int(args[2]) if len(args) >= 3 else 25
    fq_dir = os.path.join(BASE, cohort)
    fqs = sorted([f for f in os.listdir(fq_dir) if f.endswith(".fastq.gz")])
    print(f"共 {len(fqs)} 个 fastq，抽查前 {min(n_max, len(fqs))} 个")
    all_qpos = []
    for f in fqs[:min(n_max, len(fqs))]:
        fpath = os.path.join(fq_dir, f)
        rlens, qsum, qcnt = peek(fpath)
        n = len(rlens)
        mean_len = sum(rlens) / n
        qflat = []
        for s, c in zip(qsum, qcnt):
            qflat.extend([s / c] * c)
        mean_q = sum(qflat) / len(qflat) if qflat else 0
        lt20 = sum(1 for q in qflat if q < 20) / len(qflat) * 100 if qflat else 0
        print(f"===== {f} =====")
        print(f"  reads={n}  读长: mean={mean_len:.0f} median={sorted(rlens)[n//2]} min={min(rlens)} max={max(rlens)}")
        print(f"  平均Q={mean_q:.1f}  Q<20比例={lt20:.1f}%  最长读长位置数={len(qcnt)}")
        # 每 step bp 的平均 Q 曲线
        curve = []
        for i in range(0, len(qcnt), step):
            s = sum(qsum[i:i + step]) / max(1, sum(qcnt[i:i + step]))
            curve.append(s)
        all_qpos.append((f, curve))
        print(f"  每{step}bp平均Q: " + " ".join(f"{s:.0f}" for s in curve[:30]))
    # 画质量曲线 PNG
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 6))
        for f, curve in all_qpos:
            xs = [i * step for i in range(len(curve))]
            ax.plot(xs, curve, label=f.split(".")[0])
        ax.axhline(20, color="gray", ls="--", lw=0.8)
        ax.axhline(30, color="gray", ls="--", lw=0.8)
        ax.set_xlabel("Position (bp)"); ax.set_ylabel("Mean Phred Q")
        ax.set_ylim(0, 50); ax.legend()
        out = os.path.join("<REPO_DIR>/results",
                           f"{cohort}_quality_py.png")
        fig.savefig(out, dpi=110, bbox_inches="tight")
        print(f"质量图PNG已存: {out}")
    except Exception as e:
        print(f"(matplotlib 画图失败: {e})")


if __name__ == "__main__":
    main()
