# -*- coding: utf-8 -*-
"""
下载 SRA 开放队列 16S FASTQ（走 ENA，中国网络可达）。
用法: python download_sra_cohorts.py PRJNA302078 PRJNA592384 ... [--threads N]
      --threads N: 并发下载数（默认 4；ENA 单连接限速，并发可显著提速）
输出: <DATA_DIR>/sra/<PRJNA号>/<run>_1.fastq.gz / <run>_2.fastq.gz
      （原始数据不参与 WPS 云同步，放在工作空间之外）
特性: 断点续传（已存在文件跳过）、失败重试 3 次、运行日志记录、多线程并发。
"""
import os
import sys
import time
import threading
import urllib.request
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = r"<DATA_DIR>/sra"
ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = "run_accession,library_strategy,library_source,fastq_ftp"
LOG = os.path.join(BASE, "download_log.txt")
_log_lock = threading.Lock()

def log(msg):
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    with _log_lock:
        print(line)
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(line + "\n")

def fetch_runs(project):
    url = (f"{ENA_API}?accession={project}&result=read_run"
           f"&fields={urllib.parse.quote(FIELDS)}&format=tsv")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        lines = r.read().decode("utf-8").strip().splitlines()
    header = lines[0].split("\t")
    idx_run, idx_strat, idx_src, idx_fq = header.index("run_accession"), header.index("library_strategy"), header.index("library_source"), header.index("fastq_ftp")
    out = []
    for ln in lines[1:]:
        cells = ln.split("\t")
        if len(cells) <= max(idx_run, idx_fq):
            continue
        strat = cells[idx_strat] if len(cells) > idx_strat else ""
        fq = cells[idx_fq] if len(cells) > idx_fq and cells[idx_fq] else ""
        if fq and strat == "AMPLICON":  # 只取 16S 扩增子
            out.append((cells[idx_run], fq))
    return out

def download(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        log(f"  跳过已存在: {os.path.basename(dest)}")
        return True
    tmp = dest + ".part"
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            os.replace(tmp, dest)
            log(f"  OK {os.path.basename(dest)} ({os.path.getsize(dest)//1024//1024} MB)")
            return True
        except Exception as e:
            log(f"  重试{attempt+1}/3 失败: {url.split('/')[-1]} ({e})")
            time.sleep(5)
    return False

def build_tasks(proj_dir, runs):
    """把 run 列表展开为 (url, dest) 下载任务；双端只取 _1/_2。"""
    tasks = []
    for run, fq in runs:
        files = [u for u in fq.split(";") if u]
        # 双端 run：ENA 常同时提供 <run>.fastq.gz(全集) + <run>_1 + <run>_2，
        # 只取 _1/_2，避免无后缀全集被误存为 R1（导致双端错乱）。
        paired_urls = [u for u in files
                       if u.split("/")[-1].endswith("_1.fastq.gz")
                       or u.split("/")[-1].endswith("_2.fastq.gz")]
        sel = paired_urls if paired_urls else files  # 真单端(如454)取无后缀文件
        for u in sel:
            https = u.replace("ftp.sra.ebi.ac.uk", "https://ftp.sra.ebi.ac.uk")
            base = u.split("/")[-1]
            if base.endswith("_1.fastq.gz"):
                dest = os.path.join(proj_dir, f"{run}_1.fastq.gz")
            elif base.endswith("_2.fastq.gz"):
                dest = os.path.join(proj_dir, f"{run}_2.fastq.gz")
            else:
                dest = os.path.join(proj_dir, f"{run}.fastq.gz")
            tasks.append((https, dest))
    return tasks

def main(projects, threads=4):
    os.makedirs(BASE, exist_ok=True)
    for proj in projects:
        proj_dir = os.path.join(BASE, proj)
        os.makedirs(proj_dir, exist_ok=True)
        log(f"=== {proj}: 查询 ENA ===")
        try:
            runs = fetch_runs(proj)
        except Exception as e:
            log(f"查询失败: {e}")
            continue
        tasks = build_tasks(proj_dir, runs)
        log(f"{proj}: 共 {len(runs)} 个 16S 扩增子 run, {len(tasks)} 个文件, 并发 {threads}")
        with ThreadPoolExecutor(max_workers=threads) as ex:
            futs = [ex.submit(download, u, d) for u, d in tasks]
            for fut in as_completed(futs):
                fut.result()  # 内部已 log，异常已在 download 内处理
    log("=== 全部完成 ===")

if __name__ == "__main__":
    argv = sys.argv[1:]
    threads = 4
    if "--threads" in argv:
        i = argv.index("--threads")
        threads = int(argv[i + 1])
        del argv[i:i + 2]
    projs = argv or ["PRJNA302078", "PRJNA592384"]
    main(projs, threads)
