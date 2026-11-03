#!/usr/bin/env python3
"""Helper: download remaining non-gyn files in reverse manifest order (skip existing), md5-verified, atomic rename."""
import requests, hashlib, time, os
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path, paths as cfg_paths, locate_tcga
D=cfg_path('tcga_nongyn_dir'); A=pd.read_csv(D/'nongyn_full_manifest.csv')
def dl(r):
    p=locate_tcga(r.local_path)
    if p.exists(): return 'skip'
    for k in range(6):
        try:
            x=requests.get(f'https://api.gdc.cancer.gov/data/{r.file_id}',timeout=300); x.raise_for_status()
            if hashlib.md5(x.content).hexdigest()!=r.md5sum: raise ValueError('md5')
            if p.exists(): return 'skip'
            tmp=p.with_suffix('.part_fwd'); tmp.write_bytes(x.content); os.replace(tmp,p); return 'ok'
        except Exception: time.sleep(3*(k+1))
    return 'fail'
with ThreadPoolExecutor(8) as ex: res=list(ex.map(dl,A.itertuples()))
print(pd.Series(res).value_counts().to_dict())
