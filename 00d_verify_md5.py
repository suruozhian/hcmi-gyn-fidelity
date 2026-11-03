#!/usr/bin/env python3
import hashlib, pandas as pd
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path, paths as cfg_paths, locate_tcga
D=cfg_path('tcga_nongyn_dir'); A=pd.read_csv(D/'nongyn_full_manifest.csv')
def md5(p):
    h=hashlib.md5()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
A['md5_observed']=[md5(locate_tcga(p)) if locate_tcga(p).exists() else None for p in A.local_path]; A['md5_ok']=A.md5_observed==A.md5sum
A.to_csv(D/'nongyn_full_manifest.csv',index=False); print('md5 ok',int(A.md5_ok.sum()),'/',len(A)); print(A.groupby('project').md5_ok.agg(['sum','count']))
