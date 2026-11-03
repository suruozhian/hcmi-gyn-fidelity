#!/usr/bin/env python3
"""Download every public input needed to rebuild the analysis from scratch, with md5 verification.

  python3 fetch_public_data.py --data data            # all automatic downloads (GDC, figshare, GEO, ArrayExpress, HGNC, Nature ESM)
  python3 fetch_public_data.py --data data --check    # verify only (no download)

Sources and checksums: data_manifest.csv (non-GDC files) and inputs/manifests/gdc_download_list.csv (5,328 GDC open-access
files: TCGA 11 projects, TCGA pilot subsample, HCMI RNA, CPTAC; HCMI masked MAFs, ASCAT gene-level copy
number and SeSAMe methylation beta values; GDC Data Release 46.0). Rows with retrieval == "api:<name>" are small extracts from
public web APIs (cBioPortal study ucec_cptac_2020; GDC cases endpoint) written by public_api.py. Existing files with
the correct md5 are skipped, so the script can be re-run after an interruption (progress in <data>/fetch_progress.log).

Manual step (cannot be automated: the DepMap portal is behind a Cloudflare Turnstile browser check):
rows of data_manifest.csv with retrieval == "manual" (DepMap Public 26Q3 OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
and Model.csv; files of the NextGen Model Manuscript 2026 file set) must be downloaded in a browser from
https://depmap.org/portal/data_page/?tab=allData and saved under the listed local_path. This script verifies their md5.
"""
import argparse, hashlib, json, os, sys, time, zipfile, urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pandas as pd
PIPE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser()
ap.add_argument('--data',default=str(PIPE/'data'))
ap.add_argument('--check',action='store_true')
ap.add_argument('--threads',type=int,default=16)
ap.add_argument('--skip-gdc',action='store_true')
ap.add_argument('--gdc-sets',default='',help='comma-separated subset of the sets column of gdc_download_list.csv (default: all)')
a=ap.parse_args()
DATA=Path(a.data).resolve(); DATA.mkdir(parents=True,exist_ok=True)
ROOT=DATA.parent   # data_manifest local_path values start with data/
LOG=open(DATA/'fetch_progress.log','a')
def log(*x):
    s=time.strftime('%Y-%m-%d %H:%M:%S ')+' '.join(map(str,x)); print(s,flush=True); LOG.write(s+'\n'); LOG.flush()
def md5(p):
    h=hashlib.md5()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<22),b''): h.update(b)
    return h.hexdigest()
def get(url,dest,want,tries=5):
    dest=Path(dest); dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists() and md5(dest)==want: return 'ok (present)'
    if a.check: return 'MISSING' if not dest.exists() else 'MD5 MISMATCH'
    for i in range(tries):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'gyn-model-fidelity-fetch/1.0'})
            with urllib.request.urlopen(req,timeout=300) as r, open(str(dest)+'.part','wb') as f:
                while True:
                    b=r.read(1<<20)
                    if not b: break
                    f.write(b)
            os.replace(str(dest)+'.part',dest)
            if md5(dest)==want: return 'ok (downloaded)'
            log('md5 mismatch, retry',dest.name)
        except Exception as e:
            log('error',dest.name,repr(e)[:200]); time.sleep(5*(i+1))
    return 'FAILED'
status={}
DM=pd.read_csv(PIPE/'data_manifest.csv')
for r in DM.itertuples():
    dest=ROOT/r.local_path
    if r.retrieval=='auto':
        status[r.key]=get(r.url,dest,r.md5); log(r.key,status[r.key])
for r in DM[DM.retrieval.str.startswith('extract:')].itertuples():   # members of downloaded archives
    src=ROOT/DM.set_index('key').loc[r.retrieval.split(':',1)[1],'local_path']; dest=ROOT/r.local_path
    if not (dest.exists() and md5(dest)==r.md5) and src.exists() and not a.check:
        with zipfile.ZipFile(src) as z: dest.write_bytes(z.read(Path(r.local_path).name))
    status[r.key]='ok' if dest.exists() and md5(dest)==r.md5 else 'FAILED'; log(r.key,status[r.key])
for r in DM[DM.retrieval.str.startswith('api:')].itertuples():   # extracts from public web APIs (public_api.py)
    dest=ROOT/r.local_path; dest.parent.mkdir(parents=True,exist_ok=True)
    if not (dest.exists() and md5(dest)==r.md5) and not a.check:
        import public_api
        try: public_api.EXTRACT[r.retrieval.split(':',1)[1]](dest)
        except Exception as e: log('error',r.key,repr(e)[:200])
    status[r.key]=('ok (api extract verified)' if dest.exists() and md5(dest)==r.md5 else
                   'MD5 MISMATCH (api extract; the public source may have changed)' if dest.exists() else 'MISSING'); log(r.key,status[r.key])
manual_missing=[]
for r in DM[DM.retrieval=='manual'].itertuples():
    dest=ROOT/r.local_path
    ok=dest.exists() and md5(dest)==r.md5
    status[r.key]='ok (manual file verified)' if ok else ('MD5 MISMATCH (manual file)' if dest.exists() else 'MISSING (manual download required)')
    log(r.key,status[r.key])
    if not ok: manual_missing.append(f'  {r.local_path}\n      from {r.url}\n      expected md5 {r.md5}, {r.bytes:,} bytes')
if not a.skip_gdc:
    G=pd.read_csv(PIPE/'inputs/manifests/gdc_download_list.csv'); GD=DATA/'gdc'
    if a.gdc_sets: G=G[G.sets.str.split(';').apply(lambda x: bool(set(x)&set(a.gdc_sets.split(','))))]
    log('GDC files',len(G),f'{G.file_size.sum()/1e9:.1f} GB')
    def one(r): return r.file_id,get('https://api.gdc.cancer.gov/data/'+r.file_id,GD/f'{r.file_id}__{r.file_name}',r.md5sum)
    res={}; t0=time.time()
    with ThreadPoolExecutor(a.threads) as ex:
        for i,(fid,st) in enumerate(ex.map(one,list(G.itertuples())),1):
            res[fid]=st
            if i%250==0: log(f'GDC {i}/{len(G)} done ({time.time()-t0:.0f} s)')
    bad={k:v for k,v in res.items() if not v.startswith('ok')}
    status['gdc']=(f'ok ({len(res)}/{len(res)} files verified)' if not bad else f'FAILED ({len(bad)}/{len(res)} files not verified)'); log('GDC',status['gdc'],'failed:',list(bad)[:10])
json.dump(status,open(DATA/'fetch_status.json','w'),indent=1)
if manual_missing:
    log('MANUAL STEP REQUIRED (DepMap portal, Cloudflare browser check). Download in a browser and save as:\n'+'\n'.join(manual_missing))
failed=[k for k,v in status.items() if not str(v).startswith('ok')]
log('SUMMARY', 'all inputs present and verified' if not failed else f'not ready: {failed}')
sys.exit(0 if not failed else 1)
