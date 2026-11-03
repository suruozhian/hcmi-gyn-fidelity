#!/usr/bin/env python3
"""Build full-cohort manifest for 8 non-gyn TCGA projects (same rule as gyn reference:
open STAR - Counts Gene Expression Quantification, sample type Primary Tumor (code 01), one file per case),
record GDC data release, verify gyn manifest file IDs/md5 against the same release, and download with md5 check."""
import requests, json, time, hashlib, sys
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path, paths as cfg_paths, locate_tcga
D=cfg_path('tcga_nongyn_dir'); D.mkdir(parents=True,exist_ok=True)
FD=cfg_paths('tcga_file_dirs')[-1]; FD.mkdir(parents=True,exist_ok=True)   # downloads go to the last configured file directory
API='https://api.gdc.cancer.gov'
PROJ=['TCGA-BRCA','TCGA-COAD','TCGA-STAD','TCGA-CESC','TCGA-KIRC','TCGA-LUAD','TCGA-PAAD','TCGA-BLCA']
st=requests.get(API+'/status',timeout=60).json(); print('GDC status',st)
json.dump(st,open(D/'gdc_status_at_manifest.json','w'),indent=1)
def query(filters, fields, size=20000):
    r=requests.post(API+'/files',json=dict(filters=filters,fields=','.join(fields),format='JSON',size=size),timeout=300); r.raise_for_status()
    return r.json()['data']['hits']
F=['file_id','file_name','md5sum','file_size','cases.submitter_id','cases.project.project_id','cases.samples.submitter_id','cases.samples.sample_type','analysis.workflow_type','data_type','access']
filt={'op':'and','content':[
 {'op':'in','content':{'field':'cases.project.project_id','value':PROJ}},
 {'op':'=','content':{'field':'data_type','value':'Gene Expression Quantification'}},
 {'op':'=','content':{'field':'analysis.workflow_type','value':'STAR - Counts'}},
 {'op':'=','content':{'field':'access','value':'open'}},
 {'op':'=','content':{'field':'cases.samples.sample_type','value':'Primary Tumor'}}]}
hits=query(filt,F)
rows=[]
for h in hits:
    c=h['cases'][0]; s=c['samples'][0]
    rows.append(dict(file_id=h['file_id'],project=c['project']['project_id'],case=c['submitter_id'],sample=s['submitter_id'],
                     sample_type=s['sample_type'],file_name=h['file_name'],md5sum=h['md5sum'],file_size=h['file_size']))
A=pd.DataFrame(rows); A=A[A['sample'].str[13:15]=='01']
print('primary tumour STAR files',A.groupby('project').size().to_dict())
# one file per case: keep the file already used in the 40/project subsample if present; else lowest sample vial (01A<01B), then file_id
ng=pd.read_csv(cfg_path('tcga_nongyn_subsample_manifest')); prev=set(ng.file_id)
A['_prev']=~A.file_id.isin(prev)
A=A.sort_values(['case','_prev','sample','file_id']).drop_duplicates('case').drop(columns='_prev')
print('one per case',A.groupby('project').size().to_dict(),'total',len(A))
# local reuse
# local_path = file name <file_id>__<file_name>; resolved against tcga_file_dirs (files already present locally are reused)
A['local_path']=[f"{f}__{n}" for f,n in zip(A.file_id,A.file_name)]
A['gdc_release']=st['data_release']
A.to_csv(D/'nongyn_full_manifest.csv',index=False)
# verify gyn manifest under same release
G=pd.read_csv(cfg_path('tcga_gyn_full_manifest'))
gh=query({'op':'in','content':{'field':'file_id','value':G.file_id.tolist()}},['file_id','md5sum'])
gm={h['file_id']:h['md5sum'] for h in gh}
G['in_current_release']=G.file_id.isin(gm); G['md5_same']=[gm.get(f)==m for f,m in zip(G.file_id,G.md5sum)]
G[['file_id','project','case','in_current_release','md5_same']].to_csv(D/'gyn_manifest_release_check.csv',index=False)
print('gyn files present in current release',int(G.in_current_release.sum()),'/',len(G),'md5 same',int(G.md5_same.sum()))
# also check whether current release has gyn primary-tumour cases not in our gyn manifest (completeness)
gf=dict(filt); gf={'op':'and','content':[{'op':'in','content':{'field':'cases.project.project_id','value':['TCGA-OV','TCGA-UCEC','TCGA-UCS']}}]+filt['content'][1:]}
gh2=query(gf,['cases.submitter_id','cases.project.project_id','cases.samples.submitter_id'])
gc=pd.DataFrame([dict(case=h['cases'][0]['submitter_id'],project=h['cases'][0]['project']['project_id'],sample=h['cases'][0]['samples'][0]['submitter_id']) for h in gh2])
gc=gc[gc['sample'].str[13:15]=='01'].drop_duplicates('case')
print('current-release gyn primary cases',gc.groupby('project').size().to_dict(),'manifest',G.groupby('project').size().to_dict())
open(D/'release_check.txt','w').write(f"GDC {st['data_release']}\nnongyn one-per-case: {A.groupby('project').size().to_dict()}\n"
 f"gyn manifest files present in current release: {int(G.in_current_release.sum())}/{len(G)}; md5 identical: {int(G.md5_same.sum())}\n"
 f"current-release gyn primary-tumour cases: {gc.groupby('project').size().to_dict()}; manifest: {G.groupby('project').size().to_dict()}\n")
if '--manifest-only' in sys.argv: sys.exit()
need=A[~A.local_path.map(lambda p: locate_tcga(p).exists())]
print('to download',len(need))
def md5(p):
    h=hashlib.md5(); 
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()
def dl(r):
    for k in range(6):
        try:
            x=requests.get(f"{API}/data/{r.file_id}",timeout=300); x.raise_for_status()
            if hashlib.md5(x.content).hexdigest()!=r.md5sum: raise ValueError('md5 mismatch')
            (FD/Path(r.local_path).name).write_bytes(x.content); return r.file_id,'ok'
        except Exception as e:
            err=str(e); time.sleep(3*(k+1))
    return r.file_id,'fail:'+err
res=[]
with ThreadPoolExecutor(8) as ex:
    fut=[ex.submit(dl,r) for r in need.itertuples()]
    for i,f in enumerate(as_completed(fut),1):
        res.append(f.result())
        if i%100==0: print('progress',i,'/',len(need),'fail',sum(1 for _,s in res if s!='ok'),flush=True)
R=pd.DataFrame(res,columns=['file_id','status']); R.to_csv(D/'download_status.csv',index=False)
# md5 check of reused local files
A['md5_observed']=[md5(locate_tcga(p)) if locate_tcga(p).exists() else None for p in A.local_path]
A['md5_ok']=A.md5_observed==A.md5sum
A.to_csv(D/'nongyn_full_manifest.csv',index=False)
print('DONE ok',int(A.md5_ok.sum()),'/',len(A),'fail',(R.status!='ok').sum())
