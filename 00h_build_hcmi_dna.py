#!/usr/bin/env python3
"""Build the HCMI DNA inputs of the identity check (Methods 2.3) from GDC open-access files (fetch_public_data.py).

Inputs : inputs/manifests/hcmi_cn_manifest.csv   (56 ASCAT gene-level copy-number files, open access, GDC Data Release 46.0)
         inputs/manifests/hcmi_maf_manifest.csv  (60 masked somatic mutation MAFs, open access)
         GDC files <file_id>__<file_name> in tcga_file_dirs (or --gdc-dir)
Outputs: <hcmi_processed_dir>/hcmi_cn_gene.parquet  (copy_number by Ensembl gene ID; one column per file, '<sample>|<file_id[:8]>')
         <hcmi_processed_dir>/hcmi_cn_samples.csv   (manifest with column key 'col')
         <hcmi_processed_dir>/hcmi_maf_combined.csv.gz (selected MAF columns of all files, with file_id and HCMI sample)
         <hcmi_processed_dir>/hcmi_maf_samples.csv  (sample, file_id, case, role, lineage, histology)
Only open-access (masked) files are used; no controlled-access data. Same processing as the original analysis."""
import argparse, sys
from pathlib import Path
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path, paths as cfg_paths
PIPE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(); ap.add_argument('--gdc-dir',default=None); a=ap.parse_args()
dirs=[Path(a.gdc_dir)] if a.gdc_dir else cfg_paths('tcga_file_dirs')
def find(name):
    for d in dirs:
        if (d/name).exists(): return d/name
    raise FileNotFoundError(name)
OUT=cfg_path('hcmi_processed_dir'); OUT.mkdir(parents=True,exist_ok=True)
# copy number
cn=pd.read_csv(PIPE/'inputs/manifests/hcmi_cn_manifest.csv'); C={}
for x in cn.itertuples():
    t=pd.read_csv(find(f'{x.file_id}__{x.file_name}'),sep='\t'); C[x.col]=t.set_index('gene_id').copy_number
pd.DataFrame(C).to_parquet(OUT/'hcmi_cn_gene.parquet'); cn.to_csv(OUT/'hcmi_cn_samples.csv',index=False)
# masked somatic mutations
mf=pd.read_csv(PIPE/'inputs/manifests/hcmi_maf_manifest.csv')
cols=['Hugo_Symbol','Chromosome','Start_Position','Variant_Classification','Variant_Type','HGVSp_Short','Tumor_Sample_Barcode','t_depth','t_alt_count','IMPACT','Protein_position']
L=[]
for x in mf.itertuples():
    t=pd.read_csv(find(f'{x.file_id}__{x.file_name}'),sep='\t',comment='#',usecols=cols,low_memory=False); t['file_id']=x.file_id
    t['sample']=t.Tumor_Sample_Barcode.str.extract(r'^(HCM-[A-Z]+-\d+-C\d+-\w+?)-\d\d[DR]')[0]
    L.append(t)
M=pd.concat(L); M.to_csv(OUT/'hcmi_maf_combined.csv.gz',index=False)
s=M[['sample','file_id']].drop_duplicates().merge(mf[['sample','case','role','lineage','histology']].drop_duplicates(),on='sample',how='left')
assert s.role.notna().all() and len(s)==len(mf)
s.to_csv(OUT/'hcmi_maf_samples.csv',index=False)
print('HCMI copy-number files',len(cn),'(genes',len(C[next(iter(C))]),'); MAF files',len(mf),'variants',len(M))
