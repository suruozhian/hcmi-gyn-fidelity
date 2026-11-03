#!/usr/bin/env python3
"""Build the HCMI RNA matrix, the GENCODE v36 gene annotation and the relative gene lengths from the downloaded GDC
STAR gene-count files (open access; fetch_public_data.py).

Inputs : inputs/manifests/hcmi_rna_manifest.csv (55 HCMI RNA files: 33 models + 22 tumours; file_id, md5, column name)
         GDC files <file_id>__<file_name> in the configured tcga_file_dirs (or --gdc-dir)
Outputs: <hcmi_processed_dir>/hcmi_rna_tpm.parquet  (STAR tpm_unstranded; Ensembl gene IDs with version; one column per file)
         <hcmi_processed_dir>/hcmi_rna_samples.csv (the manifest, for reference)
         <gene_annotation>       gene_id -> gene_name, gene_type (GENCODE v36, as in every GDC STAR file)
         <relative_gene_lengths> per-gene median of (unstranded count / TPM), each file scaled by its median, over the first
                                 20 HCMI STAR files in file-name order, counts >= 20 (used to convert counts-only resources to TPM)
This reproduces the original processing exactly (same files, same order)."""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path, paths as cfg_paths
PIPE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(); ap.add_argument('--gdc-dir',default=None); a=ap.parse_args()
dirs=[Path(a.gdc_dir)] if a.gdc_dir else cfg_paths('tcga_file_dirs')
def find(name):
    for d in dirs:
        if (d/name).exists(): return d/name
    raise FileNotFoundError(name)
m=pd.read_csv(PIPE/'inputs/manifests/hcmi_rna_manifest.csv')
m['fname']=m.file_id+'__'+m.file_name
mats={}; counts={}
for x in m.itertuples():
    t=pd.read_csv(find(x.fname),sep='\t',comment='#'); t=t[t.gene_id.str.startswith('ENSG')]
    mats[x.col]=t.set_index('gene_id').tpm_unstranded
genes=t.set_index('gene_id')[['gene_name','gene_type']]
OUT=cfg_path('hcmi_processed_dir'); OUT.mkdir(parents=True,exist_ok=True)
pd.DataFrame(mats).to_parquet(OUT/'hcmi_rna_tpm.parquet'); m.drop(columns='fname').to_csv(OUT/'hcmi_rna_samples.csv',index=False)
ga=cfg_path('gene_annotation'); ga.parent.mkdir(parents=True,exist_ok=True); genes.to_csv(ga)
R=[]
for f in sorted(m.fname)[:20]:
    t=pd.read_csv(find(f),sep='\t',comment='#'); t=t[t.gene_id.str.startswith('ENSG')].set_index('gene_id')
    r=(t.unstranded/t.tpm_unstranded).replace([np.inf,-np.inf],np.nan); r=r[(t.unstranded>=20)]
    R.append(r/r.median())
L=pd.concat(R,axis=1).median(axis=1); L.index=L.index.str.split('.').str[0]
rl=cfg_path('relative_gene_lengths'); rl.parent.mkdir(parents=True,exist_ok=True); L.to_csv(rl)
print('HCMI RNA',len(mats),'files; genes',len(genes),'; relative lengths',int(L.notna().sum()))
