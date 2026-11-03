#!/usr/bin/env python3
"""Build the stromal/immune purity-filter gene list (genes excluded before feature-gene selection; Methods 2.2).

Genes expressed in TCGA (mean log2(TPM+1) > 1) whose expression ranks correlate (Spearman rho > 0.4) with a simple
stromal score (COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRB, FAP, THY1, COL6A3, SPARC) or immune score (PTPRC, CD3E, CD2,
CD14, CD68, LYZ, CD74, HLA-DRA, MS4A1, CD53) across the original TCGA subsample (777 primary tumours:
inputs/manifests/tcga_gyn_subsample_manifest.csv; OV, UCEC, UCS and 40/project of 8 other projects; GDC STAR files).
Scores = mean z-score of log2(TPM+1) of the marker genes. Output: <purity_filtered_genes> (one GENCODE gene name per row),
identical to the list used in the analysis (checked against inputs/intermediate/rna_purity_filtered_genes.reference.csv)."""
import argparse, hashlib, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path, paths as cfg_paths
PIPE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(); ap.add_argument('--gdc-dir',default=None); a=ap.parse_args()
dirs=[Path(a.gdc_dir)] if a.gdc_dir else cfg_paths('tcga_file_dirs')
def find(name):
    for d in dirs:
        if (d/name).exists(): return d/name
    raise FileNotFoundError(name)
m=pd.read_csv(PIPE/'inputs/manifests/tcga_gyn_subsample_manifest.csv')
T={}
for x in m.itertuples():
    t=pd.read_csv(find(f'{x.file_id}__{x.file_name}'),sep='\t',comment='#',usecols=['gene_id','tpm_unstranded'])
    T[x.case]=t[t.gene_id.str.startswith('ENSG')].set_index('gene_id').tpm_unstranded.astype('float32')
T=pd.DataFrame(T)
H=pd.read_parquet(cfg_path('hcmi_processed_dir')/'hcmi_rna_tpm.parquet')
ga=pd.read_csv(cfg_path('gene_annotation'),index_col=0)
pc=ga.index[ga.gene_type=='protein_coding']
common=pc.intersection(H.index).intersection(T.index)
LT=np.log2(T.loc[common].astype(float)+1)
sym=ga.gene_name
def score(L,genes):
    ids=[g for g in sym.index[sym.isin(genes)] if g in L.index]
    z=L.loc[ids]; z=(z.sub(z.mean(1),axis=0)).div(z.std(1)+1e-9,axis=0); return z.mean()
stromal=['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRB','FAP','THY1','COL6A3','SPARC']
immune=['PTPRC','CD3E','CD2','CD14','CD68','LYZ','CD74','HLA-DRA','MS4A1','CD53']
sS=score(LT,stromal); sI=score(LT,immune)
expr=LT.index[(LT.mean(1)>1)]
R=LT.loc[expr].T.apply(lambda c: c.rank())
rS=R.corrwith(pd.Series(rankdata(sS),index=sS.index)); rI=R.corrwith(pd.Series(rankdata(sI),index=sI.index))
purity_genes=set(expr[(rS.values>0.4)|(rI.values>0.4)])
out=cfg_path('purity_filtered_genes'); out.parent.mkdir(parents=True,exist_ok=True)
pd.Series(sorted(sym[list(purity_genes)])).to_csv(out,index=False)
ref=PIPE/'inputs/intermediate/rna_purity_filtered_genes.reference.csv'
same=hashlib.md5(out.read_bytes()).hexdigest()==hashlib.md5(ref.read_bytes()).hexdigest()
print('TCGA subsample',T.shape[1],'tumours; expressed genes',len(expr),'; purity-associated genes',len(purity_genes),'; identical to shipped reference list:',same)
sys.exit(0 if same else 1)
