#!/usr/bin/env python3
"""Ovarian histotype reference self-checks (Results 3.4), recomputed with the pipeline gene key and gene set.
Same definitions as the published histotype reference analysis: top-300 ANOVA genes on GSE295399 malignant samples (genes shared with
GSE226870, GSE121103 and the harmonised model matrix; mean log2(TPM+1) > 1); Spearman nearest histotype centroid.
  GSE295399 leave-one-out (centroids rebuilt without the held-out sample);
  GSE226870 external (CCOC -> CCC, ENOC -> EC);  GSE121103 external.
Usage: python3 03_histotype_reference_checks.py [--gene-id entrez|symbol] [--out results/histotype_reference_checks_<gene-id>.json]
"""
import sys, argparse, json, gzip
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import f_oneway
PIPE=Path(__file__).resolve().parent; sys.path.insert(0,str(PIPE))
from paths import path as cfg_path
import genemap
from fidelity_core import aligned, assert_same_order
ap=argparse.ArgumentParser(); ap.add_argument('--gene-id',choices=['entrez','symbol'],default='entrez'); ap.add_argument('--out',default=None); a=ap.parse_args(); GID=a.gene_id
ga=pd.read_csv(cfg_path('gene_annotation'),index_col=0); ga['ens']=ga.index.str.split('.').str[0]; pc=ga[ga.gene_type=='protein_coding']
L=pd.read_csv(cfg_path('relative_gene_lengths'),index_col=0).iloc[:,0]; L.index=L.index.str.split('.').str[0]
def c2t(C):
    l=L.reindex(C.index); C=C[l.notna()]; l=l[l.notna()]; rpk=C.div(l,axis=0); return np.log2(rpk/rpk.sum()*1e6+1)
def e2k(M):
    if GID=='entrez': return genemap.to_entrez(M,genemap.ens2entrez().entrez)
    m=pc.drop_duplicates('ens').set_index('ens').gene_name; M=M[M.index.isin(m.index)]; M.index=m.loc[M.index].values; return M.groupby(level=0).max()
def series_chars(f):
    d={}
    for line in gzip.open(f,'rt'):
        if line.startswith('!Sample_title'): d['title']=[x.strip('"\n') for x in line.split('\t')[1:]]
        if line.startswith('!Sample_characteristics_ch1'):
            v=[x.strip('"\n') for x in line.split('\t')[1:]]; k=v[0].split(':')[0]; d[k]=[x.split(': ',1)[-1] for x in v]
    return pd.DataFrame(d)
m1=series_chars(str(cfg_path('gse295399_series_matrix'))); m1['col']=m1.title.str.replace('Tumor_','')
C1=pd.read_csv(cfg_path('gse295399_counts'),sep='\t',index_col=0); C1.index=C1.index.str.split('.').str[0]
R1=e2k(c2t(C1)); m1=m1[m1.col.isin(R1.columns)].set_index('col')
keep=m1[m1.histotype.isin(['HGSC','LGSC','CCC','EC','MC'])]; R1=R1[keep.index]; y=keep.histotype
if GID=='entrez':
    R2=pd.read_csv(cfg_path('gse226870_matrix'),sep='\t').drop(columns='SYMBOL'); R2['ENSEMBL']=R2.ENSEMBL.astype(str).str.split('.').str[0]; R2=e2k(R2.groupby('ENSEMBL').max())
else:
    R2=pd.read_csv(cfg_path('gse226870_matrix'),sep='\t').drop(columns='ENSEMBL').groupby('SYMBOL').max()
m2=series_chars(str(cfg_path('gse226870_series_matrix'))).set_index('title')
cmap=pd.read_csv(cfg_path('gse121103_colmap')); C3=pd.read_csv(cfg_path('gse121103_counts'),sep='\t').set_index('Geneid')
C3.index=C3.index.str.split('.').str[0]; C3=C3[list(cmap.count_col)].rename(columns=dict(zip(cmap.count_col,cmap.title.str.replace('_RNA-seq',''))))
R3=e2k(c2t(C3)); y3=cmap.set_index(cmap.title.str.replace('_RNA-seq',''))['histotype_short'].reindex(R3.columns)
HARM=cfg_path('reference_dir')/'expr_entrez_log2tpm_models.parquet' if GID=='entrez' else cfg_path('harmonised_matrix')
xg=pd.read_parquet(HARM).index
genes=sorted(set(R1.index)&set(xg)&set(R2.index)&set(R3.index))
R1g=R1.loc[genes]; R1g=R1g[R1g.mean(axis=1)>1]
F=pd.Series(f_oneway(*[R1g.loc[:,y==c].values.T for c in sorted(y.unique())]).statistic,index=R1g.index).sort_values(ascending=False)
G=F.index[:300].tolist()
def cent(M,lab): return pd.concat({c:aligned(G,M)[lab.index[lab==c]].mean(axis=1) for c in sorted(lab.unique())},axis=1)
def nearest(S,Cn):
    S=aligned(G,S); assert_same_order(S,Cn); return pd.DataFrame({c:S.corrwith(Cn[c],method='spearman') for c in Cn.columns}).idxmax(axis=1)
pred=pd.Series({s:nearest(R1[[s]],cent(R1,y.drop(s))).iloc[0] for s in y.index})
Cn=cent(R1,y)
pr2=nearest(R2,Cn); t2=m2.histology.reindex(pr2.index).map({'CCOC':'CCC','ENOC':'EC'})
pr3=nearest(R3,Cn)
out=dict(gene_id=GID,genes_top300_from=len(genes),GSE295399_LOO=f'{int((pred==y).sum())}/{len(y)}',
         GSE226870_external=f'{int((pr2==t2).sum())}/{int(t2.notna().sum())}',GSE121103_external=f'{int((pr3==y3).sum())}/{int(y3.notna().sum())}',
         GSE226870_crosstab=pd.crosstab(t2,pr2).to_dict(),GSE121103_crosstab=pd.crosstab(y3,pr3).to_dict())
print(json.dumps({k:v for k,v in out.items() if 'crosstab' not in k},indent=1))
o=Path(a.out) if a.out else PIPE/'results'/f'histotype_reference_checks_{GID}.json'; json.dump(out,open(o,'w'),indent=1,default=str); print('wrote',o)
