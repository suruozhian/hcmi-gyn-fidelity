#!/usr/bin/env python3
"""Build the Entrez-keyed harmonised model matrix (all non-TCGA samples of the cohort annotation: HCMI models and
tumours, CPTAC tumours, DepMap 2D lines [24Q4 file], NextGen organoids [NextGen manuscript file set], Sanger organoids,
GSE307568, E-MTAB-10664, GSE208216) from each resource's native identifiers:
  Ensembl (HCMI, CPTAC, GSE307568, GSE208216)  -> Entrez via genemap rules E1/E2
  Entrez in header (DepMap, NextGen)           -> taken directly
  symbols only (Sanger organoids, E-MTAB-10664) -> Entrez via genemap rules S1/S2
Expression processing per resource is unchanged from the symbol-keyed matrix (log2(TPM+1); counts -> TPM with the
same relative gene lengths). Genes = Entrez IDs present in every resource except NextGen (NextGen reindexed; missing =
NaN), as for the symbol-keyed matrix.
Output: <reference_dir>/expr_entrez_log2tpm_models.parquet + gene_map_summary.json"""
import sys, json, gzip, re
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path, has as cfg_has
import genemap
OUT=cfg_path('reference_dir')
HP=cfg_path('hcmi_processed_dir')   # built by 00e_build_hcmi_rna_and_gene_annotation.py
A=pd.read_csv(cfg_path('sample_annotations'))
A=A[A.dataset!='TCGA']
need={d:list(g['sample']) for d,g in A.groupby('dataset')}
em=genemap.ens2entrez(); pc=genemap.gencode_pc()
L=pd.read_csv(cfg_path('relative_gene_lengths'),index_col=0).iloc[:,0]; L.index=L.index.str.split('.').str[0]
def counts_to_log2tpm(C):
    l=L.reindex(C.index); C=C[l.notna()]; l=l[l.notna()]; rpk=C.div(l,axis=0); return np.log2(rpk/rpk.sum()*1e6+1)
def ens_entrez(M): return genemap.to_entrez(M,em.entrez)
smap_log={}
def sym_entrez(M,name):
    sm=genemap.symbol_map(M.index); smap_log[name]=sm.rule.value_counts().to_dict(); return genemap.to_entrez(M,sm.entrez)
mats={}
# HCMI (GDC STAR TPM; Ensembl)
H=pd.read_parquet(HP/'hcmi_rna_tpm.parquet'); H.index=H.index.str.split('.').str[0]
H=H[H.index.isin(pc.index)]; H=H[~H.index.duplicated()]
mats['HCMI']=ens_entrez(np.log2(H[need['HCMI']].astype(float)+1))
# CPTAC (GDC STAR TPM; Ensembl)
cm=pd.read_csv(cfg_path('cptac_manifest')); C={}
for x in cm.itertuples():
    t=pd.read_csv(cfg_path('cptac_files_dir')/f'{x.file_id}__{x.file_name}',sep='\t',comment='#',usecols=['gene_id','tpm_unstranded'])
    t=t[t.gene_id.str.startswith('ENSG')]; t['ens']=t.gene_id.str.split('.').str[0]; t=t[t.ens.isin(pc.index)].drop_duplicates('ens').set_index('ens')
    C['CPTAC|'+x.case]=t.tpm_unstranded
C=pd.DataFrame(C); mats['CPTAC']=ens_entrez(np.log2(C[need['CPTAC']].astype(float)+1))
# DepMap 2D (24Q4 file; Entrez in header)
E=pd.read_csv(cfg_path('depmap_expr_24q4'),index_col=0)
ids=[s.split('|')[1] for s in need['DepMap_24Q4']]
Eg=E.loc[ids].T; Eg.index=genemap.depmap_entrez(Eg.index); Eg=Eg[Eg.index.notna()].groupby(level=0).max()
Eg.columns=['DepMap|'+c for c in Eg.columns]; mats['DepMap_24Q4']=Eg
# NextGen (manuscript file set; Entrez in header)
NG=pd.read_csv(cfg_path('nextgen_dir')/'next_gen_expression.csv',index_col=0)
ids=[s.split('|')[1] for s in need['DepMap_NextGen']]
Ng=NG.loc[ids].T; Ng.index=genemap.depmap_entrez(Ng.index); Ng=Ng[Ng.index.notna()].groupby(level=0).max(); Ng.columns=['NextGen|'+c for c in Ng.columns]
# Sanger organoids (symbols)
S=pd.read_csv(cfg_path('sanger_organoid_rna'),index_col=0); S.columns=S.columns.str.replace('_exp$','',regex=True)
So=S.loc[[s.split('|')[1] for s in need['Sanger_organoids']]].T; So.columns=['Sanger|'+c for c in So.columns]
mats['Sanger_organoids']=sym_entrez(So,'Sanger_organoids')
# GSE307568 (raw counts; Ensembl)
G=pd.read_csv(cfg_path('gse307568_counts'),sep='\t',index_col=0); G.index=G.index.str.split('.').str[0]
libs=[s.split('|')[1] for s in need['GSE307568']]
Gt=counts_to_log2tpm(G[libs]); Gt=Gt[Gt.index.isin(pc.index)]; Gt.columns=['GSE307568|'+c for c in Gt.columns]; mats['GSE307568']=ens_entrez(Gt)
# E-MTAB-10664 (normalised log2 counts; symbols) - same TPM-like conversion as before, then symbols -> Entrez
B=pd.read_csv(cfg_path('emtab10664_counts'),sep='\t',index_col=0)
lin=(2**B-1).clip(lower=0)
gp=pc.reset_index().rename(columns={'index':'ens'})
sl=gp.assign(L=lambda d: L.reindex(d.ens).values).dropna(subset=['L']).groupby('gene_name').L.median()
lin=lin[lin.index.isin(sl.index)]; tpm=lin.div(sl.loc[lin.index],axis=0); tpm=tpm/tpm.sum()*1e6
Bt=np.log2(tpm+1); Bt=Bt.groupby(level=0).max(); Bt.columns=['Berg2021|'+c for c in Bt.columns]
mats['Berg2021_EMTAB10664']=sym_entrez(Bt[need['Berg2021_EMTAB10664']],'E-MTAB-10664')
# GSE208216 (raw counts; Ensembl)
C8=pd.read_csv(cfg_path('gse208216_counts'),sep='\t',index_col=0); C8.index=C8.index.str.split('.').str[0]
C8=C8[[s.split('|')[1] for s in need['GSE208216']]]
G8=counts_to_log2tpm(C8); G8=G8[G8.index.isin(pc.index)]; G8.columns=['GSE208216|'+c for c in G8.columns]; mats['GSE208216']=ens_entrez(G8)
genes=sorted(set.intersection(*[set(v.index) for v in mats.values()]),key=int)
X=pd.concat([v.reindex(genes) for v in mats.values()]+[Ng.reindex(genes)],axis=1).astype('float32')
X=X[[s for s in A['sample'] if s in X.columns]]
assert set(X.columns)==set(A['sample']), set(A['sample'])^set(X.columns)
X.to_parquet(OUT/'expr_entrez_log2tpm_models.parquet')
summ=dict(genes_common=len(genes),samples=X.shape[1],per_resource_genes={k:int(v.shape[0]) for k,v in mats.items()},nextgen_genes=int(Ng.shape[0]),
          nextgen_missing_in_common=int(X[[c for c in X.columns if c.startswith('NextGen|')]].isna().all(axis=1).sum()),
          ensembl_rules=em.rule.value_counts().to_dict(),symbol_rules=smap_log)
# QC (optional; only where the earlier symbol-keyed matrix is available): identical values for genes with a 1:1 GENCODE name
if cfg_has('harmonised_matrix') and cfg_path('harmonised_matrix').exists():
  Xs=pd.read_parquet(cfg_path('harmonised_matrix'),columns=list(X.columns))
  e1=em[em.rule.isin(['E1','E2'])]; one=e1.groupby('entrez').filter(lambda d: len(d)==1)
  one=one[~one.gene_name.isin(pc.gene_name[pc.gene_name.duplicated()])]
  s2e=one.reset_index().set_index('gene_name').entrez
  qc={}
  for ds,cols in need.items():
      gs=[g for g in Xs.index if g in s2e.index and s2e[g] in X.index]
      d=np.nanmax(np.abs(Xs.loc[gs,cols].values-X.loc[s2e[gs].values,cols].values)); qc[ds]=dict(genes=len(gs),max_abs_diff=float(d))
  summ['qc_vs_symbol_matrix_1to1_genes']=qc
json.dump(summ,open(OUT/'gene_map_summary.json','w'),indent=1)
print(json.dumps(summ,indent=1))
