#!/usr/bin/env python3
"""Build the gene-symbol-keyed model matrix used by the gene-symbol sensitivity analysis (`--gene-id symbol`) and the
symbol-keyed identity check (configuration key `harmonised_matrix`).

Same resource-native processing as the original (earlier, symbol-keyed) analysis:
- GDC STAR files (HCMI, CPTAC): log2(TPM+1) of GENCODE v36 protein-coding genes, keyed by gene name (max over duplicates);
- DepMap 24Q4 (2D lines) and NextGen organoids (DepMap portal file): log2(TPM+1) as provided, header symbol (max over duplicates);
- Sanger organoids: log2(TPM+1) as provided; GSE307568 and GSE208216 raw counts -> TPM with relative gene lengths;
  E-MTAB-10664 normalised counts -> TPM-like with relative gene lengths (as in 02_build_model_matrix.py);
- gene set: protein-coding symbols present in every resource and in the TCGA reference (sorted).
Columns: all non-TCGA samples of the sample annotation, in its order. Output: <harmonised_matrix> (parquet)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path
A=pd.read_csv(cfg_path('sample_annotations')); A=A[A.dataset!='TCGA']
cols=lambda ds:[s for s in A['sample'] if A.set_index('sample').loc[s,'dataset']==ds]
AS=A.set_index('sample')
ga=pd.read_csv(cfg_path('gene_annotation'),index_col=0); ga['ens']=ga.index.str.split('.').str[0]
pc=ga[ga.gene_type=='protein_coding']; sym_of=pc.gene_name; ens2sym=pc.drop_duplicates('ens').set_index('ens').gene_name
L=pd.read_csv(cfg_path('relative_gene_lengths'),index_col=0).iloc[:,0]; L.index=L.index.str.split('.').str[0]
def counts_to_log2tpm(C):
    l=L.reindex(C.index); C=C[l.notna()]; l=l[l.notna()]; rpk=C.div(l,axis=0); return np.log2(rpk/rpk.sum()*1e6+1)
def ens_to_symbol(M):
    M=M[M.index.isin(ens2sym.index)]; M.index=ens2sym.loc[M.index].values; return M.groupby(level=0).max()
def star_symbol(T):
    T=np.log2(T.loc[T.index.isin(pc.index)].astype(float)+1); T.index=sym_of.loc[T.index].values; return T.groupby(level=0).max()
mats={}
H=pd.read_parquet(cfg_path('hcmi_processed_dir')/'hcmi_rna_tpm.parquet'); mats['HCMI']=star_symbol(H)
cm=pd.read_csv(cfg_path('cptac_manifest')); C={}
for x in cm.itertuples():
    t=pd.read_csv(cfg_path('cptac_files_dir')/f'{x.file_id}__{x.file_name}',sep='\t',comment='#',usecols=['gene_id','tpm_unstranded'])
    C['CPTAC|'+x.case]=t[t.gene_id.str.startswith('ENSG')].set_index('gene_id').tpm_unstranded
mats['CPTAC']=star_symbol(pd.DataFrame(C))
E=pd.read_csv(cfg_path('depmap_expr_24q4'),index_col=0); E.columns=E.columns.str.replace(r' \(\d+\)$','',regex=True)
dm=[s.split('|')[1] for s in cols('DepMap_24Q4')]
Eg=E.loc[E.index.intersection(dm)].T.groupby(level=0).max(); Eg.columns=['DepMap|'+c for c in Eg.columns]; mats['DepMap']=Eg
S=pd.read_csv(cfg_path('sanger_organoid_rna'),index_col=0); S.columns=S.columns.str.replace('_exp$','',regex=True)
so=[s.split('|')[1] for s in cols('Sanger_organoids')]
So=S.loc[S.index.intersection(so)].T; So.columns=['Sanger|'+c for c in So.columns]; mats['Sanger']=So
G=pd.read_csv(cfg_path('gse307568_counts'),sep='\t',index_col=0); G.index=G.index.str.split('.').str[0]
Gt=ens_to_symbol(counts_to_log2tpm(G[[s.split('|')[1] for s in cols('GSE307568')]])); Gt.columns=['GSE307568|'+c for c in Gt.columns]; mats['GSE307568']=Gt
B=pd.read_csv(cfg_path('emtab10664_counts'),sep='\t',index_col=0)
lin=(2**B-1).clip(lower=0)
sl=pc.drop_duplicates('ens').assign(L=lambda d: L.reindex(d.ens).values).dropna(subset=['L']).groupby('gene_name').L.median()
lin=lin[lin.index.isin(sl.index)]; tpm=lin.div(sl.loc[lin.index],axis=0); tpm=tpm/tpm.sum()*1e6
Bt=np.log2(tpm+1).groupby(level=0).max(); Bt.columns=['Berg2021|'+c for c in Bt.columns]
Bt=Bt[[s for s in cols('Berg2021_EMTAB10664')]]; mats['Berg']=Bt
tcga_genes=set(sym_of.values)          # TCGA reference: all GENCODE v36 protein-coding gene names
genes=set.intersection(*[set(v.index) for v in mats.values()])&tcga_genes
X=pd.concat([v.loc[sorted(genes)] for v in mats.values()],axis=1).astype('float32')
NG=pd.read_csv(cfg_path('nextgen_dir')/'next_gen_expression.csv',index_col=0); NG.columns=NG.columns.str.replace(r' \(\d+\)$','',regex=True)
ng=[s.split('|')[1] for s in cols('DepMap_NextGen')]
N=NG.loc[ng].T.groupby(level=0).max(); N.columns=['NextGen|'+c for c in N.columns]
C8=pd.read_csv(cfg_path('gse208216_counts'),sep='\t',index_col=0); C8=C8[[c for c in C8.columns if c.startswith('PDO')]]; C8.index=C8.index.str.split('.').str[0]
G8=ens_to_symbol(counts_to_log2tpm(C8)); G8.columns=['GSE208216|'+c for c in G8.columns]
genes=sorted(set(X.index)&set(G8.index)&set(ens2sym.values))
X=pd.concat([X.reindex(genes),N.reindex(genes),G8.reindex(genes)],axis=1)
X=X[list(A['sample'])]
out=cfg_path('harmonised_matrix'); out.parent.mkdir(parents=True,exist_ok=True); X.to_parquet(out)
print('symbol-keyed model matrix: genes',X.shape[0],'samples',X.shape[1])
