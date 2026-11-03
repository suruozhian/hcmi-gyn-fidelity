#!/usr/bin/env python3
"""Diagnostic for the real-data regression test (round 4): decomposes the Method B call changes of selected perturbations
into (a) the effect on the contrastive components / centring (genes removed from the DepMap file used for the fit, DepMap
unit profiles intact) and (b) the effect of the missing genes in the DepMap unit profiles themselves (fit on the full file).
Same random draws as test_gene_alignment.py (seed 20261008). Informational; not a pass/fail test.
Run: python3 tests/diagnose_regression.py --tag 26Q3 --depmap-expr <file> --cases 5:0,20:0,20:4 --out <json>
"""
import sys, argparse, json, re
from pathlib import Path
import numpy as np, pandas as pd
PIPE=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(PIPE)); sys.path.insert(0,str(Path(__file__).resolve().parent))
from fidelity_core import aligned, method_a, method_b_fit, method_b_score, centroids
from test_gene_alignment import subspace_cc
def diag(tag,depmap_expr,cases):
    from paths import path as cfg_path
    import genemap
    REF=cfg_path('reference_dir')
    T=pd.read_parquet(REF/'tcga11_full_entrez_log2tpm.parquet'); TS=pd.read_csv(REF/'tcga11_full_samples.csv').set_index('sample')
    X=pd.read_parquet(REF/'expr_entrez_log2tpm_models.parquet'); A=pd.read_csv(cfg_path('sample_annotations')).set_index('sample')
    E=pd.read_csv(depmap_expr,index_col=0,low_memory=False)
    if 'ModelID' in E.columns: E=E[E.IsDefaultEntryForModel.astype(str)=='Yes'].set_index('ModelID')
    gc=[c for c in E.columns if re.search(r'\(\d+\)$',str(c))]; E=E[gc].astype('float32').T
    E.index=genemap.depmap_entrez(E.index); E=E.groupby(level=0).max()
    top=pd.read_csv(PIPE/'results'/tag/'feature_genes_top2000.csv').iloc[:,0].astype(str).tolist()
    genes=X.index.intersection(T.index); T=T.loc[genes]; X=X.loc[genes]
    proj=TS.loc[T.columns,'project'].values
    U=pd.read_csv(PIPE/'inputs/units_215.csv')
    U['samples']=[list(A.index[(A.dataset==r.dataset)&(A.model_unit==r.model_unit)]) for r in U.itertuples()]
    cols=sorted({s for l in U.samples for s in l})
    dm=[s for s in cols if A.loc[s,'dataset'] in ('DepMap_24Q4','DepMap_NextGen')]
    mdl=pd.read_csv(cfg_path('depmap_model_csv')).set_index('ModelID')
    S=pd.read_csv(cfg_path('sanger_organoid_rna'),index_col=0); S.columns=S.columns.str.replace('_exp$','',regex=True)
    S=genemap.to_entrez(S.T,genemap.symbol_map(S.columns).entrez)
    da=pd.read_excel(cfg_path('sanger_organoid_annotation'),'data_availability').set_index('sample_ID')
    st=da.reindex(S.columns).primary_tumour_type.fillna('NA').values
    cent=centroids(T,proj,top)
    centring=pd.Series({s:('tumour' if A.loc[s,'kind'] in ('reference_tumour','tumour') else 'depmap' if A.loc[s,'dataset']=='DepMap_24Q4' else 'organoid') for s in cols})
    def run(Efit,Eprof,legacy=False):
        Xs=X[cols].copy()
        for s in dm:
            m=A.loc[s,'model_unit']
            if m in Eprof.columns: Xs[s]=Eprof[m].reindex(Xs.index).values
        dl=mdl.reindex(Efit.columns).OncotreeLineage.fillna('NA').values
        f=method_b_fit(T,proj,cent,top,Efit,dl,S,st,legacy_bug=legacy)
        CA=method_a(cent,Xs,top); CB=method_b_score(f,Xs,centring)
        na=pd.Series([CA.loc[l].mean().idxmax() for l in U.samples]); nb=pd.Series([CB.loc[l].mean().idxmax() for l in U.samples])
        return f,na,nb
    f0,na0,nb0=run(E,E)
    ds=pd.Series([A.loc[l[0],'dataset'] for l in U.samples])
    rng=np.random.default_rng(20261008); present=[g for g in top if g in E.index]; out=[]
    for k in (5,20,36):
        for rep in range(5):
            drop=list(rng.choice(present,k,replace=False))
            if (k,rep) not in cases: continue
            E1=E.drop(index=drop)
            for lab,(ef,ep) in {'both (as in test)':(E1,E1),'cPCA/centring only (DepMap profiles intact)':(E1,E),'DepMap profiles only (cPCA on full file)':(E,E1)}.items():
                f1,na1,nb1=run(ef,ep); chB=nb1!=nb0; chA=na1!=na0
                out.append(dict(k=k,rep=rep,perturbation=lab,changed_A=int(chA.sum()),changed_B=int(chB.sum()),
                                changed_B_by_dataset=ds[chB].value_counts().to_dict(),cc=round(subspace_cc(f0,f1),4)))
                print(out[-1],flush=True)
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--tag',default='26Q3'); ap.add_argument('--depmap-expr',required=True)
    ap.add_argument('--cases',default='5:0,20:0'); ap.add_argument('--out',default=None); a=ap.parse_args()
    cases={tuple(int(x) for x in c.split(':')) for c in a.cases.split(',')}
    out=diag(a.tag,a.depmap_expr,cases)
    if a.out: json.dump(out,open(a.out,'w'),indent=1)
