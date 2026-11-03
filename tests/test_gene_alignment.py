#!/usr/bin/env python3
"""Regression test for gene alignment (round 4, 2026-10-08).

Run:   python3 tests/test_gene_alignment.py --synthetic            (no input data needed; seconds)
       python3 tests/test_gene_alignment.py --real --tag 26Q3 --depmap-expr <OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv>
                                                                    (needs config_paths.json inputs and results/<tag>/feature_genes_top2000.csv)

Perturbation: k feature genes are removed from the DepMap expression file (the real failure mode: genes absent from a
DepMap release), k in {5, 20, 36}, 5 random draws each. Method A and Method B are recomputed with the current code and
with the earlier Method B code path (fidelity_core.method_b_fit(legacy_bug=True)).
Pass criteria (fixed before the test was run):
  C1  contrastive-PC subspace preserved: smallest canonical correlation between baseline and perturbed cPC subspaces
      (over the retained genes) >= 0.90;
  C2  (real data only) <= 10 of 215 units (about 5%) change nearest centroid, for Method A and for Method B.
The current code must pass C1 and C2 in every perturbation; the earlier code path must fail (the test reports it).
"""
import sys, argparse, json, re
from pathlib import Path
import numpy as np, pandas as pd
PIPE=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(PIPE))
from fidelity_core import aligned, method_a, method_b_fit, method_b_score, centroids

def subspace_cc(fit0,fit1):
    """smallest canonical correlation between the cPC subspaces of two fits, over genes used by the perturbed fit's cPCA."""
    g=[x for x in fit1['cpca_genes'] if x in fit0['genes']]; i0=[fit0['genes'].index(x) for x in g]; i1=[fit1['genes'].index(x) for x in g]
    U0=fit0['cp'][i0,:]; U1=fit1['cp'][i1,:]
    q0,_=np.linalg.qr(U0); q1,_=np.linalg.qr(U1)
    return float(np.linalg.svd(q0.T@q1,compute_uv=False).min())

def check_aligned_raises():
    a=pd.Series([1.,2.,3.],index=['b','a','c']); g=['a','b','c']
    r=aligned(g,a); assert list(r.index)==g and r['b']==1.0
    bad=pd.Series([1.,2.],index=['a','b'])
    try:
        from fidelity_core import assert_same_order
        assert_same_order(r,bad); return False
    except AssertionError: return True

def synthetic():
    rng=np.random.default_rng(1)
    G=400; genes=[str(x) for x in rng.choice(np.arange(1,200000),G,replace=False)]      # unsorted Entrez-like IDs
    P=11; base=rng.normal(0,1.5,(G,P)); proj=np.repeat([f'P{i}' for i in range(P)],40)
    micro=rng.normal(0,1,(G,4))*np.array([.26,.23,.2,.17])                                       # 4 weaker tumour-only axes (distinct variances); lineage axes shared with models dominate, as in real data
    T=pd.DataFrame(base[:,np.repeat(np.arange(P),40)]+micro@rng.normal(0,1,(4,P*40))+rng.normal(0,.5,(G,P*40)),index=genes)
    T.columns=[f't{i}' for i in range(T.shape[1])]
    load=rng.normal(0,1,(G,3))                                                                  # culture-only axes
    lin=np.repeat(np.arange(P),27); nm=len(lin)                                                # lineage-balanced models
    E=pd.DataFrame(base[:,lin]+load@rng.normal(0,.3,(3,nm))+rng.normal(0,.5,(G,nm)),index=genes)
    E.columns=[f'm{i}' for i in range(nm)]
    S=E.iloc[:,:30].copy(); st=lin[:30]
    cent=centroids(T,proj,genes)
    f0=method_b_fit(T,proj,cent,genes,E,lin,S,st)
    res=[]
    for k in (5,20,36):
        for rep in range(5):
            drop=set(rng.choice(genes,k,replace=False)); E1=E.drop(index=list(drop))
            fx=method_b_fit(T,proj,cent,genes,E1,lin,S,st)
            fl=method_b_fit(T,proj,cent,genes,E1,lin,S,st,legacy_bug=True)
            res.append(dict(k=k,rep=rep,cc_current=round(subspace_cc(f0,fx),3),cc_legacy=round(subspace_cc(f0,fl),3)))
    R=pd.DataFrame(res); print(R.to_string(index=False))
    ok_cur=bool((R.cc_current>=0.90).all()); fail_leg=bool((R.cc_legacy<0.90).any())
    return dict(aligned_raises_on_mismatch=check_aligned_raises(),current_passes_C1=ok_cur,legacy_fails_C1=fail_leg,
                legacy_fail_fraction=float((R.cc_legacy<0.90).mean()),table=R.to_dict('records'))

def real(tag,depmap_expr):
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
    def run(Edm,legacy):
        Xs=X[cols].copy()
        for s in dm:
            m=A.loc[s,'model_unit']
            if m in Edm.columns: Xs[s]=Edm[m].reindex(Xs.index).values
        dl=mdl.reindex(Edm.columns).OncotreeLineage.fillna('NA').values
        f=method_b_fit(T,proj,cent,top,Edm,dl,S,st,legacy_bug=legacy)
        CA=method_a(cent,Xs,top); CB=method_b_score(f,Xs,centring)
        na=pd.Series([CA.loc[l].mean().idxmax() for l in U.samples]); nb=pd.Series([CB.loc[l].mean().idxmax() for l in U.samples])
        return f,na,nb
    f0,na0,nb0=run(E,False)
    rng=np.random.default_rng(20261008); res=[]
    present=[g for g in top if g in E.index]
    for k in (5,20,36):
        for rep in range(5):
            drop=list(rng.choice(present,k,replace=False)); E1=E.drop(index=drop)
            for legacy in (False,True):
                f1,na1,nb1=run(E1,legacy)
                res.append(dict(k=k,rep=rep,path='earlier' if legacy else 'current',cc=round(subspace_cc(f0,f1),3),
                                changed_nearest_A=int((na1!=na0).sum()),changed_nearest_B=int((nb1!=nb0).sum())))
            print(res[-2],res[-1],flush=True)
    R=pd.DataFrame(res)
    cur=R[R.path=='current']; leg=R[R.path!='current']
    ok=lambda d: (d.cc>=0.90)&(d.changed_nearest_A<=10)&(d.changed_nearest_B<=10)
    return dict(tag=tag,depmap_expr=Path(depmap_expr).name,feature_genes_absent_in_file=len(top)-len(present),
                current_passes_all=bool(ok(cur).all()),legacy_fails_any=bool((~ok(leg)).any()),legacy_fail_fraction=float((~ok(leg)).mean()),
                current_max_changed_A=int(cur.changed_nearest_A.max()),current_max_changed_B=int(cur.changed_nearest_B.max()),current_min_cc=float(cur.cc.min()),
                legacy_max_changed_B=int(leg.changed_nearest_B.max()),legacy_min_cc=float(leg.cc.min()),table=R.to_dict('records'))

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--synthetic',action='store_true'); ap.add_argument('--real',action='store_true')
    ap.add_argument('--tag',default='26Q3'); ap.add_argument('--depmap-expr',default=None); ap.add_argument('--out',default=None)
    a=ap.parse_args(); out={}
    if a.synthetic or not a.real: out['synthetic']=synthetic()
    if a.real: out['real']=real(a.tag,a.depmap_expr)
    summ={k:{kk:vv for kk,vv in v.items() if kk!='table'} for k,v in out.items()}
    print(json.dumps(summ,indent=1))
    if a.out: json.dump(out,open(a.out,'w'),indent=1)
    passed=all(v.get('current_passes_C1',v.get('current_passes_all')) and v.get('legacy_fails_C1',v.get('legacy_fails_any')) for v in out.values())
    if 'synthetic' in out: passed=passed and out['synthetic']['aligned_raises_on_mismatch']
    print('RESULT:', 'PASS' if passed else 'FAIL'); sys.exit(0 if passed else 1)
