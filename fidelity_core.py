"""Core matrix operations for Method A and Method B with explicit gene alignment (round 4, 2026-10-08).

Every matrix operation first reindexes all operands to one explicit gene list and asserts that their gene order is
identical (`aligned`). Method B (Celligner-like) was previously exposed to a silent re-ordering: pandas aligns on labels
and returns a sorted union when the operands' gene sets differ (`legacy_bug=True` reproduces that earlier code path
for the regression test only)."""
import numpy as np, pandas as pd

def aligned(genes, *objs):
    """Reindex every Series/DataFrame to `genes` (rows) and assert identical gene order. Returns the reindexed objects."""
    genes=list(genes)
    out=[o.reindex(genes) for o in objs]
    for o in out: assert list(o.index)==genes, 'gene order mismatch'
    return out if len(out)>1 else out[0]

def assert_same_order(*objs):
    ref=list(objs[0].index)
    for o in objs[1:]: assert list(o.index)==ref, 'gene order mismatch'

def centroids(T,proj,genes):
    Tg=aligned(genes,T)
    C=Tg.T.groupby(np.asarray(proj)).mean().T
    assert list(C.index)==list(genes); return C

def method_a(cent,Xm,genes=None):
    """Spearman correlation of each model profile (columns of Xm) with each centroid, over `genes`."""
    g=list(genes) if genes is not None else list(cent.index)
    c,x=aligned(g,cent,Xm)
    rc=c.rank(); R=x.rank()
    assert_same_order(rc,R)
    return pd.DataFrame({s:rc.corrwith(R[s]) for s in R.columns}).T

def method_b_fit(T,proj,cent,genes,E_dm,dm_lineage,S_org,org_type,seed=20261006,ncp=4,legacy_bug=False):
    """Fit the Celligner-like Method B: centring vectors and `ncp` contrastive PCs (top eigenvectors of
    Cov(TCGA, 40/project) - Cov(all DepMap models)). Returns a dict."""
    g=list(genes)
    proj=pd.Series(np.asarray(proj),index=T.columns)
    Eb=E_dm.reindex(g)                                   # genes absent from the DepMap file -> NaN rows
    mu_cl=Eb.T.groupby(np.asarray(dm_lineage)).mean().mean()
    Sg=S_org.reindex(g); mu_org=Sg.T.groupby(np.asarray(org_type)).mean().mean()
    c=aligned(g,cent); mu_tum=c.mean(axis=1)
    Tg=aligned(g,T)
    mu_cl,mu_org,mu_tum=aligned(g,mu_cl,mu_org,mu_tum)
    Tc=Tg.sub(mu_tum,axis=0); assert_same_order(Tc,mu_tum)
    rng=np.random.default_rng(seed)
    pr=sorted(proj.unique()); sizes={p:min(40,int((proj==p).sum())) for p in pr}
    bal=np.concatenate([rng.choice(np.where(proj.values==p)[0],sizes[p],replace=False) for p in pr])
    if legacy_bug:   # earlier code path (pre-2026-10-08): re-sorts genes when any feature gene is all-NaN
        Ec=Eb.dropna(how='all').sub(mu_cl,axis=0).fillna(0)
        Ct=np.cov(Tc.values[:,bal]); Cm=np.cov(Ec.values)
        w,V=np.linalg.eigh(Ct-Cm); cp=V[:,np.argsort(w)[::-1][:ncp]]
        gc=g
    else:
        # contrastive PCA on the explicit list of feature genes present in the DepMap file (absent genes have no model
        # covariance and would otherwise appear tumour-specific); loadings of absent genes are 0.
        gc=[x for x in g if not Eb.loc[x].isna().all()]
        Ec=aligned(gc,Eb.sub(mu_cl,axis=0)).fillna(0); Tcc=aligned(gc,Tc); assert_same_order(Ec,Tcc)
        Ct=np.cov(Tcc.values[:,bal]); Cm=np.cov(Ec.values)
        w,V=np.linalg.eigh(Ct-Cm); cpc=V[:,np.argsort(w)[::-1][:ncp]]
        cp=pd.DataFrame(cpc,index=gc).reindex(g).fillna(0.0).values
    TcB=pd.DataFrame(Tc.values-cp@(cp.T@Tc.values),index=g,columns=Tc.columns)
    centB=TcB.T.groupby(proj.values).mean().T
    assert list(centB.index)==g
    return dict(genes=g,cpca_genes=gc,mu_cl=mu_cl,mu_org=mu_org,mu_tum=mu_tum,cp=cp,centB=centB,
                absent_from_depmap=int(Eb.isna().all(axis=1).sum()),gene_order_ok=list(Ec.index)==list(gc))

def method_b_score(fit,Xm,centring):
    """Xm: genes x samples; centring: Series sample -> 'tumour' | 'depmap' | 'organoid'."""
    g=fit['genes']; x=aligned(g,Xm); cp=fit['cp']; centB=fit['centB']
    out={}
    for s in x.columns:
        k=centring[s]
        mu={'tumour':fit['mu_tum'],'depmap':fit['mu_cl'],'organoid':fit['mu_org']}[k]
        assert_same_order(x,mu,centB)
        v=(x[s]-mu).values
        if k!='tumour': v=np.nan_to_num(v,nan=0.0)
        v=v-cp@(cp.T@v)
        out[s]=pd.Series([np.corrcoef(centB[c].values,v)[0,1] for c in centB.columns],index=centB.columns)
    return pd.DataFrame(out).T
