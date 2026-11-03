#!/usr/bin/env python3
"""One-command recompute of lineage fidelity (Method A/B), symmetric strict/loose Table 2, confusion matrices,
equal-depth and reference-composition sensitivities, CPTAC transfer, TCGA LOO, histotype calls, marker (WT1/KRT20/CDX2)
values and Table S2, for a given DepMap expression file.

Usage:
  python3 run_all.py                                   # DepMap 24Q4 (default expression; sensitivity version)
  python3 run_all.py --depmap-expr /path/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv --tag 26Q3
Then: python3 make_outputs.py --tag <tag> [--publish]   (figures, xlsx, number-audit); run_all.py calls it automatically.

Gene key (round 4): --gene-id entrez (default, primary) | symbol (sensitivity). All matrix operations go through
fidelity_core.aligned (explicit gene list + identical-order assertion).

Fixed inputs (no sample exclusions): the 215 independent expression units in inputs/units_215.csv
(dataset, model_unit, stratum). Rules: strict = ovary-labelled -> OV; EC (non-CS) -> UCEC; carcinosarcoma -> UCS.
loose = nearest in {OV, UCEC, UCS}. Headline reference = 11 TCGA projects, full open cohorts (1 primary tumour per case).
"""
import argparse, json, re, sys, subprocess, hashlib
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import f_oneway
import gzip

PIPE=Path(__file__).resolve().parent
sys.path.insert(0,str(PIPE))
from paths import path as cfg_path   # input locations come from config_paths.json (see paths.py)
import genemap
from fidelity_core import aligned, assert_same_order, centroids as core_centroids, method_a, method_b_fit, method_b_score
DEFAULT_EXPR=cfg_path('depmap_expr_24q4')   # DepMap Public 24Q4
DEFAULT_MODEL=cfg_path('depmap_model_csv')  # DepMap Public 26Q3 Model.csv
GYN=['OV','UCEC','UCS']
SEED=20261008

ap=argparse.ArgumentParser()
ap.add_argument('--depmap-expr',default=str(DEFAULT_EXPR))
ap.add_argument('--tag',default=None)
ap.add_argument('--reps',type=int,default=25)
ap.add_argument('--per-project',type=int,default=40)
ap.add_argument('--no-outputs',action='store_true')
ap.add_argument('--reference-dir',default=None)
ap.add_argument('--expect-models',type=int,default=None)
ap.add_argument('--publish',action='store_true',help='also copy figures/tables/audit into the manuscript package')
ap.add_argument('--gene-id',choices=['entrez','symbol'],default='entrez',help='gene key: Entrez Gene ID (primary) or gene symbol (sensitivity)')
a=ap.parse_args()
GID=a.gene_id
expr_path=Path(a.depmap_expr)
tag=a.tag or (('24Q4' if expr_path==DEFAULT_EXPR else expr_path.stem)+('' if GID=='entrez' else '_symbol'))
OUT=PIPE/'results'/tag; OUT.mkdir(parents=True,exist_ok=True)
LOG=[]
def p(*x):
    s=' '.join(map(str,x)); print(s,flush=True); LOG.append(s)

def md5(path):
    h=hashlib.md5()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1<<22),b''): h.update(b)
    return h.hexdigest()

# ───────────── 1. model matrix (non-TCGA columns of the harmonised matrix) ─────────────
REF=Path(a.reference_dir) if a.reference_dir else cfg_path('reference_dir')
HARM=REF/'expr_entrez_log2tpm_models.parquet' if GID=='entrez' else cfg_path('harmonised_matrix')
X=pd.read_parquet(HARM)
A=pd.read_csv(cfg_path('sample_annotations')).set_index('sample')
X=X[[c for c in X.columns if A.loc[c,'dataset']!='TCGA']].astype('float32')
A=A.loc[X.columns]

DEPMAP_QC={}
def load_depmap(path):
    head=pd.read_csv(path,nrows=2,index_col=0)
    if 'ModelID' in head.columns:            # profile-level format (26Q3): keep default entry per model only
        E=pd.read_csv(path,index_col=0,low_memory=False)
        DEPMAP_QC['rows_total']=len(E)
        assert 'IsDefaultEntryForModel' in E.columns
        E=E[E['IsDefaultEntryForModel'].astype(str)=='Yes']
        DEPMAP_QC['rows_default_model']=len(E); DEPMAP_QC['unique_models']=int(E.ModelID.nunique())
        assert E.ModelID.is_unique, 'more than one default row per model'
        if a.expect_models: assert len(E)==a.expect_models, f'expected {a.expect_models} default rows, got {len(E)}'
        E=E.set_index('ModelID')
        gene_cols=[c for c in E.columns if re.search(r'\(\d+\)$',str(c))]
        DEPMAP_QC['non_gene_columns']=[c for c in E.columns if c not in gene_cols]
        E=E[gene_cols]
    else:
        E=pd.read_csv(path,index_col=0)
    if GID=='entrez':
        sym=pd.Index(genemap.depmap_entrez(E.columns))
        DEPMAP_QC['gene_matching']='Entrez Gene ID from column header'
        keep=~sym.isna(); E=E.loc[:,keep]; sym=sym[keep]
    else:
        sym=pd.Index(E.columns).str.replace(r' \(\d+\)$','',regex=True)
        DEPMAP_QC['gene_matching']='stripped gene symbol'
    dup=sym[sym.duplicated(keep=False)]
    DEPMAP_QC['gene_columns']=len(sym); DEPMAP_QC['duplicate_gene_keys']=sorted(set(dup))
    DEPMAP_QC['duplicate_handling']='max over columns sharing a gene key (same rule as the harmonised matrix)'
    E.columns=sym
    E=E.astype('float32').T.groupby(level=0).max()       # genes x models
    return E
p(f'DepMap expression: {expr_path.name}  (md5 {md5(expr_path)})')
E_all=load_depmap(expr_path)
p('DepMap expression matrix genes x models',E_all.shape)
p('DepMap QC:',{k:v for k,v in DEPMAP_QC.items()})
REQ=['ACH-003350','ACH-003723','ACH-003709','ACH-003710','ACH-003711','ACH-002419','ACH-003903']
DEPMAP_QC['required_models_present']={m:bool(m in E_all.columns) for m in REQ}
p('required models with default rows:',DEPMAP_QC['required_models_present'])
units_dm=pd.read_csv(PIPE/'inputs/units_215.csv'); units_dm=units_dm[units_dm.dataset.isin(['DepMap_24Q4','DepMap_NextGen'])]
DEPMAP_QC['unit_models_missing_in_file']=sorted(set(units_dm.model_unit)-set(E_all.columns))
DEPMAP_QC['unit_models_present_in_file']=sorted(set(units_dm.model_unit)&set(E_all.columns))
p('DepMap/NextGen unit models absent from this file (previous values retained):',DEPMAP_QC['unit_models_missing_in_file'])
# swap DepMap 2D and NextGen sample columns where the model is present in the given file
src={}
for s in X.columns:
    ds=A.loc[s,'dataset']
    if ds in ('DepMap_24Q4','DepMap_NextGen'):
        mid=A.loc[s,'model_unit']
        if mid in E_all.columns:
            old=X[s].copy()
            X[s]=E_all[mid].reindex(X.index).values
            src[s]=f'{tag}' ; 
        else:
            src[s]='24Q4 matrix (DepMap 2D)' if ds=='DepMap_24Q4' else 'NextGen manuscript file set'
    else:
        src[s]='resource-native'
A['expr_source']=pd.Series(src)
sw=A[A.dataset.isin(['DepMap_24Q4','DepMap_NextGen'])].groupby(['dataset','expr_source']).size()
p('DepMap/NextGen expression source per sample:\n'+sw.to_string())
# sanity for default: swapped 24Q4 values identical to harmonised matrix
X0=pd.read_parquet(HARM,columns=[s for s in X.columns if A.loc[s,'dataset']=='DepMap_24Q4'])
diff=np.nanmax(np.abs(X0.values-X[X0.columns].values)); del X0
p(f'max |delta| DepMap 2D vs harmonised 24Q4 values ({GID}): {diff:.3g}')
nan_frac=X[[s for s in X.columns if A.loc[s,'dataset']=='DepMap_24Q4']].isna().mean().mean()
p(f'DepMap 2D missing-gene fraction after swap: {nan_frac:.4f}')

# ───────────── 2. TCGA references ─────────────
T=pd.read_parquet(REF/f'tcga11_full_{GID}_log2tpm.parquet')
TS=pd.read_csv(REF/'tcga11_full_samples.csv').set_index('sample')
genes=X.index.intersection(T.index)
T=T.loc[genes]; X=X.loc[genes]
assert_same_order(T,X)
proj=TS.loc[T.columns,'project']
p('TCGA reference (11 projects, full):',proj.value_counts().sort_index().to_dict(),'genes',len(genes))
pur_sym=pd.read_csv(cfg_path('purity_filtered_genes')).iloc[:,0]
pur=set(genemap.gencode_symbols_to_entrez(pur_sym)) if GID=='entrez' else set(pur_sym)
p('purity-filtered genes',len(pur_sym),'->',len(pur),GID)
E2S=genemap.entrez2symbol() if GID=='entrez' else None
def disp(g): return [E2S.get(x,x) for x in g] if GID=='entrez' else list(g)
def feature_genes(Tref):
    expd=Tref.index[Tref.mean(axis=1)>1]; keep=[g for g in expd if g not in pur]
    return Tref.loc[keep].var(axis=1).sort_values(ascending=False).index[:2000]
top=feature_genes(T)
pd.DataFrame({'gene':top,'symbol':disp(top)}).to_csv(OUT/'feature_genes_top2000.csv',index=False)
legacy_top=pd.read_csv(cfg_path('legacy_feature_genes')).iloc[:,0]
if GID=='symbol': p('feature genes',len(top),'; overlap with legacy top2000:',len(set(top)&set(legacy_top)))
dm_absent_top=int(sum(1 for g in top if g not in E_all.index))
p('feature genes absent from this DepMap file:',dm_absent_top)

U=pd.read_csv(PIPE/'inputs/units_215.csv')
assert len(U)==215
U['samples']=[list(A.index[(A.dataset==r.dataset)&(A.model_unit==r.model_unit)]) for r in U.itertuples()]
assert U.samples.map(len).min()>0
unit_samples=sorted({s for l in U.samples for s in l})
cptac=list(A.index[A.dataset=='CPTAC'])

def methodA(cent,cols,genes_=None):
    return method_a(cent,X[cols],genes_ if genes_ is not None else cent.index)

def to_units(C):
    rows=[]
    for r in U.itertuples():
        m=C.loc[r.samples].mean(); rows.append(m)
    M=pd.DataFrame(rows,index=U.index); return M

def calls(M):
    near=M.idxmax(axis=1)
    strict=np.where(U.stratum.isin(['HGSC','serous_NOS','nonHGSC']),near=='OV',
            np.where(U.stratum=='EC',near=='UCEC',np.where(U.stratum=='CS',near=='UCS',False)))
    loose=near.isin(GYN).values
    return near,pd.Series(strict,index=U.index),pd.Series(loose,index=U.index)

def centroids(Tref,pr,g): return core_centroids(Tref,np.asarray(pr),g)

# ───────────── 3. Headline: Method A on 11-project full reference ─────────────
cent=centroids(T,proj,top)
CA=methodA(cent,unit_samples+cptac)
MA=to_units(CA); nearA,strictA,looseA=calls(MA)
U['nearest_A']=nearA; U['strict_A']=strictA; U['loose_A']=looseA
for c in cent.columns: U['rA_'+c]=MA[c].round(4)
p('Method A done')

# ───────────── 4. Method B (Celligner-like) ─────────────
mdl=pd.read_csv(DEFAULT_MODEL).set_index('ModelID')
dl=mdl.reindex(E_all.columns).OncotreeLineage.fillna('NA')
S=pd.read_csv(cfg_path('sanger_organoid_rna'),index_col=0); S.columns=S.columns.str.replace('_exp$','',regex=True)
S=S.T
if GID=='entrez': S=genemap.to_entrez(S,genemap.symbol_map(S.index).entrez)
da=pd.read_excel(cfg_path('sanger_organoid_annotation'),'data_availability').set_index('sample_ID')
st=da.reindex(S.columns).primary_tumour_type.fillna('NA')
FITB=method_b_fit(T,proj.values,cent,top,E_all,dl.values,S,st.values)
METHODB_QC=dict(gene_order_fixed=True,gene_order_ok=FITB['gene_order_ok'],feature_genes_absent_from_depmap=FITB['absent_from_depmap'])
def centring(s):
    k=A.loc[s,'kind']; d=A.loc[s,'dataset']
    return 'tumour' if k in ('reference_tumour','tumour') else ('depmap' if d=='DepMap_24Q4' else 'organoid')
cols_=unit_samples+cptac
CB=method_b_score(FITB,X[cols_],pd.Series({s:centring(s) for s in cols_}))
MB=to_units(CB); nearB,strictB,looseB=calls(MB)
U['nearest_B']=nearB; U['strict_B']=strictB; U['loose_B']=looseB
p('Method B done')

# ───────────── 5. Table 2 ─────────────
U['is_org']=U.is_org.astype(bool)
def row(label,m):
    d=U[m]; o=d[d.is_org]
    return dict(stratum=label,n=len(d),strict_A=int(d.strict_A.sum()),strict_A_pct=round(100*d.strict_A.mean(),1),
                loose_A=int(d.loose_A.sum()),loose_A_pct=round(100*d.loose_A.mean(),1),strict_B=int(d.strict_B.sum()),loose_B=int(d.loose_B.sum()),
                organoid_n=len(o),organoid_strict_A=int(o.strict_A.sum()),organoid_loose_A=int(o.loose_A.sum()),
                twoD_n=len(d)-len(o),twoD_strict_A=int(d.strict_A.sum()-o.strict_A.sum()),twoD_loose_A=int(d.loose_A.sum()-o.loose_A.sum()))
hg=U.stratum.isin(['HGSC','serous_NOS'])
T2=pd.DataFrame([
 row('HGSC/serous (incl. serous NOS)',hg),
 row('  OV_HGSC',U.histology=='OV_HGSC'),
 row('  HCMI/Sanger OV_serous',U.histology=='OV_serous'),
 row('  DepMap serous NOS',U.stratum=='serous_NOS'),
 row('  HGSC/serous organoids',hg&U.is_org),
 row('  HGSC/serous 2D lines',hg&~U.is_org),
 row('Endometrial carcinoma (non-CS)',U.stratum=='EC'),
 row('  EC organoids',(U.stratum=='EC')&U.is_org),
 row('  EC 2D lines',(U.stratum=='EC')&~U.is_org),
 row('Uterine carcinosarcoma',U.stratum=='CS'),
 row('Ovarian non-HGSC (coverage gap)',U.stratum=='nonHGSC'),
 row('All ovary',U.lineage=='ovary'),
 row('All uterus',U.lineage=='uterus'),
 row('All 215 units',U.index==U.index)])
T2.to_csv(OUT/'table2_strict_loose.csv',index=False)
p('\nTable 2 (Method A strict|loose; B strict|loose):\n'+T2[['stratum','n','strict_A','loose_A','strict_B','loose_B','organoid_n','organoid_strict_A','twoD_n','twoD_strict_A']].to_string(index=False))

# confusion matrices
for tag_,col in [('A','nearest_A'),('B','nearest_B')]:
    cm=pd.crosstab(U.histology,U[col],margins=True,margins_name='All')
    cm=cm.reindex(columns=[c for c in sorted(cent.columns) if c in cm.columns]+['All'],fill_value=0)
    cm.to_csv(OUT/f'confusion_matrix_method{tag_}.csv')
    cs=pd.crosstab(U.stratum,U[col]); cs.to_csv(OUT/f'confusion_by_stratum_method{tag_}.csv')
# extra descriptive numbers
ov=U[U.lineage=='ovary']
extra=dict(
 ovary_nearest_UCS=int((ov.nearest_A=='UCS').sum()),ovary_nearest_UCEC=int((ov.nearest_A=='UCEC').sum()),ovary_n=len(ov),
 CS_nearest=U[U.stratum=='CS'].nearest_A.value_counts().to_dict(),
 EC2D_loose_ids=';'.join(U[(U.stratum=='EC')&~U.is_org&U.loose_A].model_unit),
 EC2D_strict_ids=';'.join(U[(U.stratum=='EC')&~U.is_org&U.strict_A].model_unit),
 serousNOS_nearest=U[U.stratum=='serous_NOS'].nearest_A.value_counts().to_dict())
U[(U.stratum=='EC')&~U.is_org][['dataset','model_unit','histology','nearest_A','strict_A','loose_A','nearest_B']].to_csv(OUT/'EC_2D_units_calls.csv',index=False)

# ───────────── 6. CPTAC transfer + TCGA LOO ─────────────
Cc=CA.loc[cptac]; nc=Cc.idxmax(axis=1); lin=A.loc[cptac,'lineage']
cpt=dict(OV_n=int((lin=='ovary').sum()),OV_nearestOV=int(((lin=='ovary')&(nc=='OV')).sum()),
         UCEC_n=int((lin=='uterus').sum()),UCEC_nearest_UCEC_or_UCS=int(((lin=='uterus')&nc.isin(['UCEC','UCS'])).sum()),
         UCEC_nearest_UCEC=int(((lin=='uterus')&(nc=='UCEC')).sum()))
p('CPTAC transfer',cpt)
M=aligned(top,T); assert_same_order(M,cent); RM=M.rank().values; rc=cent.rank()
ok={}
for pp in cent.columns:
    idx=np.where(proj.values==pp)[0]; tot=M.iloc[:,idx].sum(axis=1).values; n=len(idx)
    base=rc.values
    for i in idx:
        cc=cent.values.copy(); j=list(cent.columns).index(pp); cc[:,j]=(tot-M.iloc[:,i].values)/(n-1)
        rcc=pd.DataFrame(cc).rank().values; x=RM[:,i]
        r=[np.corrcoef(rcc[:,k],x)[0,1] for k in range(cc.shape[1])]
        ok[T.columns[i]]=cent.columns[int(np.argmax(r))]==pp
loo=pd.Series(ok); looA=loo.groupby(proj.reindex(loo.index)).mean().round(3)
p('TCGA LOO A (11 full):',looA.to_dict())
looA.to_csv(OUT/'tcga_loo_methodA.csv')

# ───────────── 7. Equal-depth sensitivity (40/project x 25, all 11 projects, full cohorts) ─────────────
rng=np.random.default_rng(SEED); reps=[]; strata_def={'HGSC_serous':hg,'EC_nonCS':U.stratum=='EC','EC_organoid':(U.stratum=='EC')&U.is_org,
  'EC_2D':(U.stratum=='EC')&~U.is_org,'HGSC_organoid':hg&U.is_org,'HGSC_2D':hg&~U.is_org,'CS':U.stratum=='CS','nonHGSC':U.stratum=='nonHGSC'}
for k in range(a.reps):
    idx=np.concatenate([rng.choice(np.where(proj.values==pp)[0],a.per_project,replace=False) for pp in sorted(proj.unique())])
    c_=centroids(T.iloc[:,idx],proj.iloc[idx],top)
    n_,s_,l_=calls(to_units(methodA(c_,unit_samples)))
    r={'rep':k}
    for nm,m in strata_def.items():
        r[nm+'_n']=int(m.sum()); r[nm+'_strict']=int(s_[m].sum()); r[nm+'_loose']=int(l_[m].sum())
    reps.append(r)
R=pd.DataFrame(reps); R.to_csv(OUT/'equal_depth_reps.csv',index=False)
summ=[]
for nm in strata_def:
    for kind in ['strict','loose']:
        v=R[f'{nm}_{kind}']; n=R[f'{nm}_n'].iloc[0]
        summ.append(dict(stratum=nm,rule=kind,n=n,mean_count=v.mean(),min_count=v.min(),max_count=v.max(),
                         mean_pct=round(100*v.mean()/n,1),min_pct=round(100*v.min()/n,1),max_pct=round(100*v.max()/n,1)))
ES=pd.DataFrame(summ); ES.to_csv(OUT/'equal_depth_summary.csv',index=False)
p('\nEqual-depth 40/project x',a.reps,':\n'+ES.to_string(index=False))

# ───────────── 8. Reference-composition sensitivity (same 215 units, same rules, same feature genes) ─────────────
ng40=pd.read_csv(cfg_path('tcga_nongyn_subsample_manifest')); ng40_cases=set(ng40.case)
ph1=pd.read_csv(cfg_path('tcga_gyn_subsample_manifest')); ph1_gyn=set(ph1[ph1.project.isin(['TCGA-OV','TCGA-UCEC','TCGA-UCS'])].case)
case=TS.loc[T.columns,'case']; isg=proj.isin(GYN)
configs={
 'Full 11 projects (headline)':np.ones(len(proj),bool),
 'Gyn full + non-gyn 40/project':(isg|case.isin(ng40_cases)).values,
 'Gyn ~200/project + non-gyn full':((isg&case.isin(ph1_gyn))|~isg).values,
}
sens=[]
def add_sens(label,s_,l_,ref_desc,n_ref):
    r=dict(reference=label,ref_tumours=n_ref,ref_desc=ref_desc)
    for nm,m in [('HGSC_serous',hg),('EC_nonCS',U.stratum=='EC'),('EC_organoid',(U.stratum=='EC')&U.is_org),('EC_2D',(U.stratum=='EC')&~U.is_org),('CS',U.stratum=='CS'),('nonHGSC',U.stratum=='nonHGSC'),('uterus_all',U.lineage=='uterus')]:
        r[nm+'_n']=int(m.sum()); r[nm+'_strict']=float(s_[m].sum()); r[nm+'_loose']=float(l_[m].sum())
    sens.append(r)
for label,mask in configs.items():
    c_=centroids(T.loc[:,mask],proj[mask],top)
    n_,s_,l_=calls(to_units(methodA(c_,unit_samples)))
    desc=proj[mask].value_counts().sort_index().to_dict()
    add_sens(label,s_,l_,json.dumps(desc),int(mask.sum()))
# equal-depth mean as a row
r=dict(reference=f'Equal 40/project (mean of {a.reps})',ref_tumours=40*11,ref_desc='40 per project, random')
for nm in ['HGSC_serous','EC_nonCS','EC_organoid','EC_2D','CS','nonHGSC']:
    r[nm+'_n']=int(R[nm+'_n'].iloc[0]); r[nm+'_strict']=R[nm+'_strict'].mean(); r[nm+'_loose']=R[nm+'_loose'].mean()
sens.append(r)
SE=pd.DataFrame(sens); SE.to_csv(OUT/'reference_composition_sensitivity.csv',index=False)
p('\nReference composition sensitivity:\n'+SE.drop(columns=['ref_desc']).to_string(index=False))

# legacy reproduction check: previous headline reference (gyn full + non-gyn 40) with legacy feature genes
mask=configs['Gyn full + non-gyn 40/project']
if GID=='symbol':
  lt=[g for g in legacy_top if g in T.index]
  c_=centroids(T.loc[:,mask],proj[mask],lt)
  n_,s_,l_=calls(to_units(methodA(c_,unit_samples,lt)))
  leg=dict(HGSC_serous=f'{int(s_[hg].sum())}/{int(hg.sum())} | {int(l_[hg].sum())}',EC=f"{int(s_[U.stratum=='EC'].sum())}/83 | {int(l_[U.stratum=='EC'].sum())}",
         EC_org=f"{int(s_[(U.stratum=='EC')&U.is_org].sum())}/57",EC_2D=f"{int(s_[(U.stratum=='EC')&~U.is_org].sum())}/26 | {int(l_[(U.stratum=='EC')&~U.is_org].sum())}",
         CS=f"{int(s_[U.stratum=='CS'].sum())}/18 | {int(l_[U.stratum=='CS'].sum())}")
  p('Legacy reproduction (gyn full + non-gyn 40, legacy genes):',leg)
else: leg='not applicable (legacy feature genes are symbol-keyed; see the symbol-matching sensitivity run)'

# ───────────── 9. Histotype calls (GSE295399 centroids; three centring regimes) ─────────────
ga=pd.read_csv(cfg_path('gene_annotation'),index_col=0); ga['ens']=ga.index.str.split('.').str[0]
pc=ga[ga.gene_type=='protein_coding']
L=pd.read_csv(cfg_path('relative_gene_lengths'),index_col=0).iloc[:,0]; L.index=L.index.str.split('.').str[0]
def c2t(C):
    l=L.reindex(C.index); C=C[l.notna()]; l=l[l.notna()]; rpk=C.div(l,axis=0); return np.log2(rpk/rpk.sum()*1e6+1)
def e2s(Mx):
    if GID=='entrez': return genemap.to_entrez(Mx,genemap.ens2entrez().entrez)
    m=pc.drop_duplicates('ens').set_index('ens').gene_name; Mx=Mx[Mx.index.isin(m.index)]; Mx.index=m.loc[Mx.index].values; return Mx.groupby(level=0).max()
def series_chars(f):
    d={}
    for line in gzip.open(f,'rt'):
        if line.startswith('!Sample_title'): d['title']=[x.strip('"\n') for x in line.split('\t')[1:]]
        if line.startswith('!Sample_characteristics_ch1'):
            v=[x.strip('"\n') for x in line.split('\t')[1:]]; k=v[0].split(':')[0]; d[k]=[x.split(': ',1)[-1] for x in v]
    return pd.DataFrame(d)
m1=series_chars(str(cfg_path('gse295399_series_matrix'))); m1['col']=m1.title.str.replace('Tumor_','')
C1=pd.read_csv(cfg_path('gse295399_counts'),sep='\t',index_col=0); C1.index=C1.index.str.split('.').str[0]
R1=e2s(c2t(C1)); m1=m1[m1.col.isin(R1.columns)].set_index('col')
keep=m1[m1.histotype.isin(['HGSC','LGSC','CCC','EC','MC'])]; R1=R1[keep.index]; y=keep.histotype
if GID=='entrez':
    R2=pd.read_csv(cfg_path('gse226870_matrix'),sep='\t').drop(columns='SYMBOL'); R2['ENSEMBL']=R2.ENSEMBL.astype(str).str.split('.').str[0]
    R2=e2s(R2.groupby('ENSEMBL').max())
else:
    R2=pd.read_csv(cfg_path('gse226870_matrix'),sep='\t').drop(columns='ENSEMBL').groupby('SYMBOL').max()
cmap=pd.read_csv(cfg_path('gse121103_colmap')); C3=pd.read_csv(cfg_path('gse121103_counts'),sep='\t').set_index('Geneid')
C3.index=C3.index.str.split('.').str[0]; C3=C3[list(cmap.count_col)].rename(columns=dict(zip(cmap.count_col,cmap.title.str.replace('_RNA-seq',''))))
R3=e2s(c2t(C3))
Xh=pd.read_parquet(HARM,columns=list(X.columns)).astype('float32')
# apply the same DepMap swap on the full harmonised gene space
for s in X.columns:
    if A.loc[s,'dataset'] in ('DepMap_24Q4','DepMap_NextGen') and A.loc[s,'model_unit'] in E_all.columns:
        Xh[s]=E_all[A.loc[s,'model_unit']].reindex(Xh.index).values
hg_genes=sorted(set(R1.index)&set(Xh.index)&set(R2.index)&set(R3.index))
R1g=R1.loc[hg_genes]; R1g=R1g[R1g.mean(axis=1)>1]
Fst=pd.Series(f_oneway(*[R1g.loc[:,y==c].values.T for c in sorted(y.unique())]).statistic,index=R1g.index).sort_values(ascending=False)
G=Fst.index[:300].tolist()
Cn=pd.concat({c:R1.loc[G,y.index[y==c]].mean(axis=1) for c in sorted(y.unique())},axis=1)
labmap={'OV_HGSC':'HGSC','OV_serous':'HGSC','OV_serous_NOS':'HGSC','OV_clear_cell':'CCC','OV_endometrioid':'EC','OV_mucinous':'MC','OV_LGSC':'LGSC'}
O=U[U.lineage=='ovary'].copy(); O['label']=O.histology.map(labmap).fillna('other/unknown')
V=pd.concat({i:aligned(G,Xh[r.samples]).mean(axis=1) for i,r in O.iterrows()},axis=1); Cn=aligned(G,Cn); assert_same_order(V,Cn)
def near(Vm,Cm_): return pd.DataFrame({c:Vm.corrwith(Cm_[c],method='spearman') for c in Cm_.columns}).idxmax(axis=1)
O['hist_uncentred']=near(V,Cn)
mu_c=Cn.mean(axis=1); O['hist_ordinary']=near(V.sub(mu_c,axis=0),Cn.sub(mu_c,axis=0))
res=[]
for ds,g in O.groupby('dataset'):
    gm=pd.concat([V[g.index[g.label==c]].mean(axis=1) for c in g.label.unique()],axis=1).mean(axis=1)
    res.append(near(V[g.index].sub(gm,axis=0),Cn.sub(mu_c,axis=0)))
O['hist_label_balanced']=pd.concat(res)
hist_summary={}
for col in ['hist_uncentred','hist_ordinary','hist_label_balanced']:
    nh=O[O.label.isin(['CCC','EC','MC','LGSC'])]; hh=O[O.label=='HGSC']
    hist_summary[col]=dict(nonHGSC=f'{int((nh.label==nh[col]).sum())}/{len(nh)}',HGSC=f"{int((hh[col]=='HGSC').sum())}/{len(hh)}")
p('Histotype:',hist_summary)
ppv=O.groupby('hist_uncentred').apply(lambda d: pd.Series(dict(true_positives=int((d.label==d.name).sum()),n_called=len(d)))).reset_index().rename(columns={'hist_uncentred':'called_as'})
ppv['PPV']=(ppv.true_positives/ppv.n_called).round(4); ppv.to_csv(OUT/'histotype_call_PPV.csv',index=False)
for col in ['hist_uncentred','hist_ordinary','hist_label_balanced']:
    pd.crosstab(O.label,O[col],margins=True,margins_name='All').to_csv(OUT/f'histotype_confusion_{col.split("_",1)[1]}.csv')
U=U.join(O[['label','hist_uncentred','hist_ordinary','hist_label_balanced']])

# ───────────── 10. markers (WT1/KRT20/CDX2/MLH1) for all units ─────────────
S2E=genemap.hgnc().dropna(subset=['entrez_id']).set_index('symbol').entrez_id
def gkey(sym_): return S2E[sym_] if GID=='entrez' else sym_
MARKER_KEYS={g_:gkey(g_) for g_ in ['WT1','KRT20','CDX2','MLH1','PAX8']}
p('marker gene keys:',MARKER_KEYS)
for gname in ['WT1','KRT20','CDX2','MLH1','PAX8']:
    gk=MARKER_KEYS[gname]
    U[gname+'_log2tpm']=[float(Xh.loc[gk,r.samples].mean()) if gk in Xh.index else np.nan for r in U.itertuples()]
U['WT1_status']=np.where(U.WT1_log2tpm>=3,'retained (>=3)','low (<3)')
mu_=U[U.histology=='OV_mucinous'].copy()
U['mucinous_GI_like']=np.where(U.histology=='OV_mucinous',((U.KRT20_log2tpm>=3)|(U.CDX2_log2tpm>=3))&(U.WT1_log2tpm<3)&U.nearest_A.isin(['STAD','COAD']),np.nan)
markers=dict(WT1_low_HGSC_serous=f"{int((U[hg].WT1_log2tpm<3).sum())}/{int(hg.sum())}",
             mucinous_GI_like=';'.join(U[U.mucinous_GI_like==True].model_unit),
             mucinous_GI_like_n=f"{int((U.mucinous_GI_like==True).sum())}/{int((U.histology=='OV_mucinous').sum())}",
             mucinous_other=U[(U.histology=='OV_mucinous')&(U.mucinous_GI_like!=True)][['model_unit','nearest_A','WT1_log2tpm','KRT20_log2tpm','CDX2_log2tpm']].round(2).to_dict('records'))
p('Markers:',markers)
# MLH1 threshold-crossing check (DNA-surrogate MMRd criterion uses MLH1 <= 3.095 for DepMap/NextGen lines)
X24=pd.read_parquet(HARM,columns=list(X.columns))
U['MLH1_log2tpm_24Q4']=[float(X24.loc[MARKER_KEYS['MLH1'],r.samples].mean()) for r in U.itertuples()]; del X24
cross=U[(U.lineage=='uterus')&U.dataset.isin(['DepMap_24Q4','DepMap_NextGen'])&((U.MLH1_log2tpm<=3.095)!=(U.MLH1_log2tpm_24Q4<=3.095))]
markers['MLH1_threshold_crossings_vs_24Q4']=';'.join(cross.model_unit) or 'none'
p('MLH1 threshold crossings vs 24Q4 (DepMap/NextGen uterine):',markers['MLH1_threshold_crossings_vs_24Q4'])


# ───────────── 11. save unit-level table + key numbers ─────────────
U.drop(columns=['samples']).to_csv(OUT/'units_215_calls.csv',index=False)
key=dict(tag=tag,gene_id=GID,feature_genes_absent_from_depmap=dm_absent_top,depmap_expr=expr_path.name,depmap_expr_md5=md5(expr_path),gdc_release=json.load(open(cfg_path('gdc_status_json')))['data_release'],
         reference_counts=proj.value_counts().sort_index().to_dict(),extra=extra,cptac=cpt,loo_A=looA.to_dict(),legacy_reproduction=leg,
         histotype=hist_summary,markers={k:v for k,v in markers.items() if k!='mucinous_other'},
         expr_source_counts={f'{k[0]}|{k[1]}':int(v) for k,v in sw.items()},depmap_qc=DEPMAP_QC,methodB_qc=METHODB_QC)
json.dump(key,open(OUT/'key_numbers.json','w'),indent=1,default=str)
open(OUT/'run_log.txt','w').write('\n'.join(LOG)+'\n')
p('wrote',OUT)
if not a.no_outputs:
    cmd=[sys.executable,str(PIPE/'make_outputs.py'),'--tag',tag]+(['--publish'] if a.publish else [])
    subprocess.run(cmd,check=True)
