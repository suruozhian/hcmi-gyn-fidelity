#!/usr/bin/env python3
"""NextGen-HCMI multi-omic identity check (Methods 2.3), gene-keyed by Entrez Gene ID (default) or symbol.
(i)  expression residual correlation (feature genes; gene mean across all models subtracted) vs every cohort model;
     threshold = 99th percentile of all non-NextGen model pairs.
(ii) ASCAT gene-level copy-number correlation vs HCMI models; threshold = maximum between different HCMI cases.
(iii) shared damaging/truncating genes (excluding TP53) vs HCMI models; evidence only when the count exceeds the observed
     maximum between different HCMI cases in the same cohort (an operational rule set from negative controls before the
     NextGen pairs were examined, not a pre-registered constant) and is specific to one HCMI case (Jaccard; other NextGen
     models); the specificity part has no numeric threshold and is judged from the reported inputs, as in the original check.
Usage: python3 identity_test.py --tag 26Q3 [--depmap-expr <26Q3 expression csv>] [--gene-id entrez|symbol]
"""
import sys, argparse, json, re, gzip
from pathlib import Path
import numpy as np, pandas as pd
PIPE=Path(__file__).resolve().parent; sys.path.insert(0,str(PIPE))
from paths import path as cfg_path
import genemap
from fidelity_core import aligned, assert_same_order
ap=argparse.ArgumentParser(); ap.add_argument('--tag',required=True); ap.add_argument('--depmap-expr',default=None)
ap.add_argument('--gene-id',choices=['entrez','symbol'],default='entrez'); ap.add_argument('--features',default=None)
ap.add_argument('--out',default=None)
a=ap.parse_args(); GID=a.gene_id
OUT=Path(a.out) if a.out else PIPE/'results'/a.tag/'identity'; OUT.mkdir(parents=True,exist_ok=True)
LOG=[]
def p(*x): s=' '.join(map(str,x)); print(s,flush=True); LOG.append(s)
OFFICIAL={'ACH-003350':'HCM-WCMC-0673-C56','ACH-003723':'HCM-WCMC-0673-C56','ACH-003709':'HCM-CSHL-0806-C54','ACH-003710':'HCM-CSHL-0810-C54','ACH-003711':'HCM-CSHL-0811-C55'}
TP53={'entrez':'7157','symbol':'TP53'}[GID]
ND=cfg_path('nextgen_dir'); P1=cfg_path('hcmi_processed_dir')
md=pd.read_csv(ND/'model_metadata.csv'); NG=md[md.OncotreeLineage.isin(['Ovary/Fallopian Tube','Uterus'])&(md.IsNextGen==True)].set_index('ModelID'); ids=list(NG.index)
def hdr_keys(cols):
    return genemap.depmap_entrez(cols) if GID=='entrez' else [str(c).split(' (')[0] for c in cols]
# ---------- (i) expression
A=pd.read_csv(cfg_path('sample_annotations')); Am=A[A.kind=='model']
HARM=cfg_path('reference_dir')/'expr_entrez_log2tpm_models.parquet' if GID=='entrez' else cfg_path('harmonised_matrix')
X=pd.read_parquet(HARM,columns=list(Am['sample']))
src='DepMap 24Q4 (2D lines) + NextGen manuscript file set (organoids)'
if a.depmap_expr:
    E=pd.read_csv(a.depmap_expr,index_col=0,low_memory=False)
    if 'ModelID' in E.columns: E=E[E.IsDefaultEntryForModel.astype(str)=='Yes'].set_index('ModelID')
    gc=[c for c in E.columns if re.search(r'\(\d+\)$',str(c))]; E=E[gc].T; E.index=hdr_keys(E.index); E=E[pd.notna(E.index)].groupby(level=0).max()
    n=0
    for s in X.columns:
        d=Am.set_index('sample').loc[s,'dataset']; m=s.split('|')[1]
        if d in ('DepMap_24Q4','DepMap_NextGen') and m in E.columns: X[s]=E[m].reindex(X.index).values; n+=1
    src=f'{Path(a.depmap_expr).name} for DepMap 2D lines and NextGen organoids ({n} profiles swapped)'
p('expression source:',src,'| gene key:',GID)
feat=pd.read_csv(a.features or PIPE/'results'/a.tag/'feature_genes_top2000.csv').iloc[:,0].astype(str).tolist()
ngc=[c for c in X.columns if c.startswith('NextGen|')]; other=[c for c in X.columns if c not in ngc]
genes=[g for g in feat if g in X.index and X.loc[g].notna().all()]
p(f'feature genes {len(feat)}; present without missing values in all {X.shape[1]} model profiles: {len(genes)}')
M=aligned(genes,X); R=M.sub(M.mean(axis=1),axis=0); C=R.corr()
anno=Am.set_index('sample')
rows=[]
for c in ngc:
    s=C.loc[c,other].sort_values(ascending=False); z=(s.iloc[0]-s.iloc[1:].mean())/s.iloc[1:].std(); b=anno.loc[s.index[0]]
    rows.append(dict(nextgen=c.split('|')[1],best_sample=s.index[0],best_dataset=b.dataset,best_unit=b.model_unit,best_r=round(s.iloc[0],3),second_r=round(s.iloc[1],3),z_best=round(z,2)))
I1=pd.DataFrame(rows)
ngC=C.loc[ngc,ngc].where(~np.eye(len(ngc),dtype=bool))
I1['best_other_nextgen_r']=ngC.max(axis=1).round(3).values; I1['best_other_nextgen']=ngC.idxmax(axis=1).str.split('|').str[1].values
allpairs=C.loc[other,other].values[np.triu_indices(len(other),1)]
thr_expr=float(np.percentile(allpairs,99))
hc=Am[Am.dataset=='HCMI']; same=[]
for case,g in hc.groupby('patient'):
    s_=list(g['sample']); same+=[C.loc[s_[i],s_[j]] for i in range(len(s_)) for j in range(i+1,len(s_))]
p(f'expression residual r: HCMI same-case model pairs median {np.median(same):.3f} (min {np.min(same):.3f}, n={len(same)}); all non-NextGen model pairs 99th pct {thr_expr:.3f}')
# ---------- (ii) copy number
ga=pd.read_csv(cfg_path('gene_annotation'),index_col=0)
hcn=pd.read_parquet(P1/'hcmi_cn_gene.parquet'); hs=pd.read_csv(P1/'hcmi_cn_samples.csv').set_index('col')
hcn=hcn[[c for c in hcn.columns if hs.loc[c,'role']=='model']]
if GID=='entrez':
    em=genemap.ens2entrez(); k=pd.Series(hcn.index.str.split('.').str[0],index=hcn.index).map(em.entrez)
else:
    k=pd.Series(ga.reindex(hcn.index).gene_name.values,index=hcn.index)
hcn=hcn[k.notna().values]; hcn.index=k.dropna().values; hcn=hcn.groupby(level=0).median()
hl=np.log2((hcn+0.5).div(hcn.median(axis=0)+0.5,axis=1))
CN=pd.read_csv(ND/'next_gen_copy_number.csv',index_col=0)
nc=CN.loc[[i for i in ids if i in CN.index]].T; nc.index=hdr_keys(nc.index); nc=nc[pd.notna(nc.index)].groupby(level=0).mean(); nl=np.log2(nc.clip(lower=0.05))
g=sorted(set(hl.index)&set(nl.index)); g=[x for x in g if hl.loc[x].notna().all() and nl.loc[x].notna().all()]
hlg,nlg=aligned(g,hl,nl)
cc=pd.DataFrame({h:nlg.corrwith(hlg[h]) for h in hlg.columns})
I2=pd.DataFrame({'nextgen':cc.index,'best_hcmi_cn':cc.idxmax(axis=1).values,'best_cn_r':cc.max(axis=1).round(3).values,'second_cn_r':cc.apply(lambda r:r.sort_values().iloc[-2],axis=1).round(3).values})
hcase=hs.loc[hlg.columns,'case']; hcc=hlg.corr()
pos=[hcc.iloc[i,j] for i in range(len(hcc)) for j in range(i+1,len(hcc)) if hcase.iloc[i]==hcase.iloc[j]]
neg=[hcc.iloc[i,j] for i in range(len(hcc)) for j in range(i+1,len(hcc)) if hcase.iloc[i]!=hcase.iloc[j]]
thr_cn=float(np.max(neg))
p(f'CN genes compared: {len(g)}; HCMI same-case model pairs r median {np.median(pos):.3f} (n={len(pos)}); different-case max {thr_cn:.3f}')
# ---------- (iii) truncating mutations
maf=pd.read_csv(P1/'hcmi_maf_combined.csv.gz',low_memory=False); msm=pd.read_csv(P1/'hcmi_maf_samples.csv')
maf=maf[maf['sample'].isin(msm[msm.role=='model']['sample'])]
tr=maf[maf.Variant_Classification.isin(['Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site','Nonstop_Mutation'])].copy()
if GID=='entrez':   # Ensembl gene of each variant from the raw GDC MAF -> Entrez (rules E1/E2)
    fd=cfg_path('hcmi_files_dir'); raw=[]
    for fid in tr.file_id.unique():
        f=[x for x in fd.glob(f'{fid}__*.maf.gz')][0]
        r=pd.read_csv(f,sep='\t',comment='#',usecols=['Hugo_Symbol','Chromosome','Start_Position','Variant_Classification','Gene'],low_memory=False); r['file_id']=fid; raw.append(r)
    raw=pd.concat(raw).drop_duplicates(['file_id','Chromosome','Start_Position','Hugo_Symbol','Variant_Classification'])
    tr=tr.merge(raw,on=['file_id','Chromosome','Start_Position','Hugo_Symbol','Variant_Classification'],how='left')
    tr['key']=tr.Gene.astype(str).str.split('.').str[0].map(em.entrez)
    p(f'truncating HCMI variants {len(tr)}; with Entrez via Ensembl {int(tr.key.notna().sum())}')
else: tr['key']=tr.Hugo_Symbol
hset=tr.dropna(subset=['key']).groupby('sample').key.apply(lambda s:set(s)-{TP53})
DM=pd.read_csv(ND/'next_gen_damaging.csv',index_col=0)
dmg=DM.loc[[i for i in ids if i in DM.index]]; dmg.columns=hdr_keys(dmg.columns); dmg=dmg.loc[:,pd.notna(dmg.columns)]
nset={i:set(dmg.columns[dmg.loc[i]>0])-{TP53} for i in dmg.index}
hk=list(hset.index); case_of=msm.set_index('sample').case
neg3=[len(hset[x]&hset[y]) for n_,x in enumerate(hk) for y in hk[n_+1:] if case_of[x]!=case_of[y]]
thr_mut=int(max(neg3))
r3=[]
for i,s in nset.items():
    sh={h:len(s&t) for h,t in hset.items()}; b=max(sh,key=sh.get); j=len(s&hset[b])/max(len(s|hset[b]),1)
    others=max([len(nset[o]&hset[b]) for o in nset if o!=i] or [0])
    r3.append(dict(nextgen=i,n_damaging_genes=len(s),best_hcmi_mut=b,shared_damaging_genes=sh[b],jaccard=round(j,3),max_shared_same_hcmi_other_nextgen=others))
I3=pd.DataFrame(r3)
p(f'HCMI truncating genes: different-case max shared {thr_mut}')
I=I1.merge(I2,on='nextgen',how='outer').merge(I3,on='nextgen',how='outer')
I['OncotreeCode']=I.nextgen.map(NG.OncotreeCode)
I['crit_i']=I.best_r>thr_expr; I['crit_ii']=I.best_cn_r>thr_cn; I['crit_iii_count']=I.shared_damaging_genes>thr_mut
def case_of_cn(s): return hs.loc[s,'case'] if isinstance(s,str) and s in hs.index else None
I['case_i']=np.where(I.crit_i&(I.best_dataset=='HCMI'),I.best_unit,None)
I['case_ii']=[case_of_cn(s) if ok else None for s,ok in zip(I.best_hcmi_cn,I.crit_ii)]
I['case_iii']=[case_of.get(s) if ok else None for s,ok in zip(I.best_hcmi_mut,I.crit_iii_count)]
I['official_hcmi_case']=I.nextgen.map(OFFICIAL)
I.to_csv(OUT/'nextgen_identity_check.csv',index=False)
summ=dict(tag=a.tag,gene_id=GID,expression_source=src,feature_genes=len(feat),feature_genes_used=len(genes),thr_expr_99pct=round(thr_expr,3),
          thr_cn_max_diff_case=round(thr_cn,3),thr_mut_max_diff_case=thr_mut,cn_genes=len(g))
pd.set_option('display.width',250)
p(I[['nextgen','best_dataset','best_unit','best_r','crit_i','best_hcmi_cn','best_cn_r','crit_ii','best_hcmi_mut','shared_damaging_genes','jaccard','max_shared_same_hcmi_other_nextgen','crit_iii_count','official_hcmi_case']].to_string(index=False))
json.dump(summ,open(OUT/'identity_summary.json','w'),indent=1)
open(OUT/'nextgen_identity_summary.txt','w').write('\n'.join(LOG)+'\n'+json.dumps(summ,indent=1)+'\n')
