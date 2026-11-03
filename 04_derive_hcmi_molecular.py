#!/usr/bin/env python3
"""HCMI molecular derivations from GDC open-access files (Methods 2.7; Fig. 4; Supplementary Fig. S7; HCMI rows of the
DNA surrogate classes in Table 3, Table S2 and Table S10).

Inputs : <hcmi_case_diagnoses>  GDC cases-API extract (case, primary_site, primary_diagnosis; fetch_public_data.py)
         inputs/manifests/hcmi_{rna,cn,maf,meth}_manifest.csv and <hcmi_processed_dir> (written by 00e and 00h)
         HCMI SeSAMe level-3 beta files (60 files; GDC open access) in tcga_file_dirs
         <epic_manifest>        EPIC hg38 manifest, GENCODE v36 annotation (Zhou lab InfiniumAnnotationV1)
         TCGA subsample STAR files (inputs/manifests/tcga_gyn_subsample_manifest.csv; 777 tumours) and <gene_annotation>
Outputs in <molecular_dir>:
         hcmi_case_labels.csv               histology group per HCMI gynaecological case (rule grp() below; leiomyosarcoma excluded)
         hcmi_meth_beta.parquet             beta values, one column per file ('<sample>|<file_id[:8]>')
         mlh1_promoter_probes.csv           EPIC probes annotated to MLH1 with distance to TSS -1000..+500 bp
         hcmi_sample_molecular_table.csv    per-sample DNA features and endometrial surrogate class (uterine cases)
         rna_paired_purityfilt.csv          model vs own / other tumours, Spearman rho over the 2,000 purity-filtered feature genes
         paired_mutation_concordance.csv    model vs own tumour, shared non-synonymous variants (Hugo:chromosome:start)
         paired_methylation_concordance.csv model vs own / other tumours, Pearson r over the 10,000 most variable autosomal CpGs
Surrogate class (WHO hierarchy POLE > MMRd > p53abn > NSMP; not IHC): POLE exonuclease-domain hotspot -> POLEmut; MLH1
promoter beta >= 0.3 or truncating MLH1/MSH2/MSH6/PMS2 variant -> MMRd; any non-synonymous TP53 variant -> p53abn; otherwise NSMP
('incomplete (no MAF)' without mutation data; 'no data' without mutation and methylation data)."""
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path, paths as cfg_paths
PIPE=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(); ap.add_argument('--gdc-dir',default=None); a=ap.parse_args()
dirs=[Path(a.gdc_dir)] if a.gdc_dir else cfg_paths('tcga_file_dirs')
def find(name):
    for d in dirs:
        if (d/name).exists(): return d/name
    raise FileNotFoundError(name)
HP=cfg_path('hcmi_processed_dir'); OUT=cfg_path('molecular_dir'); OUT.mkdir(parents=True,exist_ok=True)
MAN=PIPE/'inputs/manifests'

# ---- 1. case histology groups (GDC primary diagnosis)
dx=pd.read_csv(cfg_path('hcmi_case_diagnoses'))
def grp(site,d):
    d=d.lower()
    if site=='Ovary':
        if 'borderline' in d: return 'OV_borderline'
        if 'mucinous' in d: return 'OV_mucinous'
        if 'clear cell' in d: return 'OV_clear_cell'
        if 'endometrioid' in d: return 'OV_endometrioid'
        if 'mullerian' in d: return 'OV_mullerian_mixed'
        if 'serous' in d: return 'OV_serous'
        return 'OV_other'
    if 'carcinosarcoma' in d: return 'UT_carcinosarcoma'
    if 'serous' in d: return 'UT_serous'
    if 'endometrioid' in d: return 'UT_endometrioid'
    if 'mixed' in d: return 'UT_mixed'
    if 'undifferentiated' in d: return 'UT_undifferentiated'
    return 'UT_other'
dx['histology']=[grp(s,d) for s,d in zip(dx.primary_site,dx.primary_diagnosis)]
dx['lineage']=np.where(dx.primary_site=='Ovary','ovary','uterus')
dx=dx[~dx.primary_diagnosis.str.contains('Leiomyosarcoma')]
dx.to_csv(OUT/'hcmi_case_labels.csv',index=False)
lab=dx.set_index('case')
for k in ['rna','cn','maf','meth']:   # the shipped file manifests carry the same labels
    mf=pd.read_csv(MAN/f'hcmi_{k}_manifest.csv'); j=mf.join(lab[['histology','lineage']],on='case',rsuffix='_gdc')
    assert (j.histology==j.histology_gdc).all() and (j.lineage==j.lineage_gdc).all(), f'hcmi_{k}_manifest labels differ from GDC'

# ---- 2. methylation beta matrix (SeSAMe level-3 betas)
bs=pd.read_csv(MAN/'hcmi_meth_manifest.csv')
B=pd.DataFrame({x.col:pd.read_csv(find(f'{x.file_id}__{x.file_name}'),sep='\t',header=None,index_col=0).iloc[:,0].astype('float32') for x in bs.itertuples()})
B.index.name='probeID'
B.to_parquet(OUT/'hcmi_meth_beta.parquet'); bs=bs.set_index('col')

# ---- 3. paired RNA concordance (purity-filtered 2,000 most variable TCGA genes; Spearman)
H=pd.read_parquet(HP/'hcmi_rna_tpm.parquet'); hs=pd.read_csv(HP/'hcmi_rna_samples.csv').set_index('col')
tm=pd.read_csv(MAN/'tcga_gyn_subsample_manifest.csv'); T={}
for x in tm.itertuples():
    t=pd.read_csv(find(f'{x.file_id}__{x.file_name}'),sep='\t',comment='#',usecols=['gene_id','tpm_unstranded'])
    T[x.case]=t[t.gene_id.str.startswith('ENSG')].set_index('gene_id').tpm_unstranded.astype('float32')
T=pd.DataFrame(T)
ga=pd.read_csv(cfg_path('gene_annotation'),index_col=0); sym=ga.gene_name
pc=ga.index[ga.gene_type=='protein_coding']; common=pc.intersection(H.index).intersection(T.index)
LH=np.log2(H.loc[common]+1); LT=np.log2(T.loc[common].astype(float)+1)
def score(L,genes):
    ids=[g for g in sym.index[sym.isin(genes)] if g in L.index]
    z=L.loc[ids]; z=(z.sub(z.mean(axis=1),axis=0)).div(z.std(axis=1)+1e-9,axis=0); return z.mean()
sS=score(LT,['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRB','FAP','THY1','COL6A3','SPARC'])
sI=score(LT,['PTPRC','CD3E','CD2','CD14','CD68','LYZ','CD74','HLA-DRA','MS4A1','CD53'])
expr=LT.index[(LT.mean(axis=1)>1)]
R=LT.loc[expr].T.apply(lambda c: c.rank())
rS=R.corrwith(pd.Series(rankdata(sS),index=sS.index)); rI=R.corrwith(pd.Series(rankdata(sI),index=sI.index))
purity=set(expr[(rS.values>0.4)|(rI.values>0.4)])
top_f=LT.loc[[g for g in expr if g not in purity]].var(axis=1).sort_values(ascending=False).index[:2000]
Bx=LH.loc[top_f]; out=hs[['case','sample','sample_type','role','lineage','histology']]
mo=out[out.role=='model']; tum=out[out.role=='tumour']; P=[]
for s,r in mo.iterrows():
    own=tum[tum.case==r.case]
    if len(own)==0: continue
    cors=pd.Series({t:spearmanr(Bx[s],Bx[t])[0] for t in tum.index})
    for t in own.index:
        P.append(dict(model=s,case=r.case,histology=r.histology,tumour=t,tumour_type=tum.loc[t,'sample_type'],rho_own=cors[t],
            rho_other_median=cors.drop(own.index).median(),rho_other_max=cors.drop(own.index).max(),
            rank_own=int((cors>cors[t]).sum())+1,n_tumours=len(cors)))
pd.DataFrame(P).to_csv(OUT/'rna_paired_purityfilt.csv',index=False)

# ---- 4. per-sample DNA features and endometrial surrogate class
NONSYN=['Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Splice_Site','Nonstop_Mutation','Translation_Start_Site']
TRUNC=['Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site','Nonstop_Mutation']
POLE_HOT={'p.P286R','p.V411L','p.S297F','p.A456P','p.S459F','p.F367S','p.L424I','p.M295R','p.P436R','p.M444K','p.D368Y'}
M=pd.read_csv(HP/'hcmi_maf_combined.csv.gz',low_memory=False); ms=pd.read_csv(HP/'hcmi_maf_samples.csv')
M=M[M.Variant_Classification.isin(NONSYN)].copy()
M['key']=M.Hugo_Symbol+':'+M.Chromosome.astype(str)+':'+M.Start_Position.astype(str)
feat=[]
for s,g in M.groupby('sample'):
    pole=g[g.Hugo_Symbol=='POLE']
    mmr=g[g.Hugo_Symbol.isin(['MLH1','MSH2','MSH6','PMS2'])&g.Variant_Classification.isin(TRUNC)]
    tp53=g[g.Hugo_Symbol=='TP53']
    feat.append(dict(sample=s,n_nonsyn=len(g),POLE_hotspot=';'.join(sorted(set(pole.HGVSp_Short)&POLE_HOT)),
        MMR_truncating=';'.join(mmr.Hugo_Symbol+' '+mmr.HGVSp_Short.fillna('')),TP53=';'.join(tp53.HGVSp_Short.fillna(tp53.Variant_Classification))))
F=pd.DataFrame(feat).merge(ms,on='sample',how='right')
cn=pd.read_parquet(HP/'hcmi_cn_gene.parquet'); cs=pd.read_csv(HP/'hcmi_cn_samples.csv').set_index('col'); cnb={}
for c in cn.columns:
    v=cn[c].dropna(); mode=v.mode().iloc[0]; cnb[cs.loc[c,'sample']]=dict(cn_frac_altered=float((v!=mode).mean()))
CNB=pd.DataFrame(cnb).T; CNB.index.name='sample'
ann=pd.read_csv(cfg_path('epic_manifest'),sep='\t',usecols=['CpG_chrm','probeID','geneNames','distToTSS'],low_memory=False)
a_=ann[ann.geneNames.fillna('').str.contains(r'(?:^|;)MLH1(?:;|$)')].copy()
def mind(r):
    dd=[int(d) for g,d in zip(r.geneNames.split(';'),str(r.distToTSS).split(';')) if g=='MLH1' and d not in ('NA','nan')]
    return min(dd,key=abs) if dd else np.nan
a_['d']=a_.apply(mind,axis=1)
prom=[p for p in a_[(a_.d>=-1000)&(a_.d<=500)].probeID if p in B.index]
pd.Series(prom).to_csv(OUT/'mlh1_promoter_probes.csv',index=False)
mlh=B.loc[prom].mean().rename('MLH1_prom_beta'); mlh.index=bs.loc[mlh.index,'sample']
allsamples=pd.concat([ms[['sample','case','role']],bs[['sample','case','role']],cs[['sample','case','role']],out[['sample','case','role']]]).drop_duplicates('sample')
S=allsamples.merge(F.drop(columns=['case','role','lineage','histology','file_id']),on='sample',how='left').merge(CNB,left_on='sample',right_index=True,how='left')\
    .merge(mlh.groupby(level=0).mean(),left_on='sample',right_index=True,how='left')
S=S.join(lab[['primary_site','primary_diagnosis','histology','lineage']],on='case')
def ecclass(r):
    if pd.isna(r.n_nonsyn) and pd.isna(r.MLH1_prom_beta): return 'no data'
    if isinstance(r.POLE_hotspot,str) and r.POLE_hotspot: return 'POLEmut (surrogate)'
    if (r.MLH1_prom_beta>=0.3) or (isinstance(r.MMR_truncating,str) and r.MMR_truncating): return 'MMRd (surrogate)'
    if isinstance(r.TP53,str) and r.TP53: return 'p53abn (surrogate)'
    if pd.isna(r.n_nonsyn): return 'incomplete (no MAF)'
    return 'NSMP (surrogate)'
S['EC_class']=np.where(S.lineage=='uterus',S.apply(ecclass,axis=1),'')
S=S.sort_values(['lineage','histology','case','role','sample'])
S.to_csv(OUT/'hcmi_sample_molecular_table.csv',index=False)

# ---- 5. paired mutation concordance
Q=[]
for case,g in S.groupby('case'):
    for m_ in g[g.role=='model']['sample']:
        for t_ in g[g.role=='tumour']['sample']:
            x=set(M[M['sample']==m_].key); y=set(M[M['sample']==t_].key)
            if x and y: Q.append(dict(case=case,model=m_,tumour=t_,kind='mutation',n_model=len(x),n_tumour=len(y),shared=len(x&y),jaccard=len(x&y)/len(x|y),tumour_muts_in_model=len(x&y)/len(y)))
pd.DataFrame(Q).join(lab[['histology']],on='case').to_csv(OUT/'paired_mutation_concordance.csv',index=False)

# ---- 6. paired methylation concordance (autosomal probes without missing values; 10,000 most variable; Pearson)
chrom=ann.set_index('probeID').CpG_chrm
auto=[p for p in B.index if chrom.get(p) not in ('chrX','chrY',None) and isinstance(chrom.get(p),str)]
Bm=B.loc[auto].dropna(); Bm.columns=bs.loc[Bm.columns,'sample']; Bm=Bm.T.groupby(level=0).mean().T
C=Bm.loc[Bm.var(axis=1).sort_values(ascending=False).index[:10000]].corr()
alltum=[s for s in S[S.role=='tumour']['sample'] if s in C]; Q=[]
for case,g in S.groupby('case'):
    mods=[s for s in g[g.role=='model']['sample'] if s in C]; tums=[s for s in g[g.role=='tumour']['sample'] if s in C]
    for m_ in mods:
        for t_ in tums:
            oth=C.loc[m_,[t for t in alltum if t not in tums]]
            Q.append(dict(case=case,model=m_,tumour=t_,r_own=C.loc[m_,t_],r_other_median=oth.median(),r_other_max=oth.max(),rank_own=int((C.loc[m_,alltum]>C.loc[m_,t_]).sum())+1,n_tumours=len(alltum)))
pd.DataFrame(Q).join(lab[['histology']],on='case').to_csv(OUT/'paired_methylation_concordance.csv',index=False)
u=S[(S.lineage=='uterus')&(S.role=='model')]
print('HCMI cases',len(dx),'; methylation files',B.shape[1],'; MLH1 promoter probes',len(prom),'; RNA pairs',len(P),'; uterine model surrogate classes',
      u.drop_duplicates('case').EC_class.value_counts().to_dict())
