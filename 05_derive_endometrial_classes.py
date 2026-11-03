#!/usr/bin/env python3
"""Endometrial DNA surrogate classes (Methods 2.7; Fig. 3, Table 3, Supplementary Tables S9 and S10, DNA columns of Table S2).

1. CPTAC UCEC evaluation (Table S9, Fig. 3a): the surrogate rule applied to cBioPortal ucec_cptac_2020 mutations and MLH1
   methylation, compared with the published genomic subtype (POLE -> POLEmut, MSI-H -> MMRd, CNV_high -> p53abn, CNV_low -> NSMP).
   MLH1 expression threshold: Youden index of MLH1 log2(TPM+1) (GDC CPTAC STAR, harmonised matrix) for MLH1-methylated MSI-H
   tumours (in-sample; used for the expression proxy and the DepMap/NextGen MMRd criterion).
2. Classes of the independent uterine units (inputs/units_215.csv) plus the DNA-only NextGen uterine organoids
   (NextGen model_metadata, not the same patient as an HCMI model per the official NextGen-HCMI pairs in the crosswalk):
   HCMI: surrogate class of the first model of the case (04_derive_hcmi_molecular.py; DNA + MLH1 promoter methylation);
   DepMap 2D (DepMap 24Q4 OmicsSomaticMutations / OmicsSignatures): POLE hotspot -> POLEmut; MSIScore >= 10 or MLH1
   log2(TPM+1) <= threshold -> MMRd; non-synonymous TP53 -> p53abn; else NSMP ('no mutation data' without mutation records);
   NextGen organoids (NextGen file set hotspot/damaging matrices and MSIScore): same hierarchy (MMR: MSIScore >= 10 or MLH1 low);
   GSE307568 and E-MTAB-10664 organoids (no DNA data): MLH1 expression proxy only (unit mean MLH1 <= threshold: MLH1-low,
   else MLH1-retained); this is not a molecular class.
Outputs in <molecular_dir>: cptac_rule_evaluation.csv, uterine_units_classes.csv (dataset, model_unit, histology, cls, organoid),
mlh1_threshold.json. Unit histology is checked (the script stops otherwise) against the sample annotations of every
model profile of the unit, against inputs/units_215.csv and, for HCMI units, against the GDC case labels."""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path
PIPE=Path(__file__).resolve().parent; OUT=cfg_path('molecular_dir'); OUT.mkdir(parents=True,exist_ok=True)
HOT={'P286R','V411L','S297F','A456P','S459F','F367S','L424I','M295R','P436R','M444K','D368Y'}
MMR=['MLH1','MSH2','MSH6','PMS2']; TRUNC={'Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site','Nonstop_Mutation'}
def cls(pole,mmr,p53): return 'POLEmut' if pole else 'MMRd' if mmr else 'p53abn' if p53 else 'NSMP'
X=pd.read_parquet(cfg_path('harmonised_matrix'))   # symbol-keyed model matrix (00g); MLH1 row
# ---- 1. CPTAC UCEC
cl=pd.read_csv(cfg_path('cptac_ucec_genomic_subtype')); mu=pd.read_csv(cfg_path('cptac_ucec_mutations')); me=pd.read_csv(cfg_path('cptac_ucec_mlh1_methylation'))
truth=cl.set_index('sampleId').GENOMICS_SUBTYPE.map({'POLE':'POLEmut','MSI-H':'MMRd','CNV_high':'p53abn','CNV_low':'NSMP'})
mutd=set(mu.sampleId); meth=me.set_index('sampleId').MLH1_meth; rows=[]
for s in truth.index:
    if s not in mutd or s not in meth.index: continue
    g=mu[mu.sampleId==s]
    rows.append(dict(sample=s,truth=truth[s],POLE_hot=bool(set(g[g.gene=='POLE'].proteinChange)&HOT),MMR_trunc=bool(len(g[g.gene.isin(MMR)&g.mutationType.isin(TRUNC)])),
        MLH1_meth=meth[s],TP53_mut=bool(len(g[g.gene=='TP53'])),MLH1_expr=X.loc['MLH1','CPTAC|'+s] if 'CPTAC|'+s in X.columns else np.nan,n_mut=len(g)))
E=pd.DataFrame(rows)
E['rule_meth']=[cls(r.POLE_hot,(r.MLH1_meth>=0.3) or r.MMR_trunc,r.TP53_mut) for r in E.itertuples()]
ee=E.dropna(subset=['MLH1_expr']); y=(ee.truth=='MMRd')&(ee.MLH1_meth>=0.3)
best=max(np.sort(ee.MLH1_expr.unique()),key=lambda t: ((ee.MLH1_expr<=t)&y).sum()/max(y.sum(),1)-((ee.MLH1_expr<=t)&~y).sum()/max((~y).sum(),1))
E['rule_expr']=[cls(r.POLE_hot,(r.MLH1_expr<=best) or r.MMR_trunc,r.TP53_mut) if not np.isnan(r.MLH1_expr) else None for r in E.itertuples()]
E['rule_noDNA']=[('MMRd' if r.MLH1_expr<=best else 'not-MMRd(MLH1 expr)') if not np.isnan(r.MLH1_expr) else None for r in E.itertuples()]
E.to_csv(OUT/'cptac_rule_evaluation.csv',index=False)
json.dump(dict(mlh1_log2tpm_threshold=float(best),mlh1_log2tpm_threshold_rounded=round(float(best),3),n_cptac=len(E),n_cptac_rna=len(ee),
               rule_meth_concordant=int((E.rule_meth==E.truth).sum())),open(OUT/'mlh1_threshold.json','w'),indent=1)
# ---- 2. uterine units
U=pd.read_csv(PIPE/'inputs/units_215.csv'); A=pd.read_csv(cfg_path('sample_annotations'))
u=U[U.lineage=='uterus'][['dataset','model_unit','histology','is_org']].rename(columns={'is_org':'organoid'}).copy()
# unit histology from the per-sample annotations of the unit's model profiles; must be unique per unit and equal to units_215
# (and therefore to Table S2 / S10, checked again in make_outputs.py); E-MTAB-10664 units are single organoid lines (e.g. OEC-04-CC)
hs_=A[A.kind=='model'].groupby(['dataset','model_unit']).histology.agg(lambda x: ';'.join(sorted(set(x))))
u['histology_annotation']=[hs_.get((r.dataset,r.model_unit)) for r in u.itertuples()]
bad=u[u.histology_annotation!=u.histology]
if len(bad): raise SystemExit('unit histology differs between sample annotations and inputs/units_215.csv:\n'+bad.to_string())
lab_=pd.read_csv(OUT/'hcmi_case_labels.csv').set_index('case').histology   # HCMI: GDC primary diagnosis (04_derive_hcmi_molecular.py)
hb=u[(u.dataset=='HCMI')&(u.model_unit.map(lab_)!=u.histology)]
if len(hb): raise SystemExit('HCMI unit histology differs from the GDC case labels:\n'+hb.to_string())
u=u.drop(columns='histology_annotation')
H=pd.read_csv(OUT/'hcmi_sample_molecular_table.csv'); h1=H[H.role=='model'].drop_duplicates('case').set_index('case')
# DepMap 2D lines (DepMap 24Q4 mutation calls and MSIScore)
dmu=u[u.dataset=='DepMap_24Q4'].model_unit.tolist()
M=pd.read_csv(cfg_path('depmap24q4_somatic_mutations'),usecols=['ModelID','HugoSymbol','ProteinChange','VariantInfo','LikelyLoF'],low_memory=False)
M=M[M.ModelID.isin(dmu)]; sig=pd.read_csv(cfg_path('depmap24q4_signatures'),index_col=0)
NONSYN=M.VariantInfo.astype(str).str.contains('missense|stop_gained|frameshift|splice_acceptor|splice_donor|inframe|start_lost|stop_lost')
dm={}
for mid in dmu:
    if mid not in set(M.ModelID): dm[mid]='no mutation data'; continue
    g=M[(M.ModelID==mid)&NONSYN]; msi=sig.MSIScore.get(mid,np.nan)
    dm[mid]=cls(bool(set(g[g.HugoSymbol=='POLE'].ProteinChange.astype(str).str.replace('p.',''))&HOT),
                (msi>=10 if not np.isnan(msi) else False) or (X.loc['MLH1','DepMap|'+mid]<=best),bool(len(g[g.HugoSymbol=='TP53'])))
# NextGen organoids
ND=cfg_path('nextgen_dir'); md=pd.read_csv(ND/'model_metadata.csv').set_index('ModelID')
XW=pd.read_csv(PIPE/'inputs/crosswalk/crosswalk_hcmi_ach_sidm.csv'); same=set(XW[XW.official_hcmi_from_nextgen_supp.notna()].ACH)
HS=pd.read_csv(ND/'next_gen_hotspot.csv',index_col=0); DM=pd.read_csv(ND/'next_gen_damaging.csv',index_col=0); SG=pd.read_csv(ND/'next_gen_omics_signatures.csv',index_col=0)
for T_ in (HS,DM): T_.columns=T_.columns.str.split(' \\(').str[0]
ng={}
ngu=md[(md.OncotreeLineage=='Uterus')&(md.IsNextGen==True)&~md.index.isin(same)].index
for i in ngu:
    mlh1=X.loc['MLH1','NextGen|'+i] if 'NextGen|'+i in X.columns else np.nan; msi=SG.MSIScore.get(i,np.nan)
    ng[i]=cls(HS.loc[i,'POLE']>0,msi>=10 or (not np.isnan(mlh1) and mlh1<=best),(HS.loc[i,'TP53']>0) or (DM.loc[i,'TP53']>0))
HMAP={'UCS':'UT_carcinosarcoma','UEC':'UT_endometrioid','USC':'UT_serous','UMEC':'UT_mixed'}
dna_only=[dict(dataset='DepMap_NextGen',model_unit=i,histology=HMAP[md.loc[i,'OncotreeCode']],organoid=True) for i in ngu if i not in set(u.model_unit)]
u=pd.concat([u,pd.DataFrame(dna_only)],ignore_index=True)
# expression-only resources: unit mean MLH1 (samples of the unit in the harmonised matrix)
A=A[A.dataset.isin(['GSE307568','Berg2021_EMTAB10664'])].copy(); A['MLH1']=X.loc['MLH1',A['sample']].values
mlh1u=A.groupby(['dataset','model_unit']).MLH1.mean()
def unit_cls(r):
    if r.dataset=='HCMI': return str(h1.EC_class.get(r.model_unit,'NA')).replace(' (surrogate)','')
    if r.dataset=='DepMap_24Q4': return dm[r.model_unit]
    if r.dataset=='DepMap_NextGen': return ng[r.model_unit]
    return 'MLH1-low' if mlh1u[(r.dataset,r.model_unit)]<=best else 'MLH1-retained'
u['cls']=[unit_cls(r) for r in u.itertuples()]
chk=u.merge(U[['dataset','model_unit','histology']],on=['dataset','model_unit'],how='inner',suffixes=('','_u215'))
assert len(chk)==int((U.lineage=='uterus').sum()) and (chk.histology==chk.histology_u215).all(), 'uterine unit histology != inputs/units_215.csv'
u[['dataset','model_unit','histology','cls','organoid']].to_csv(OUT/'uterine_units_classes.csv',index=False)
print('MLH1 threshold',best,'; CPTAC concordant',int((E.rule_meth==E.truth).sum()),'/',len(E),'; uterine units',len(u),'(DNA-only NextGen',len(dna_only),')')
print(pd.crosstab(u.cls,u.dataset,margins=True).to_string())
