#!/usr/bin/env python3
"""Figures, xlsx tables, Table S2, number-audit (with manuscript line numbers) from results/<tag>.
--publish copies the package files into the manuscript folder (legacy duplicates are moved to figures/archive, never deleted)."""
import argparse, json, re, shutil, datetime, sys
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

PIPE=Path(__file__).resolve().parent
sys.path.insert(0,str(PIPE))
from paths import path as cfg_path, has as cfg_has   # input locations from config_paths.json (see paths.py)
MS=cfg_path('manuscript_dir') if cfg_has('manuscript_dir') else None
DELIV=cfg_path('deliverables_dir') if cfg_has('deliverables_dir') else None
ap=argparse.ArgumentParser(); ap.add_argument('--tag',required=True); ap.add_argument('--publish',action='store_true')
ap.add_argument('--manuscript',default=str(MS/'manuscript.md') if MS else '')
ap.add_argument('--write-expected',action='store_true',help='with --publish: refresh tests/expected (maintainers only)'); a=ap.parse_args()
STATIC=PIPE/'inputs/static'; XWALK=PIPE/'inputs/crosswalk'
# supplementary numbering (order of first citation in the manuscript; written by number_supplements.py)
SN=json.load(open(PIPE/'inputs/supp_numbering.json',encoding='utf-8'))
def STn(k): return SN['tables'][k]['n']
def SFn(k): return SN['figures'][k]['n']
def ST(k): return f"Supplementary Table S{STn(k)}"
def SF(k): return f"Supplementary Fig. S{SFn(k)}"
def TSTEM(k,suffix=''): return f"TableS{STn(k)}_{SN['tables'][k]['stem']}{suffix}"
def FSTEM(k): return f"FigS{SFn(k)}_{SN['figures'][k]['stem']}"
RES=PIPE/'results'/a.tag; FIG=RES/'figures'; TAB=RES/'tables'; FIG.mkdir(exist_ok=True); TAB.mkdir(exist_ok=True)
K=json.load(open(RES/'key_numbers.json'))
# Method B validity: runs of the fixed run_all.py carry methodB_qc.gene_order_fixed; older runs carry methodB_validity.json
if K.get('methodB_qc',{}).get('gene_order_fixed'): BV=True; BV_NOTE='Method B computed with the corrected gene ordering.'
elif (RES/'methodB_validity.json').exists():
    _bv=json.load(open(RES/'methodB_validity.json')); BV=bool(_bv['valid']); BV_NOTE=_bv['note']
else: BV=False; BV_NOTE='Method B validity unknown for this run; Method B values withheld.'
PENDING='pending (Method B rerun)'
print('Method B valid for',a.tag,':',BV)
U=pd.read_csv(RES/'units_215_calls.csv'); T2=pd.read_csv(RES/'table2_strict_loose.csv').set_index('stratum')
ES=pd.read_csv(RES/'equal_depth_summary.csv'); ER=pd.read_csv(RES/'equal_depth_reps.csv'); SE=pd.read_csv(RES/'reference_composition_sensitivity.csv')
MOL=cfg_path('molecular_dir')   # molecular derivations from public data (04/05/06 scripts)
MTH1=json.load(open(MOL/'mlh1_threshold.json'))['mlh1_log2tpm_threshold_rounded']   # MLH1 log2(TPM+1) threshold (Youden on CPTAC UCEC, in-sample; 05_derive_endometrial_classes.py)
DS_LABEL={'HCMI':'HCMI','DepMap_NextGen':'DepMap NextGen','DepMap_24Q4':'DepMap (2D lines)','CPTAC':'CPTAC','Sanger_organoids':'Sanger organoids',
          'GSE307568':'GSE307568','Berg2021_EMTAB10664':'E-MTAB-10664 (Berg)','GSE208216':'GSE208216'}
EXPR_LABEL='DepMap '+a.tag.split('_')[0] if a.tag.split('_')[0] in ('24Q4','26Q3') else a.tag
GID=K.get('gene_id','symbol'); GENE_MATCH='Entrez Gene ID' if GID=='entrez' else 'gene symbol'
# public dataset labels: the internal key 'DepMap_24Q4' denotes DepMap 2D lines, whatever the expression release
DS_PUBLIC={'DepMap_24Q4':'DepMap_2D','DepMap_NextGen':'DepMap_NextGen'}
CODE_NOTE=('Codes: histology OV_HGSC = ovarian high-grade serous; OV_serous = ovarian serous, grade not stated (HCMI/Sanger); OV_serous_NOS = DepMap ovarian serous, not otherwise specified; '
           'OV_LGSC = low-grade serous; OV_clear_cell, OV_endometrioid, OV_mucinous, OV_mixed, OV_mullerian_mixed (mixed Mullerian), OV_borderline, OV_unknown = other ovarian histologies; '
           'UT_endometrioid, UT_serous, UT_clear_cell, UT_carcinosarcoma, UT_mixed, UT_undifferentiated, UT_dedifferentiated, UT_adenosquamous, UT_EC_NOS (endometrial carcinoma NOS) = uterine histologies. '
           'stratum HGSC = OV_HGSC + OV_serous; serous_NOS = OV_serous_NOS; nonHGSC = other ovarian histologies; EC = endometrial carcinoma (non-carcinosarcoma); CS = uterine carcinosarcoma. '
           'TCGA projects: OV ovarian serous cystadenocarcinoma; UCEC uterine corpus endometrial carcinoma; UCS uterine carcinosarcoma; BRCA, COAD, STAD, CESC, KIRC, LUAD, PAAD, BLCA = breast, colon, stomach, cervical, clear-cell renal, lung adeno-, pancreatic and bladder carcinoma. '
           'dataset DepMap_2D = DepMap 2D cell lines (internal key DepMap_24Q4, whatever the expression release); DepMap_NextGen = NextGen organoids.')
def relabel(df,col='dataset'):
    df=df.copy()
    if col in df.columns: df[col]=df[col].replace(DS_PUBLIC)
    return df
plt.rcParams.update({'font.size':8,'axes.titlesize':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':6.5,
                     'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,stem):
    for ext in ('png','pdf'): fig.savefig(FIG/f'{stem}.{ext}',dpi=300,bbox_inches='tight')
    plt.close(fig)
def panel(ax,letter,x=-0.08): ax.text(x,1.06,letter,transform=ax.transAxes,fontsize=11,fontweight='bold',va='bottom',ha='right')
def frac(k,n): return f'{int(round(k))}/{int(n)}'

# ───────── Figure 1 ─────────
F=pd.read_csv(PIPE/'inputs/units_215.csv')
HG={'OV_HGSC','OV_serous','OV_serous_NOS'}
F['grp']=np.where(F.lineage=='uterus','Uterine',np.where(F.histology.isin(HG),'Ovarian HGSC/serous','Ovarian non-HGSC'))
order=['HCMI','DepMap_NextGen','DepMap_24Q4','Sanger_organoids','GSE307568','Berg2021_EMTAB10664','GSE208216']
ct=pd.crosstab(F.dataset,F.grp).reindex(order).fillna(0)[['Ovarian HGSC/serous','Ovarian non-HGSC','Uterine']]
IDDIR=RES/'identity' if (RES/'identity').exists() else PIPE/'results/26Q3/identity'   # sensitivity runs reuse the headline identity check
I=pd.read_csv(IDDIR/'nextgen_identity_check.csv'); IDS=json.load(open(IDDIR/'identity_summary.json'))
TH_E=IDS['thr_expr_99pct']; TH_C=IDS['thr_cn_max_diff_case']; TH_M=IDS['thr_mut_max_diff_case']
official={'ACH-003350','ACH-003723','ACH-003709','ACH-003710','ACH-003711'}
_i9=I.set_index('nextgen')
fig=plt.figure(figsize=(11.2,6.4)); gs=fig.add_gridspec(2,2,height_ratios=[1.05,1.0],hspace=0.45,wspace=0.3)
ax=fig.add_subplot(gs[0,:]); cols=['#1f77b4','#9ecae1','#d62728']; y=np.arange(len(ct)); left=np.zeros(len(ct))
for col,c in zip(ct.columns,cols): ax.barh(y,ct[col].values,left=left,color=c,label=col,height=0.72); left+=ct[col].values
for i,(idx,row) in enumerate(ct.iterrows()): ax.text(row.sum()+0.6,i,str(int(row.sum())),va='center',fontsize=7)
ax.set_yticks(y); ax.set_yticklabels([DS_LABEL[d] for d in ct.index]); ax.invert_yaxis(); ax.set_xlim(0,left.max()+12)
n_org=int(F.is_org.sum()); ax.set_xlabel('Independent expression units after patient-level deduplication')
ax.set_title(f'{len(F)} independent expression units ({n_org} organoids, {len(F)-n_org} 2D lines)',loc='left')
ax.legend(loc='lower right',frameon=False,ncol=3); panel(ax,'a')
c=['#d62728' if n in official else '#7f7f7f' for n in I.nextgen]
ax=fig.add_subplot(gs[1,0]); ax.scatter(I.best_cn_r,I.best_r,c=c,s=28,zorder=3,edgecolors='white',linewidths=0.3)
OFF_B={'ACH-003710':(-34,6),'ACH-003711':(6,-12),'ACH-003350':(4,4),'ACH-003723':(4,4),'ACH-003709':(4,4)}
for r in I.itertuples():
    if r.nextgen in official:
        ax.annotate(r.nextgen.replace('ACH-',''),(r.best_cn_r,r.best_r),textcoords='offset points',xytext=OFF_B.get(r.nextgen,(3,3)),fontsize=6,color='#d62728',
                    arrowprops=dict(arrowstyle='-',color='#d62728',lw=0.4,shrinkA=0,shrinkB=2))
ax.axvline(TH_C,ls='--',c='grey',lw=0.8); ax.axhline(TH_E,ls=':',c='grey',lw=0.8)
ax.set_xlabel(f'Best copy-number profile r vs HCMI models\n(dashed: maximum between different HCMI cases, {TH_C:.3f})')
ax.set_ylabel(f'Best expression residual r vs any cohort model\n(dotted: 99th percentile of all model pairs, {TH_E:.3f})')
ax.set_title('Expression and copy-number similarity (red: official NextGen–HCMI pairs)',loc='left')
ax.legend(handles=[Patch(color='#d62728',label='Official same-patient pair (5/5 recovered)'),Patch(color='#7f7f7f',label='Not called (incl. ACH-003903)')],loc='upper left',frameon=False); panel(ax,'b',-0.2)
ax=fig.add_subplot(gs[1,1]); ax.scatter(I.n_damaging_genes,I.shared_damaging_genes,c=c,s=28,zorder=3,edgecolors='white',linewidths=0.3); ax.set_xscale('log')
OFF_C={'ACH-003350':(-38,22),'ACH-003723':(-8,30),'ACH-003711':(14,22),'ACH-003710':(10,8),'ACH-003709':(-40,-4),'ACH-003903':(-12,8)}
for r in I.itertuples():
    if r.nextgen in official or r.nextgen=='ACH-003903':
        ax.annotate(r.nextgen.replace('ACH-',''),(r.n_damaging_genes,r.shared_damaging_genes),textcoords='offset points',xytext=OFF_C.get(r.nextgen,(3,3)),fontsize=6,
                    color='#d62728' if r.nextgen in official else '#333333',arrowprops=dict(arrowstyle='-',color='#999999',lw=0.4,shrinkA=0,shrinkB=2))
ax.set_xlabel('Damaging/truncating genes in NextGen model (excluding TP53)'); ax.set_ylabel('Shared truncating genes with best HCMI model')
r9=_i9.loc['ACH-003709']; r3=_i9.loc['ACH-003903']
ax.set_title(f"ACH-003709 shares {int(r9.shared_damaging_genes)} genes with {r9.best_hcmi_mut.split('-C5')[0]} (Jaccard {r9.jaccard:.3f});\nACH-003903 shares {int(r3.shared_damaging_genes)} with a hypermutated HCMI model (Jaccard {r3.jaccard:.3f}): not called",loc='left')
panel(ax,'c',-0.15)
save(fig,'Fig1_cohort_and_identity')

# ───────── Figure 2 ─────────
fig,axs=plt.subplots(1,2,figsize=(11,3.9),gridspec_kw={'width_ratios':[1.5,1]})
rows=[('HGSC/serous (incl. serous NOS)','HGSC/serous\n(all)'),('  OV_HGSC','HGSC'),('  DepMap serous NOS','DepMap\nserous NOS'),
      ('Endometrial carcinoma (non-CS)','Endometrial\ncarcinoma'),('Uterine carcinosarcoma','Carcino-\nsarcoma')]
ax=axs[0]; x=np.arange(len(rows)); w=0.38
for j,(key,lab) in enumerate(rows):
    r=T2.loc[key]; s=100*r.strict_A/r.n; l=100*r.loose_A/r.n
    ax.bar(j-w/2,s,w,color='#1f77b4',label='Strict' if j==0 else None); ax.bar(j+w/2,l,w,color='#aec7e8',label='Loose' if j==0 else None)
    ax.text(j-w/2,s+1.5,frac(r.strict_A,r.n),ha='center',fontsize=6,rotation=90); ax.text(j+w/2,l+1.5,frac(r.loose_A,r.n),ha='center',fontsize=6,rotation=90)
ax.set_xticks(x); ax.set_xticklabels([lab for _,lab in rows]); ax.set_ylim(0,118); ax.set_ylabel('Method A match rate (%)')
ax.set_title('Strict (site-exact) and loose (any gynaecological centroid: OV, UCEC or UCS) matching',loc='left'); ax.legend(frameon=False,loc='upper right',ncol=2); panel(ax,'a')
ax=axs[1]; cats=[('  EC organoids','EC\norganoids'),('  EC 2D lines','EC\n2D lines'),('  HGSC/serous organoids','HGSC/serous\norganoids'),('  HGSC/serous 2D lines','HGSC/serous\n2D lines')]
for j,(key,lab) in enumerate(cats):
    r=T2.loc[key]; s=100*r.strict_A/r.n; l=100*r.loose_A/r.n
    ax.bar(j-w/2,s,w,color='#2ca02c',label='Strict' if j==0 else None); ax.bar(j+w/2,l,w,color='#98df8a',label='Loose' if j==0 else None)
    ax.text(j-w/2,s+1.5,frac(r.strict_A,r.n),ha='center',fontsize=6,rotation=90); ax.text(j+w/2,l+1.5,frac(r.loose_A,r.n),ha='center',fontsize=6,rotation=90)
ax.axvline(1.5,ls='--',c='#dddddd',lw=0.8); ax.set_xticks(range(4)); ax.set_xticklabels([lab for _,lab in cats]); ax.set_ylim(0,118)
ax.set_ylabel('Method A match rate (%)'); ax.set_title('Organoids versus 2D lines',loc='left'); ax.legend(frameon=False,loc='upper right',ncol=2); panel(ax,'b')
plt.tight_layout(); save(fig,'Fig2_lineage_matching')

# ───────── Figure 3 (endometrial classes; source = Table S2 class columns) ─────────
S2prev=pd.read_csv(PIPE/'inputs/TableS2_prev.csv')
Ev=pd.read_csv(cfg_path('cptac_rule_evaluation'))
ordr=['POLEmut','MMRd','p53abn','NSMP']
fig,axs=plt.subplots(1,2,figsize=(10,3.6),gridspec_kw={'width_ratios':[1,1.6]})
m=pd.crosstab(Ev.truth,Ev.rule_meth).reindex(index=ordr,columns=ordr,fill_value=0)
axs[0].imshow(m.values,cmap='Greens',vmin=0)
for i in range(4):
    for j in range(4): axs[0].text(j,i,int(m.values[i,j]),ha='center',va='center',fontsize=8)
axs[0].set_xticks(range(4)); axs[0].set_xticklabels(ordr,rotation=30,ha='right'); axs[0].set_yticks(range(4)); axs[0].set_yticklabels(ordr)
axs[0].set_xlabel('DNA/methylation surrogate rule'); axs[0].set_ylabel('CPTAC UCEC genomic subtype')
axs[0].set_title(f'CPTAC UCEC: {int((Ev.rule_meth==Ev.truth).sum())}/{len(Ev)} concordant',loc='left'); axs[0].spines[:].set_visible(True); panel(axs[0],'a',-0.3)
uc=pd.read_csv(cfg_path('uterine_units_classes'))
uc['cls2']=uc.cls.replace({'MLH1-low':'MLH1-low (expression proxy)','MLH1-retained':'MLH1-retained (expression proxy)'})
dsu=['HCMI','DepMap_NextGen','DepMap_24Q4','GSE307568','Berg2021_EMTAB10664']
ctU=pd.crosstab(uc.dataset,uc.cls2).reindex(dsu).fillna(0)
cl=[c for c in ['POLEmut','MMRd','p53abn','NSMP','MLH1-low (expression proxy)','MLH1-retained (expression proxy)'] if c in ctU.columns]
pal={'POLEmut':'#2ca02c','MMRd':'#7570b3','p53abn':'#d95f02','NSMP':'#999999','MLH1-low (expression proxy)':'#fdd0a2','MLH1-retained (expression proxy)':'#f0f0f0'}
left=np.zeros(len(ctU))
for c_ in cl:
    axs[1].barh(range(len(ctU)),ctU[c_].values,left=left,color=pal[c_],edgecolor='#bbbbbb' if 'retained' in c_ else None,linewidth=0.4,label=c_); left+=ctU[c_].values
axs[1].set_yticks(range(len(ctU))); axs[1].set_yticklabels([DS_LABEL[d] for d in ctU.index]); axs[1].invert_yaxis()
axs[1].set_xlabel('Independent uterine units'); axs[1].legend(frameon=False,fontsize=6,loc='upper left',bbox_to_anchor=(1.01,1.0))
axs[1].set_title('Surrogate class by resource (DNA-classified units)',loc='left'); panel(axs[1],'b')
plt.tight_layout(); save(fig,'Fig3_endometrial_classes')

# ───────── Figure 4 (paired concordance; unchanged data, panel letters added) ─────────
P=pd.read_csv(cfg_path('paired_rna')).sort_values('rho_own',ascending=False).drop_duplicates('model')
Q=pd.read_csv(cfg_path('paired_methylation')).sort_values('r_own',ascending=False).drop_duplicates('model')
Mt=pd.read_csv(cfg_path('paired_mutation')).sort_values('tumour_muts_in_model',ascending=False).drop_duplicates('model')
fig,axs=plt.subplots(1,3,figsize=(9.6,3.4))
for ax,(d,own,oth,t),L in zip(axs[:2],[(P,'rho_own','rho_other_median','RNA (Spearman ρ, 2,000 genes)'),(Q,'r_own','r_other_median','Methylation (Pearson r, top 10k CpGs)')],'ab'):
    cc=['#3b6ea8' if h.startswith('OV') else '#c0504d' for h in d.histology]
    for i,(x0,x1) in enumerate(zip(d[oth],d[own])): ax.plot([0,1],[x0,x1],color=cc[i],alpha=0.6,lw=0.8,marker='o',ms=3)
    ax.set_xticks([0,1]); ax.set_xticklabels(['median vs\nother tumours','vs own\ntumour']); ax.set_xlim(-0.3,1.3)
    ax.set_title(t+f'\nn={len(d)} models / {d.case.nunique()} cases'); panel(ax,L)
ax=axs[2]; cc=['#3b6ea8' if h.startswith('OV') else '#c0504d' for h in Mt.histology]
ax.scatter(np.log10(Mt.n_tumour),Mt.tumour_muts_in_model,c=cc,s=14); ax.set_ylim(0,1.05)
ax.set_xlabel('log10 tumour coding mutations'); ax.set_ylabel('Fraction of tumour mutations\nfound in model')
ax.set_title(f'Mutations (masked MAF)\nn={len(Mt)} models / {Mt.case.nunique()} cases'); panel(ax,'c')
axs[0].legend(handles=[Patch(color='#3b6ea8',label='Ovary'),Patch(color='#c0504d',label='Uterus')],frameon=False)
plt.tight_layout(); save(fig,'Fig4_paired_concordance')
paired=dict(rna_own_median=round(float(P.rho_own.median()),3),rna_other_median=round(float(P.rho_other_median.median()),3))

# ───────── Supplementary figure: ovarian histotype, three centring regimes (number from supp_numbering.json) ─────────
O=U[U.lineage=='ovary'].copy(); O['label']=O.label.replace({'other/unknown':'other'})
rws=['HGSC','LGSC','CCC','EC','MC','other']; cls_=['HGSC','LGSC','CCC','EC','MC']
fig,axs=plt.subplots(1,3,figsize=(11,3.5))
for ax,(t,col),L in zip(axs,[('Uncentred','hist_uncentred'),('Ordinary centring','hist_ordinary'),('Label-balanced centring (optimistic)','hist_label_balanced')],'abc'):
    mm=pd.crosstab(O.label,O[col]).reindex(index=rws,columns=cls_,fill_value=0)
    ax.imshow(mm.values,cmap='Blues'); ax.spines[:].set_visible(True)
    for i in range(len(rws)):
        for j in range(len(cls_)): ax.text(j,i,int(mm.values[i,j]),ha='center',va='center',fontsize=7)
    ax.set_xticks(range(len(cls_))); ax.set_xticklabels(cls_); ax.set_yticks(range(len(rws))); ax.set_yticklabels(rws)
    nh=O[O.label.isin(['CCC','EC','MC','LGSC'])]; hh=O[O.label=='HGSC']
    ax.set_xlabel('Nearest GSE295399 histotype centroid'); ax.set_ylabel('Model label')
    ax.set_title(f"{t}\nnon-HGSC {int((nh.label==nh[col]).sum())}/{len(nh)}; HGSC {int((hh[col]=='HGSC').sum())}/{len(hh)}",loc='left'); panel(ax,L)
plt.tight_layout(); save(fig,FSTEM('histotype'))

# ───────── Supplementary Figs S1 (Method A confusion heatmap) and S2 (Method B heatmap; numbers from supp_numbering.json) ─────────
proj_order=['OV','UCEC','UCS','BRCA','CESC','STAD','COAD','PAAD','LUAD','KIRC','BLCA']
hist_order=['OV_HGSC','OV_serous','OV_serous_NOS','OV_LGSC','OV_clear_cell','OV_endometrioid','OV_mucinous','OV_mixed','OV_mullerian_mixed','OV_borderline','OV_unknown',
            'UT_endometrioid','UT_EC_NOS','UT_serous','UT_clear_cell','UT_adenosquamous','UT_mixed','UT_dedifferentiated','UT_undifferentiated','UT_carcinosarcoma']
def heat(col,title,cmap,stem):
    fig,ax=plt.subplots(1,1,figsize=(6.4,6.2))
    mm=pd.crosstab(U.histology,U[col]).reindex(index=hist_order,columns=proj_order,fill_value=0)
    ax.imshow(mm.values,cmap=cmap,aspect='auto'); ax.spines[:].set_visible(True)
    for i in range(mm.shape[0]):
        for j in range(mm.shape[1]):
            v=int(mm.values[i,j])
            if v: ax.text(j,i,v,ha='center',va='center',fontsize=6,color='white' if v>12 else 'black')
    ax.set_xticks(range(len(proj_order))); ax.set_xticklabels(proj_order,rotation=45,ha='right'); ax.set_yticks(range(len(hist_order))); ax.set_yticklabels([h.replace('_',' ').replace('mullerian','Müllerian') for h in hist_order])
    ax.set_xlabel('Nearest TCGA project centroid (11 projects, full cohorts)'); ax.set_title(title,loc='left')
    plt.tight_layout(); save(fig,stem)
heat('nearest_A','Method A (Spearman; primary)','Purples',FSTEM('confusion'))
# Method B (Celligner-like): shown for comparison with Method A; no conclusions drawn
if BV: heat('nearest_B','Method B (Celligner-like) assignments; failed stability criterion C2','Greys',FSTEM('methodB'))

# ───────── Supplementary figure: equal-depth subsampling (number from supp_numbering.json) ─────────
st=[('HGSC_serous','HGSC/serous (incl. NOS)','HGSC/serous (incl. serous NOS)'),('EC_nonCS','EC (non-CS)','Endometrial carcinoma (non-CS)'),('EC_organoid','EC organoids','  EC organoids'),
    ('EC_2D','EC 2D lines','  EC 2D lines'),('HGSC_organoid','HGSC/serous organoids','  HGSC/serous organoids'),('HGSC_2D','HGSC/serous 2D','  HGSC/serous 2D lines'),('CS','Carcinosarcoma','Uterine carcinosarcoma')]
fig,axs=plt.subplots(1,2,figsize=(11,3.6),sharey=True)
for ax,rule,L in zip(axs,['strict','loose'],'ab'):
    data=[100*ER[f'{k}_{rule}']/ER[f'{k}_n'] for k,_,_ in st]
    ax.boxplot(data,widths=0.5,showfliers=False,medianprops=dict(color='#333333'))
    for i,d in enumerate(data): ax.scatter(np.random.default_rng(i).normal(i+1,0.06,len(d)),d,s=6,color='#1f77b4',alpha=0.6,zorder=3)
    full=[100*T2.loc[k3][f'{rule}_A']/T2.loc[k3].n for _,_,k3 in st]
    ax.scatter(range(1,len(st)+1),full,marker='D',color='#d62728',s=22,zorder=4,label='Full 11-project reference')
    ax.set_xticks(range(1,len(st)+1)); ax.set_xticklabels([l for _,l,_ in st],rotation=35,ha='right'); ax.set_ylim(0,105)
    ax.set_ylabel('Method A match rate (%)'); ax.set_title(f'{rule.capitalize()} rule: 40 tumours/project × {len(ER)} repeats',loc='left'); panel(ax,L)
axs[1].tick_params(axis='y',labelleft=True)
axs[0].legend(frameon=False,loc='lower left'); plt.tight_layout(); save(fig,FSTEM('equal_depth'))

# ───────── Supplementary figure: reference composition on the same 215 units (number from supp_numbering.json) ─────────
lab_ref={'Full 11 projects (headline)':'Full cohorts, 11 projects','Gyn full + non-gyn 40/project':'Gynaecological full + non-gynaecological 40/project',
         'Gyn ~200/project + non-gyn full':'Gynaecological ~200/project + non-gynaecological full'}
SE['label']=[lab_ref.get(r,r.replace(' (mean of ',' (mean, ').replace(')',' repeats)') if 'Equal' in r else r) for r in SE.reference]
sts=[('HGSC_serous','HGSC/serous'),('EC_nonCS','EC (non-CS)'),('EC_organoid','EC organoids'),('EC_2D','EC 2D lines'),('CS','Carcinosarcoma')]
fig,axs=plt.subplots(1,2,figsize=(11.5,3.8),sharey=True); colors=['#1f77b4','#ff7f0e','#2ca02c','#9467bd']
for ax,rule,L in zip(axs,['strict','loose'],'ab'):
    w=0.8/len(SE)
    for k,(i,r) in enumerate(SE.iterrows()):
        vals=[100*r[f'{s}_{rule}']/r[f'{s}_n'] for s,_ in sts]
        ax.bar(np.arange(len(sts))+(k-(len(SE)-1)/2)*w,vals,w,color=colors[k%4],label=r.label)
    ax.set_xticks(range(len(sts))); ax.set_xticklabels([l for _,l in sts]); ax.set_ylim(0,105)
    ax.set_ylabel('Method A match rate (%)'); ax.set_title(f'{rule.capitalize()} rule; same 215 units and feature genes',loc='left'); panel(ax,L)
axs[1].tick_params(axis='y',labelleft=True)
h_,l_=axs[0].get_legend_handles_labels()
fig.legend(h_,l_,frameon=False,fontsize=6.5,loc='lower center',ncol=len(l_),bbox_to_anchor=(0.5,-0.01))
plt.tight_layout(rect=(0,0.07,1,1)); save(fig,FSTEM('refcomp'))

# ───────── Supplementary Fig: published Celligner context (redrawn from the same inputs; panel letters added) ─────────
Rk=pd.read_csv(STATIC/'hcmi_S12_gyn_own_tumour_rank.csv')
fig,axs=plt.subplots(1,2,figsize=(9,3.3))
axs[0].hist(Rk.own_rank,bins=[1,2,5,10,25,50,100],color='#6baed6',edgecolor='k'); axs[0].set_xscale('log')
axs[0].set_xlabel('Rank of own tumour among 311 HCMI tumours\n(HCMI Supplementary Table S12, Celligner distance)'); axs[0].set_ylabel('Gynaecological model files')
axs[0].set_title(f'Own tumour ranked first in {int((Rk.own_rank==1).sum())}/{len(Rk)} model files',loc='left'); panel(axs[0],'a')
# published nearest-neighbour counts (NextGen paper Celligner nearest-neighbour table), as in the earlier export
vals={'NextGen ovary':(3,13),'NextGen uterus':(2,4),'2D ovary':(21,66),'2D uterus':(19,43)}
axs[1].bar(range(len(vals)),[k/n for k,n in vals.values()],color=['#9e9ac8','#9e9ac8','#bdbdbd','#bdbdbd'])
for i,(k,n) in enumerate(vals.values()): axs[1].text(i,k/n+0.02,f'{k}/{n}',ha='center',fontsize=7)
axs[1].set_xticks(range(len(vals))); axs[1].set_xticklabels(list(vals)); axs[1].set_ylim(0,1)
axs[1].set_ylabel('Rank-1 Celligner neighbour\nin the same lineage (fraction)'); axs[1].set_title('NextGen paper nearest-neighbour table (published Celligner)',loc='left'); panel(axs[1],'b')
plt.tight_layout(); save(fig,FSTEM('celligner'))
# HCMI endometrial surrogate-feature heatmap (06_fig_hcmi_surrogate.py)
for ext in ('png','pdf'): shutil.copy2(MOL/f'hcmi_endometrial_surrogate.{ext}',FIG/f"{FSTEM('surrogate')}.{ext}")

# ───────── Table S2 (model selection annotation) ─────────
M26=pd.read_csv(cfg_path('depmap_model_csv')).set_index('ModelID')
prev=S2prev.set_index(['dataset','model_unit'])
S2=U.copy()
UC=pd.read_csv(cfg_path('uterine_units_classes'))   # DNA surrogate classes / MLH1 expression proxy per uterine unit (05_derive_endometrial_classes.py)
S2=S2.join(UC.set_index(['dataset','model_unit']).cls.rename('DNA_surrogate_class'),on=['dataset','model_unit']).join(prev[['HCMIID','ACH','SIDM']],on=['dataset','model_unit'])
XO=UC[~UC.set_index(['dataset','model_unit']).index.isin(U.set_index(['dataset','model_unit']).index)]   # DNA-only units outside the 215 expression units (NextGen ACH-003239)
_h=S2[S2.lineage=='uterus'].join(UC.set_index(['dataset','model_unit']).histology.rename('histology_classes'),on=['dataset','model_unit'])
if _h.histology_classes.isna().any() or (_h.histology_classes!=_h.histology).any():   # every uterine unit: Table S2/S10 histology == derivation (fail loudly)
    raise SystemExit('uterine unit histology differs between Table S2 and uterine_units_classes.csv:\n'+_h[_h.histology_classes!=_h.histology][['dataset','model_unit','histology','histology_classes']].to_string())
S2['MLH1_expression_proxy']=np.where(S2.DNA_surrogate_class.isin(['MLH1-low','MLH1-retained']),S2.DNA_surrogate_class,
    np.where(S2.dataset.isin(['GSE307568','Berg2021_EMTAB10664']),np.where(S2.MLH1_log2tpm<=MTH1,'MLH1-low','MLH1-retained'),''))
S2['DNA_surrogate_class']=S2.DNA_surrogate_class.where(~S2.DNA_surrogate_class.isin(['MLH1-low','MLH1-retained']),'')
S2['DNA_surrogate_class']=S2.DNA_surrogate_class.fillna('')
S2.loc[S2.lineage=='ovary','MLH1_expression_proxy']=''
def ach(r):
    if isinstance(r.ACH,str) and r.ACH.startswith('ACH-'): return r.ACH
    if str(r.model_unit).startswith('ACH-'): return r.model_unit
    if isinstance(r.HCMIID,str):
        hit=M26.index[M26.HCMIID.astype(str).str.startswith(r.model_unit)]
        if len(hit): return ';'.join(hit)
    return ''
S2['ACH']=[ach(r) for r in S2.itertuples()]
def cat(achs,field):
    out=[]
    for x in [t for t in str(achs).split(';') if t.startswith('ACH-')]:
        if x in M26.index and pd.notna(M26.loc[x,field]): out.append(str(M26.loc[x,field]).strip())
    return ';'.join(out)
S2['source_repository_26Q3']=[cat(x,'SourceType') for x in S2.ACH]
S2['catalog_number_26Q3']=[cat(x,'CatalogNumber') for x in S2.ACH]
S2['ATCC_catalog']=[';'.join(c.upper() if c.lower().startswith('crl') else c for c,s in zip(cn.split(';'),sr.split(';')) if s=='ATCC') if cn else ''
                   for cn,sr in zip(S2.catalog_number_26Q3,S2.source_repository_26Q3)]
S2['RRID']=[cat(x,'RRID') for x in S2.ACH]
S2['catalog_number_status']=np.where(S2.catalog_number_26Q3!='','listed in DepMap 26Q3 Model.csv CatalogNumber',
    np.where(S2.ACH!='','not listed in DepMap 26Q3 Model.csv CatalogNumber','no DepMap model ID (not applicable)'))
S2['expression_source']=np.where(S2.dataset=='DepMap_24Q4',EXPR_LABEL,np.where(S2.dataset=='DepMap_NextGen','NextGen manuscript file set' if a.tag.startswith('24Q4') else EXPR_LABEL+' (default profile)','resource-native'))
S2['resource']=S2.dataset.map(DS_LABEL)
keep=['resource','model_unit','patient','lineage','histology','stratum','model_format','is_org','nearest_A','strict_A','loose_A',
      'hist_uncentred','hist_ordinary','hist_label_balanced','DNA_surrogate_class','MLH1_expression_proxy','MLH1_log2tpm','WT1_log2tpm','WT1_status','KRT20_log2tpm','CDX2_log2tpm',
      'mucinous_GI_like','HCMIID','ACH','SIDM','RRID','source_repository_26Q3','catalog_number_26Q3','ATCC_catalog','catalog_number_status','expression_source','indep_key']
S2=S2[keep].rename(columns={'nearest_A':'nearest_TCGA_MethodA','nearest_B':'nearest_TCGA_MethodB','strict_A':'strict_match_A','loose_A':'loose_match_A','strict_B':'strict_match_B','loose_B':'loose_match_B',
    'hist_uncentred':'histotype_uncentred','hist_ordinary':'histotype_ordinary_centring','hist_label_balanced':'histotype_label_balanced','is_org':'organoid'})
for c in ['MLH1_log2tpm','WT1_log2tpm','KRT20_log2tpm','CDX2_log2tpm']: S2[c]=S2[c].round(3)
S2.to_csv(TAB/(TSTEM('selection')+'.csv'),index=False)
s2_stats=dict(ATCC_filled=int((S2.ATCC_catalog!='').sum()),any_catalog=int((S2.catalog_number_26Q3!='').sum()),
              units_with_ACH=int((S2.ACH!='').sum()),WT1_filled=int(S2.WT1_log2tpm.notna().sum()),
              DNA_class=S2.DNA_surrogate_class.replace('', 'none').value_counts().to_dict(),MLH1_proxy=S2.MLH1_expression_proxy.replace('', 'none').value_counts().to_dict())

# ───────── xlsx tables ─────────
def three_line(path):
    """Three-line table style (round 6, formatting only; cell values unchanged): thick top rule above the header,
    thin rule below the header, thick bottom rule below the last row; no vertical or other horizontal rules.
    The 'Notes' sheet holds the table footnotes and is left without rules."""
    from openpyxl import load_workbook
    from openpyxl.styles import Border, Side
    wb=load_workbook(path); thick=Side(style='medium'); thin=Side(style='thin')
    for ws in wb.worksheets:
        for row in ws.iter_rows():
            for c in row: c.border=Border()
        if ws.title=='Notes' or ws.max_row<1: continue
        last=ws.max_row
        for col in range(1,ws.max_column+1):
            ws.cell(1,col).border=Border(top=thick,bottom=thin if last>1 else thick)
            if last>1: ws.cell(last,col).border=Border(bottom=thick)
    wb.save(path)
def xl(name,sheets,notes=None):
    with pd.ExcelWriter(TAB/name) as w:
        for sh,df in sheets.items(): df.to_excel(w,sheet_name=sh[:31],index=False)
        if notes: pd.DataFrame({'Notes':notes}).to_excel(w,sheet_name='Notes',index=False)
    three_line(TAB/name)
t1=F.assign(Resource=F.dataset.map(DS_LABEL),Lineage=F.lineage.str.capitalize()).groupby(['Resource','Lineage'],sort=False).agg(Independent_units=('model_unit','count'),Organoid=('is_org','sum')).reset_index()
t1['Model_type']=np.where(t1.Organoid==t1.Independent_units,'Organoid','2D'); t1=t1.drop(columns='Organoid')
t1=pd.concat([t1,pd.DataFrame([dict(Resource='Total',Lineage='',Independent_units=len(F),Model_type=f'{int(F.is_org.sum())} organoids / {int((~F.is_org).sum())} 2D')])])
xl('Table1_cohort_composition.xlsx',{'Table1':t1},['Independent expression units after patient-level deduplication (one unit per patient).',
   'Five HCMI-derived NextGen organoids (HCM-WCMC-0673-C56-A/-B, HCM-CSHL-0806/-0810/-0811) are counted under HCMI; ACH-002419 is counted under HCMI (HCM-BROD-1117-C56).',
   'DNA-only NextGen ACH-003239 is outside the 215 expression units.'])
from table2 import build_table2
t2=build_table2(T2,K); cp=K['cptac']
rc=K['reference_counts']; esc=K['expr_source_counts']
n2d_prof=sum(v for k,v in esc.items() if k.startswith('DepMap_24Q4|')); nng_prof=sum(v for k,v in esc.items() if k.startswith('DepMap_NextGen|'))
n2d_unit=int((F.dataset=='DepMap_24Q4').sum()); nng_unit=int((F.dataset=='DepMap_NextGen').sum())
xl('Table2_lineage_fidelity.xlsx',{'Table2':t2},[
   f'Reference: full open TCGA primary tumours from 11 projects (GDC {K["gdc_release"]}), one primary tumour per case, {sum(rc.values()):,} tumours: '+', '.join(f'{k} {v:,}' for k,v in rc.items())+'.',
   f'Model expression: {EXPR_LABEL} for DepMap 2D lines'+('; NextGen organoids from the NextGen manuscript file set' if a.tag.startswith('24Q4') else ' and NextGen organoids (default profile per model)')+f'; other resources: resource-native processed data. Gene key: {GENE_MATCH}.',
   f'DepMap expression profiles in the harmonised matrix: {n2d_prof} 2D-line profiles and {nng_prof} NextGen organoid profiles. They map to {n2d_unit} independent DepMap 2D units and {nng_unit} independent NextGen units; the remaining {n2d_prof-n2d_unit} 2D profiles are sibling lines of a patient already represented (6) or a line counted under its HCMI unit (ACH-002419), and the remaining {nng_prof-nng_unit} NextGen profiles are the five HCMI-derived organoids counted under HCMI. Profiles that are not units are not scored.',
   'Strict: ovary-labelled -> OV; endometrial carcinoma (non-carcinosarcoma) -> UCEC; carcinosarcoma -> UCS. Loose: nearest centroid in {OV, UCEC, UCS} (any gynaecological TCGA project). No sample exclusions.',
   f'Ovarian non-HGSC units (n=35) are a TCGA-OV coverage gap and are reported in {ST("confusion")}, not here.'])
dna=S2[S2.DNA_surrogate_class.isin(['POLEmut','MMRd','p53abn','NSMP'])]
t3=pd.DataFrame({'Class':['MMRd','p53abn','NSMP','POLEmut']})
xo4=XO[XO.cls.isin(t3.Class)]
t3[f'All DNA-classified units (n={len(dna)+len(xo4)})']=[int((dna.DNA_surrogate_class==c).sum())+int((xo4.cls==c).sum()) for c in t3.Class]
t3[f'Organoid DNA-classified (n={int(dna.organoid.sum())+int(xo4.organoid.sum())})']=[int(((dna.DNA_surrogate_class==c)&dna.organoid).sum())+int(((xo4.cls==c)&xo4.organoid).sum()) for c in t3.Class]
xl('Table3_endometrial_classes.xlsx',{'Table3':t3},['DNA/methylation-derived surrogates (not IHC). Counts include the DNA-only NextGen organoid '+', '.join(f'{r.model_unit} ({r.cls})' for r in xo4.itertuples())+', outside the 215 expression units.',
   f'CPTAC UCEC rule accuracy {int((Ev.rule_meth==Ev.truth).sum())}/{len(Ev)}.'])
s1=pd.DataFrame([
 ['HCMI (HCMI-CMDC)','NCI GDC, open access',f'RNA STAR counts, SeSAMe EPIC betas, masked MAF, ASCAT gene-level CN (GDC {K["gdc_release"]})'],
 ['DepMap Public 26Q3 Model.csv','DepMap portal (downloaded 2026-10-07)','Model metadata; HCMIID / SangerModelID crosswalk; catalogue numbers'],
 [f'DepMap model expression ({EXPR_LABEL})','figshare 27993248 (DepMap Public 24Q4)' if a.tag.startswith('24Q4') else f'DepMap portal, DepMap Public 26Q3 OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv (downloaded 2026-10-08; IsDefaultEntryForModel = Yes)',f'OmicsExpression TPM log2(x+1), protein-coding genes; matched to the harmonised gene space by {GENE_MATCH}; md5 '+K['depmap_expr_md5']],
 ['HGNC complete set','HGNC (storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt; file dated 2026-10-06, downloaded 2026-10-08, md5 af43fd56a64ad57c0fcc72d6ef65e516)','Official Ensembl -> Entrez and symbol -> Entrez gene-identifier mapping (rules in Methods 2.2)'],
 ['DepMap Public 24Q4 expression (sensitivity analyses; confirmatory identity re-run)','figshare 27993248','OmicsExpressionProteinCodingGenesTPMLogp1.csv'],
 ['DepMap NextGen Model Manuscript 2026','DepMap portal file set; figshare 29472362','19 gynaecological NextGen organoids (expression, mutations, CN); mutations and CN used for the identity check; expression used in the 24Q4 sensitivity analysis and the confirmatory identity re-run'],
 ['NextGen Nature Supplementary Table 1','Nature ESM (downloaded 2026-10-07)','Official HCMI Cell Line Names for NextGen models'],
 ['Sanger tumour-derived organoids','figshare 28339340','Ovarian organoid RNA and annotations'],
 ['GSE307568','GEO (Chung et al., Cancer Res Commun 2026)','Endometrial PDO bulk RNA'],['E-MTAB-10664','ArrayExpress (Berg et al., Commun Med 2021)','Endometrial PDO RNA'],['GSE208216','GEO (Vias et al., eLife 2023)','HGSC PDO RNA (11 PDOs; fallopian-tube organoids excluded)'],
 ['GSE295399','GEO (Ittner et al., Sci Rep 2025)','Ovarian multi-histotype tumour RNA (centroids)'],['GSE226870','GEO (Beddows et al., Cancer Res 2024)','External CCC/EC tumour RNA check'],['GSE121103','GEO (Corona et al., Nat Commun 2020)','External multi-histotype tumour RNA check'],
 ['TCGA OV / UCEC / UCS',f'NCI GDC, open STAR - Counts ({K["gdc_release"]})','Full cohorts, one primary tumour per case: OV 426, UCEC 545, UCS 57'],
 ['TCGA BRCA, COAD, STAD, CESC, KIRC, LUAD, PAAD, BLCA',f'NCI GDC, open STAR - Counts ({K["gdc_release"]})','Full cohorts, one primary tumour per case: '+', '.join(f'{k} {v}' for k,v in K['reference_counts'].items() if k not in ('OV','UCEC','UCS'))],
 ['CPTAC UCEC (CPTAC-3)','NCI GDC open STAR counts; cBioPortal ucec_cptac_2020 (Dou et al., Cell 2020)','Transfer check; genomic subtype for surrogate-rule evaluation'],
 ['CPTAC OV (CPTAC-2, dbGaP phs000892)','NCI GDC open STAR counts (McDermott et al., Cell Rep Med 2020)','Transfer check (71 tumours)'],
 ['Cell Model Passports model_list','cellmodelpassports.sanger.ac.uk','SIDM crosswalk only (RNA matrices returned HTTP 403)']],columns=['Resource','Access / version','Content used'])
xl(TSTEM('sources','.xlsx'),{'Sources':s1})
S2x=S2.copy()
xl(TSTEM('selection','.xlsx'),{'Selection':S2x},['One row per independent expression unit (n=215). DNA_surrogate_class = DNA/methylation surrogate (not IHC); MLH1_expression_proxy = expression-only MLH1 status for resources without DNA data (not a molecular class).',
   f'MLH1-low threshold log2(TPM+1) <= {MTH1} (Youden index on CPTAC UCEC; in-sample). WT1 status threshold log2(TPM+1) >= 3.',
   'source_repository_26Q3 / catalog_number_26Q3 / RRID from DepMap 26Q3 Model.csv (SourceType, CatalogNumber, RRID); ATCC_catalog lists ATCC numbers only. HCMI units take the catalogue number of the DepMap model carrying the same HCMIID where one exists.',
   'Empty catalog_number_26Q3 / ATCC_catalog cells: CatalogNumber not listed in DepMap 26Q3 Model.csv (column catalog_number_status); no catalogue number was inferred from other sources. HCMI models are identified by their public HCMI IDs (model_unit, HCMIID).',
   CODE_NOTE,
   f'Expression-derived columns use {EXPR_LABEL} for DepMap units (expression_source column).',f'Method B (Celligner-like) per-unit calls are in {ST("methodB")} only (supplementary; no conclusions drawn).'])
I3=I.drop(columns=['best_sample']).copy()
I3['best_dataset']=I3.best_dataset.replace({'DepMap_24Q4':'DepMap 2D','DepMap_NextGen':'DepMap NextGen'})
xl(TSTEM('identity','.xlsx'),{'Identity':I3},[f"Expression data: {IDS['expression_source']}; other resources resource-native. Gene key: {'Entrez Gene ID' if IDS['gene_id']=='entrez' else 'gene symbol'} (expression, copy number and mutations).",
   f"Expression residual correlation over the headline feature genes ({IDS['feature_genes_used']} of {IDS['feature_genes']} present without missing values in all model profiles; gene mean across all models subtracted).",
   f"Criterion (i): expression residual r > {TH_E:.3f} (99th percentile of all non-NextGen model pairs). Criterion (ii): CN r > {TH_C:.3f} (maximum between different HCMI cases; {IDS['cn_genes']:,} genes).",
   f"Criterion (iii): shared truncating genes > {TH_M}, the observed maximum between different HCMI cases in the same cohort (an operational rule set from negative controls before the NextGen pairs were examined, not a pre-registered constant), interpreted together with mutation burden (Jaccard) and specificity (max_shared_same_hcmi_other_nextgen) - see Methods 2.3.",
   'Re-run 2026-10-08 with Entrez gene keys and DepMap 26Q3 expression (unchanged criteria). The earlier run (DepMap 24Q4 expression, NextGen manuscript file set, symbol matching, earlier feature genes; thresholds 0.599 and 0.597) is superseded; it made identical calls, as did re-runs with the original expression inputs and Entrez keys (archived: source_data/archive/24Q4/identity_*) and with DepMap 26Q3 expression and symbol keys (source_data/26Q3_symbol/identity_*).'])
uu=S2[S2.lineage=='uterus'][['resource','model_unit','histology','organoid','DNA_surrogate_class','MLH1_expression_proxy','MLH1_log2tpm']]
uu=pd.concat([uu,pd.DataFrame([dict(resource=DS_LABEL[r.dataset],model_unit=r.model_unit,histology=r.histology,organoid=bool(r.organoid),DNA_surrogate_class=r.cls,MLH1_expression_proxy='',MLH1_log2tpm=np.nan) for r in XO.itertuples()])])
xl(TSTEM('uterine','.xlsx'),{'Uterine':uu},[f'Independent uterine units ({int((S2.lineage=="uterus").sum())} expression units plus the DNA-only NextGen organoid '+', '.join(XO.model_unit)+').'])
O2=U[U.lineage=='ovary'][['dataset','model_unit','histology','label','hist_uncentred','hist_ordinary','hist_label_balanced']].assign(resource=lambda d:d.dataset.map(DS_LABEL)).drop(columns='dataset')
xl(TSTEM('histotype','.xlsx'),{'calls':O2,'PPV_uncentred':pd.read_csv(RES/'histotype_call_PPV.csv'),
   'confusion_uncentred':pd.read_csv(RES/'histotype_confusion_uncentred.csv'),'confusion_ordinary':pd.read_csv(RES/'histotype_confusion_ordinary.csv'),
   'confusion_label_balanced':pd.read_csv(RES/'histotype_confusion_label_balanced.csv')},
   ['Histotype codes in calls/PPV/confusion sheets: HGSC, LGSC, CCC (clear cell), EC (endometrioid), MC (mucinous), other/unknown; label = model histology mapped to these codes.',CODE_NOTE])
xl(TSTEM('cptac','.xlsx'),{'CPTAC_UCEC':Ev})
cmA=pd.read_csv(RES/'confusion_matrix_methodA.csv')
T2x=T2.reset_index(); EC2=relabel(pd.read_csv(RES/'EC_2D_units_calls.csv'))
T2x=T2x.drop(columns=[c for c in T2x.columns if c.endswith('_B')]); EC2=EC2.drop(columns=[c for c in EC2.columns if c.endswith('_B')])
sheets={'MethodA_histology':cmA,'MethodA_stratum':pd.read_csv(RES/'confusion_by_stratum_methodA.csv'),'Table2_full':T2x,'EC_2D_units':EC2}
xl(TSTEM('confusion','.xlsx'),sheets,[f'Model histology x nearest of 11 TCGA project centroids; expression {EXPR_LABEL} for DepMap units; dataset DepMap_2D = DepMap 2D lines.',f'Method B (Celligner-like) confusion matrices are in {ST("methodB")} (supplementary only).',CODE_NOTE])
RSYM=PIPE/'results'/(a.tag+'_symbol')
GM=None
if GID=='entrez' and (RSYM/'table2_strict_loose.csv').exists():
    Ts=pd.read_csv(RSYM/'table2_strict_loose.csv').set_index('stratum'); Ks=json.load(open(RSYM/'key_numbers.json'))
    gm=[]
    for key in ['HGSC/serous (incl. serous NOS)','  OV_HGSC','  DepMap serous NOS','  HGSC/serous organoids','  HGSC/serous 2D lines','Endometrial carcinoma (non-CS)','  EC organoids','  EC 2D lines','Uterine carcinosarcoma','Ovarian non-HGSC (coverage gap)']:
        e=T2.loc[key]; y=Ts.loc[key]
        gm.append({'Stratum':key.strip().replace('OV_HGSC','HGSC'),'Units':int(e.n),
                   'Strict A (Entrez, headline)':f'{int(e.strict_A)}/{int(e.n)}','Strict A (symbol)':f'{int(y.strict_A)}/{int(y.n)}',
                   'Loose A (Entrez, headline)':f'{int(e.loose_A)}/{int(e.n)}','Loose A (symbol)':f'{int(y.loose_A)}/{int(y.n)}'})
    GM=pd.DataFrame(gm)
    us=pd.read_csv(RSYM/'units_215_calls.csv'); ch=dict(nearest_A_changed=int((us.nearest_A!=U.nearest_A).sum()),nearest_B_changed=int((us.nearest_B!=U.nearest_B).sum()),
       strict_A_changed=int((us.strict_A!=U.strict_A).sum()),strict_B_changed=int((us.strict_B!=U.strict_B).sum()),
       absent_feature_genes_symbol=Ks.get('feature_genes_absent_from_depmap'),absent_feature_genes_entrez=K.get('feature_genes_absent_from_depmap'),
       units_nearest_A_changed=us[us.nearest_A!=U.nearest_A].model_unit.tolist())
    json.dump(ch,open(RES/'gene_matching_sensitivity.json','w'),indent=1)
_sens={'equal_depth_summary':ES,'equal_depth_repeats':ER,'reference_composition':SE.drop(columns=['label'])}
if GM is not None: _sens['gene_matching_symbol']=GM
xl(TSTEM('sensitivity','.xlsx'),_sens,
   ['Equal depth: 40 primary tumours per project drawn at random from the full cohorts of all 11 projects, 25 repeats (seed 20261008); feature genes fixed to the headline set.',
    'Reference composition: same 215 units, same rules and feature genes; gynaecological ~200/project = OV 200 / UCEC 200 / UCS 57 cases of the original subsample.']+
   (['Gene matching (sheet gene_matching_symbol): the full analysis repeated with gene symbols as the gene key (TCGA: GENCODE v36 gene names; DepMap: symbol of the column header; feature genes reselected on the symbol-keyed reference), same DepMap 26Q3 expression, units and rules.'] if GM is not None else []))
tmp=pd.read_csv(STATIC/'tableS_TMP_vs_DNA_surrogate_consistency.csv')
xl(TSTEM('tmp','.xlsx'),{'TMP':tmp},['HCMI Nature Supplementary Table S3 TMP predictions; UCEC TMP fields are blank in the source table.','The empty UCEC TMP fields are a limitation of the data source (the published HCMI table provides no UCEC TMP prediction); they are not missing values introduced by this analysis.'])


# ───────── unit-level DepMap release comparison (vs 24Q4) ─────────
R24=PIPE/'results'/('24Q4' if GID=='entrez' else '24Q4_symbol')
if not a.tag.startswith('24Q4') and (R24/'units_215_calls.csv').exists():
    K24=json.load(open(R24/'key_numbers.json'))
    u0=pd.read_csv(R24/'units_215_calls.csv'); u1=U
    dm=u1.dataset.isin(['DepMap_24Q4','DepMap_NextGen'])
    C=u0[dm][['dataset','model_unit','histology','stratum','is_org','nearest_A','strict_A','loose_A','nearest_B','strict_B','loose_B']].merge(
      u1[dm][['dataset','model_unit','nearest_A','strict_A','loose_A','nearest_B','strict_B','loose_B']],on=['dataset','model_unit'],suffixes=('_24Q4','_'+a.tag))
    in24=set(K24['depmap_qc'].get('unit_models_present_in_file',[])); innew=set(K['depmap_qc'].get('unit_models_present_in_file',[]))
    C['in_24Q4_expression_file']=C.model_unit.isin(in24); C[f'in_{a.tag}_default_rows']=C.model_unit.isin(innew)
    C['expression_24Q4_run']=np.where(C.in_24Q4_expression_file,'DepMap 24Q4',np.where(C.dataset=='DepMap_NextGen','NextGen manuscript file set','—'))
    C[f'expression_{a.tag}_run']=np.where(C[f'in_{a.tag}_default_rows'],'DepMap '+a.tag,'not in '+a.tag+' (previous values retained)')
    for k in ['nearest_A','strict_A','loose_A','nearest_B','strict_B','loose_B']:
        C['changed_'+k]=C[k+'_24Q4']!=C[k+'_'+a.tag]
    C['any_A_change']=C[['changed_nearest_A','changed_strict_A','changed_loose_A']].any(axis=1)
    C['resource']=C.dataset.map(DS_LABEL); C=relabel(C)
    Cpub=C.copy()
    Cpub=Cpub.drop(columns=[c for c in Cpub.columns if c.endswith('_B') or c.endswith('_B_24Q4') or '_B_' in c])   # Method B: supplementary Method B table only
    Cpub.to_csv(TAB/(TSTEM('release')+'.csv'),index=False)
    # Method B: the contrastive components come from the DepMap file, so every unit (not only DepMap units) can change
    CB=u0[['dataset','model_unit','stratum','is_org','nearest_B','strict_B','loose_B']].merge(u1[['dataset','model_unit','nearest_B','strict_B','loose_B']],on=['dataset','model_unit'],suffixes=('_24Q4','_'+a.tag))
    CB['changed_nearest_B']=CB['nearest_B_24Q4']!=CB['nearest_B_'+a.tag]; CB['resource']=CB.dataset.map(DS_LABEL); CB=relabel(CB)
    relsheets={'units':Cpub}
    xl(TSTEM('release','.xlsx'),relsheets,[f'All DepMap units (2D lines and NextGen organoids) among the 215; same reference (11 TCGA projects, full), same rules; only DepMap expression differs (24Q4 vs {a.tag}). dataset DepMap_2D = DepMap 2D lines.',
       f'in_24Q4_expression_file / in_{a.tag}_default_rows: model present in that release ({a.tag}: IsDefaultEntryForModel == Yes). NextGen organoids were taken from the NextGen manuscript file set in the 24Q4 run.',
       f'Gene matching: {K["depmap_qc"].get("gene_matching","stripped gene symbol")} in both runs (Methods 2.2).',
       f'Method B (Celligner-like) release changes for all 215 units are in {ST("methodB")} (supplementary only).',CODE_NOTE])
    chg=dict(n=len(C),nearest_A_changed=int(C.changed_nearest_A.sum()),strict_A_changed=int(C.changed_strict_A.sum()),loose_A_changed=int(C.changed_loose_A.sum()),
             strict_B_changed=(int(C.changed_strict_B.sum()) if BV else PENDING),changed_units=C[C.any_A_change][['model_unit','stratum','nearest_A_24Q4','nearest_A_'+a.tag,'strict_A_24Q4','strict_A_'+a.tag,'loose_A_24Q4','loose_A_'+a.tag]].to_dict('records'),
             only_in_24Q4=sorted(in24-innew),only_in_new=sorted((innew-in24)&set(C.model_unit)),
             nearest_B_changed_all_units=(int(CB.changed_nearest_B.sum()) if BV else PENDING),nearest_B_changed_by_resource=(CB[CB.changed_nearest_B].resource.value_counts().to_dict() if BV else PENDING))
    json.dump(chg,open(RES/'release_change_summary.json','w'),indent=1,default=str)
    print('release changes:',{k:v for k,v in chg.items() if k!='changed_units'})

# ───────── Method B workbook (round 4b: supplementary only; no conclusions drawn) ─────────
MB_STRATA=[('HGSC/serous (incl. serous NOS)','HGSC/serous (incl. serous NOS)'),('  OV_HGSC','  High-grade serous (HGSC)'),('  DepMap serous NOS','  DepMap serous NOS'),
           ('  HGSC/serous organoids','  HGSC/serous organoids'),('  HGSC/serous 2D lines','  HGSC/serous 2D lines'),('Endometrial carcinoma (non-CS)','Endometrial carcinoma (non-CS)'),
           ('  EC organoids','  EC organoids'),('  EC 2D lines','  EC 2D lines'),('Uterine carcinosarcoma','Uterine carcinosarcoma'),('Ovarian non-HGSC (coverage gap)','Ovarian non-HGSC (coverage gap)')]
MB_RUNS=[('',T2,U,f'{EXPR_LABEL}, {GENE_MATCH}')]
if GID=='entrez' and (RSYM/'table2_strict_loose.csv').exists(): MB_RUNS.append((' (symbol key)',Ts,pd.read_csv(RSYM/'units_215_calls.csv'),f'{EXPR_LABEL}, gene symbol'))
if not a.tag.startswith('24Q4') and (R24/'table2_strict_loose.csv').exists(): MB_RUNS.append((' (DepMap 24Q4)',pd.read_csv(R24/'table2_strict_loose.csv').set_index('stratum'),pd.read_csv(R24/'units_215_calls.csv'),'DepMap 24Q4, '+GENE_MATCH))
mbs=[]
for key,lab in MB_STRATA:
    r=T2.loc[key]; row={'Stratum':lab,'Units':int(r.n),'Strict A (primary; for comparison)':f'{int(r.strict_A)}/{int(r.n)}'}
    for suf,TT,_,_ in MB_RUNS:
        q=TT.loc[key]; row['Strict B'+suf]=f'{int(q.strict_B)}/{int(q.n)}'; row['Loose B'+suf]=f'{int(q.loose_B)}/{int(q.n)}'
    mbs.append(row)
MBS=pd.DataFrame(mbs)
kirc=[]
for suf,_,UU,desc in MB_RUNS:
    kk=UU[UU.nearest_B=='KIRC']; d={'Run':desc,'Units with KIRC nearest (Method B)':f'{len(kk)}/{len(UU)}'}
    d.update({DS_LABEL.get(k,k):int(v) for k,v in kk.dataset.value_counts().items()}); kirc.append(d)
MBK=pd.DataFrame(kirc).fillna(0)
for c in MBK.columns[2:]: MBK[c]=MBK[c].astype(int)
chg_rows=[]
if len(MB_RUNS)>1 and MB_RUNS[1][0]==' (symbol key)':
    us_=MB_RUNS[1][2]; m_=us_.set_index(['dataset','model_unit']).nearest_B.reindex(U.set_index(['dataset','model_unit']).index).values!=U.nearest_B.values
    chg_rows.append({'Comparison':'Gene key: symbol vs Entrez (same expression)','Nearest-centroid changes (Method B)':f'{int(m_.sum())}/215','Strict-call changes (Method B)':f"{int((us_.set_index(['dataset','model_unit']).strict_B.reindex(U.set_index(['dataset','model_unit']).index).values!=U.strict_B.values).sum())}/215",'By resource':json.dumps(U[m_].dataset.map(DS_LABEL).value_counts().to_dict())})
if 'CB' in globals():
    chg_rows.append({'Comparison':f'DepMap expression: 24Q4 vs {a.tag}','Nearest-centroid changes (Method B)':f"{int(CB.changed_nearest_B.sum())}/215",'Strict-call changes (Method B)':f"{int((CB['strict_B_24Q4']!=CB['strict_B_'+a.tag]).sum())}/215",'By resource':json.dumps(CB[CB.changed_nearest_B].resource.value_counts().to_dict())})
MBC=pd.DataFrame(chg_rows)
uu_=U[['dataset','model_unit','histology','stratum','is_org','nearest_B','strict_B','loose_B']].copy()
for suf,_,UU,_ in MB_RUNS[1:]:
    uu_=uu_.merge(UU[['dataset','model_unit','nearest_B']].rename(columns={'nearest_B':'nearest_B'+suf.replace(' ','_')}),on=['dataset','model_unit'],how='left')
uu_['resource']=uu_.dataset.map(DS_LABEL); uu_=relabel(uu_)
mbsheets={'Stratum_summary':MBS,'KIRC_nearest':MBK}
if len(MBC): mbsheets['Unit_changes']=MBC
mbsheets.update({'Confusion_histology':pd.read_csv(RES/'confusion_matrix_methodB.csv'),'Confusion_stratum':pd.read_csv(RES/'confusion_by_stratum_methodB.csv'),'Units_215':uu_})
if 'CB' in globals(): mbsheets['Release_24Q4_vs_'+a.tag]=CB
if (RES/'regression_test_real.json').exists():
    rt=json.load(open(RES/'regression_test_real.json')); rt=rt.get('real',rt); RT=pd.DataFrame(rt['table'])
    RT['C2_exceeded_methodB (>10/215)']=RT.changed_nearest_B>10; mbsheets['Stability_test']=RT
    if (RES/'regression_test_real_diagnostic.json').exists():
        DG=pd.DataFrame(json.load(open(RES/'regression_test_real_diagnostic.json'))); DG['changed_B_by_dataset']=DG.changed_B_by_dataset.map(lambda d: json.dumps({DS_LABEL.get(k,k)+(f' ({EXPR_LABEL} expression)' if k in ('DepMap_24Q4','DepMap_NextGen') else ''):v for k,v in d.items()})); mbsheets['Stability_diagnostic']=DG
_sg=pd.read_csv(cfg_path('sanger_organoid_rna'),index_col=0,usecols=[0]).index.str.replace('_exp$','',regex=True)
SANGER_TYPES=pd.read_excel(cfg_path('sanger_organoid_annotation'),'data_availability').set_index('sample_ID').reindex(_sg).primary_tumour_type.fillna('NA').value_counts()
SANGER_COMP=', '.join(f'{k.lower()} {v}' for k,v in SANGER_TYPES.items())
xl(TSTEM('methodB','.xlsx'),mbsheets,[
   'Method B (Celligner-like) is reported here as a supplementary analysis only; no conclusions are drawn from it (Methods 2.4). All conclusions rest on Method A.',
   'Method: tumour profiles centred on the mean of the 11 TCGA centroids, DepMap 2D lines on the lineage-balanced mean of all DepMap models, other models on the tumour-type-balanced Sanger organoid mean; 4 contrastive PCs (estimated on feature genes present in the DepMap file) projected out; Pearson nearest centroid. No mutual nearest-neighbour correction (not full Celligner).',
   'Stability (sheet Stability_test; pre-specified criteria): C1 smallest canonical correlation of contrastive subspaces >= 0.90; C2 at most 10/215 changed nearest-centroid calls. Dropping 5, 20 or 36 random feature genes (5 repeats each) changed 2-25/215 Method B calls with the current code, exceeding 10/215 in 9 of 15 runs: C2 NOT MET. path = earlier rows show the original gene-order defect (113-117/215 changes).',
   'Method B did not reproduce the organoid advantage within endometrial carcinoma seen with Method A (sheet Stratum_summary). KIRC nearest-centroid assignments (sheet KIRC_nearest) occur across all resources; their cause has not been established and they are not interpreted.',
   'Column suffixes: no suffix = headline expression and gene key; (symbol key) = gene-symbol sensitivity; (DepMap 24Q4) = DepMap 24Q4 expression sensitivity.',
   f'Centring reference for all non-DepMap models (HCMI, NextGen, public organoids): tumour-type-balanced mean of the {int(SANGER_TYPES.sum())} Sanger organoids in the expression file ({SANGER_COMP}); it contains no endometrial organoids.',
   'Method B unit-level calls must not be used to select models. The endometrial organoid advantage over 2D lines holds only under Method A; no biological interpretation is given for the KIRC assignments.',
   'Sheet Stability_diagnostic: resource keys as in the regression-test output; DepMap 2D lines and NextGen organoids carry the expression release of the run (labelled).',
   CODE_NOTE])

# ───────── number-audit with manuscript line numbers ─────────
ms_all=open(a.manuscript,encoding='utf-8').read().split('\n') if a.manuscript and Path(a.manuscript).exists() else []
# audit searches the manuscript body only (title to Tables; the reference list and data citations are excluded)
body_end=next((i for i,l in enumerate(ms_all) if l.strip()=='## References'),len(ms_all))
ms_lines=ms_all[:body_end]
def find(pats,ctx=None,nctx=None):
    hits=[]
    for i,l in enumerate(ms_lines,1):
        if any(p in l for p in pats) and (ctx is None or any(c in l for c in ([ctx] if isinstance(ctx,str) else ctx))) and (nctx is None or nctx not in l): hits.append(i)
    return ';'.join(map(str,hits))
rows=[]
_XL={}
def _sheet(tkey,sheet):
    k=(tkey,sheet)
    if k not in _XL: _XL[k]=pd.read_excel(TAB/TSTEM(tkey,'.xlsx'),sheet_name=sheet)
    return _XL[k]
def sloc(tkey,sheet,col,where,expect):
    """locate a value in a supplementary workbook; returns 'file | sheet | Excel row(s) | column' after verifying the cell content"""
    try:
        df=_sheet(tkey,sheet); m=pd.Series(True,index=df.index)
        for c,v in (where or {}).items():
            if c=='_nonempty': m&=df[v].fillna('').astype(str).str.strip().str.len()>0
            else: m&=df[c].astype(str).str.strip()==str(v).strip()
        idx=list(df.index[m])
        if callable(expect): ok=expect(df.loc[idx])
        else:
            assert len(idx)==1; cell=df.at[idx[0],col]
            ok=(abs(float(cell)-float(expect))<1e-6) if isinstance(expect,(int,float)) else (str(expect) in str(cell))
        if not ok or not idx: return None
        rr=(f'row {idx[0]+2}' if len(idx)==1 else 'rows '+','.join(str(i+2) for i in idx))
        return f"tables/{TSTEM(tkey,'.xlsx')} | sheet {sheet} | {rr} | column {col}"
    except Exception as e:
        print('sloc failed',tkey,sheet,col,where,e); return None
def add(claim,value,pats,src,note='',ctx=None,nctx=None,ev=None,supp=None):
    ln=find(pats,ctx,nctx) if pats else ''
    loc=sloc(*supp) if (supp and not ln) else ''
    st_=('FOUND' if ln else ('SUPP-LOCATED' if loc else 'NOT FOUND'))
    rows.append(dict(claim=claim,value=value,manuscript_lines=ln,status=st_,supplement_location=loc or '',source_file=src,notes=note,expression_version=ev or EXPR_LABEL))
T2KEY={'HGSC_serous':'HGSC/serous (incl. serous NOS)','OV_HGSC':'  OV_HGSC','serous_NOS':'  DepMap serous NOS','EC_nonCS':'Endometrial carcinoma (non-CS)','EC_organoid':'  EC organoids',
       'EC_2D':'  EC 2D lines','CS':'Uterine carcinosarcoma','nonHGSC':'Ovarian non-HGSC (coverage gap)','HGSC_organoid':'  HGSC/serous organoids','HGSC_2D':'  HGSC/serous 2D lines'}
MBLAB={k:dict(MB_STRATA)[v] for k,v in T2KEY.items()}
def t2loc(cl,col,val): return ('confusion','Table2_full',col,{'stratum':T2KEY[cl]},val)
def mbloc(cl,col,val): return ('methodB','Stratum_summary',col,{'Stratum':MBLAB[cl]},val)
NOT24='DepMap 24Q4 expression'
for key,cl in [('HGSC/serous (incl. serous NOS)','HGSC_serous'),('  OV_HGSC','OV_HGSC'),('  DepMap serous NOS','serous_NOS'),('Endometrial carcinoma (non-CS)','EC_nonCS'),
               ('  EC organoids','EC_organoid'),('  EC 2D lines','EC_2D'),('Uterine carcinosarcoma','CS'),('Ovarian non-HGSC (coverage gap)','nonHGSC'),
               ('  HGSC/serous organoids','HGSC_organoid'),('  HGSC/serous 2D lines','HGSC_2D')]:
    r=T2.loc[key]
    add(cl+'_strict_A',f'{int(r.strict_A)}/{int(r.n)} ({r.strict_A_pct:.1f}%)',[f'{int(r.strict_A)}/{int(r.n)}'],f'source_data/{a.tag}/table2_strict_loose.csv',nctx=NOT24,supp=t2loc(cl,'strict_A',int(r.strict_A)))
    add(cl+'_loose_A',f'{int(r.loose_A)}/{int(r.n)} ({r.loose_A_pct:.1f}%)',[f'{int(r.loose_A)}/{int(r.n)}'] if cl not in ('HGSC_organoid','HGSC_2D') else None,f'source_data/{a.tag}/table2_strict_loose.csv',nctx=NOT24,supp=t2loc(cl,'loose_A',int(r.loose_A)))
    if cl=='CS':
        add('CS_organoid_strict_A',f'{int(r.organoid_strict_A)}/{int(r.organoid_n)}',[f'{int(r.organoid_strict_A)}/{int(r.organoid_n)}'],f'source_data/{a.tag}/table2_strict_loose.csv')
        add('CS_2D_strict_A',f'{int(r.twoD_strict_A)}/{int(r.twoD_n)}',[f'{int(r.twoD_strict_A)}/{int(r.twoD_n)}'],f'source_data/{a.tag}/table2_strict_loose.csv')
for key,cl in [('HGSC/serous (incl. serous NOS)','HGSC_serous'),('Endometrial carcinoma (non-CS)','EC_nonCS'),('  EC organoids','EC_organoid'),('  EC 2D lines','EC_2D'),('Uterine carcinosarcoma','CS')]:
    r=T2.loc[key]
    for rule in ['strict','loose']:
        v=f"{int(r[rule+'_B'])}/{int(r.n)}"
        add(f'{cl}_{rule}_B',v,None,f'source_data/{a.tag}/table2_strict_loose.csv','Method B (supplementary only)',supp=mbloc(cl,('Strict B' if rule=='strict' else 'Loose B'),v))
if BV:
    v=f"{int((U.nearest_B=='KIRC').sum())}/{len(U)}"
    add('nearest_B_KIRC',v,[v],f'source_data/{a.tag}/units_215_calls.csv','by resource: '+json.dumps(U[U.nearest_B=='KIRC'].dataset.replace(DS_PUBLIC).value_counts().to_dict()),ctx='KIRC',supp=('methodB','KIRC_nearest','Units with KIRC nearest (Method B)',{'Run':MB_RUNS[0][3]},v))
ex=K['extra']
add('EC_2D_loose_unit_IDs',ex['EC2D_loose_ids'],None,f'source_data/{a.tag}/EC_2D_units_calls.csv',f"n={len(ex['EC2D_loose_ids'].split(';'))}",supp=('confusion','EC_2D_units','model_unit',{'loose_A':True},lambda d:set(d.model_unit)==set(ex['EC2D_loose_ids'].split(';'))))
add('EC_2D_strict_unit_IDs',ex['EC2D_strict_ids'],None,f'source_data/{a.tag}/EC_2D_units_calls.csv',f"n={len(ex['EC2D_strict_ids'].split(';'))}",supp=('confusion','EC_2D_units','model_unit',{'strict_A':True},lambda d:set(d.model_unit)==set(ex['EC2D_strict_ids'].split(';'))))
ec2=pd.read_csv(RES/'EC_2D_units_calls.csv'); nU=int((ec2.nearest_A=='UCS').sum())
add('EC_2D_nearest_UCS',f'{nU}/{len(ec2)}',[f'nearest to UCS ({nU}/{len(ec2)})',f'UCS ({nU}/{len(ec2)})'],f'source_data/{a.tag}/EC_2D_units_calls.csv')
csn=ex['CS_nearest']; add('CS_nearest_counts',json.dumps(csn),[f"{csn.get('UCS',0)} nearest UCS"],f'source_data/{a.tag}/units_215_calls.csv')
sn=U[U.stratum=='serous_NOS']; v=f"{int((sn.nearest_A=='UCS').sum())}/{len(sn)}"
add('serous_NOS_nearest_UCS',v,[f'{v} nearest UCS'],f'source_data/{a.tag}/units_215_calls.csv')
add('ovary_nearest_UCS',f"{ex['ovary_nearest_UCS']}/{ex['ovary_n']}",[f"{ex['ovary_nearest_UCS']}/{ex['ovary_n']}"],f'source_data/{a.tag}/units_215_calls.csv')
add('ovary_nearest_UCEC',f"{ex['ovary_nearest_UCEC']}/{ex['ovary_n']}",[f"{ex['ovary_nearest_UCEC']}/{ex['ovary_n']}"],f'source_data/{a.tag}/units_215_calls.csv')
pu=100*(ex['ovary_nearest_UCS']+ex['ovary_nearest_UCEC'])/ex['ovary_n']; add('ovary_nearest_uterine_pct',f'{pu:.1f}%',[f'{pu:.1f}%'],f'source_data/{a.tag}/units_215_calls.csv')
add('CPTAC_OV_nearest_OV',f"{cp['OV_nearestOV']}/{cp['OV_n']}",[f"{cp['OV_nearestOV']}/{cp['OV_n']}"],f'source_data/{a.tag}/key_numbers.json')
add('CPTAC_UCEC_nearest_UCEC_or_UCS',f"{cp['UCEC_nearest_UCEC_or_UCS']}/{cp['UCEC_n']}",[f"{cp['UCEC_nearest_UCEC_or_UCS']}/{cp['UCEC_n']}"],f'source_data/{a.tag}/key_numbers.json')
for k,v in K['loo_A'].items():
    if k in ('OV','UCEC','UCS'): add(f'TCGA_LOO_A_{k}',f'{v:.3f}',[f'{k} {v:.3f}'],f'source_data/{a.tag}/tcga_loo_methodA.csv')
for _,r in ES.iterrows():
    if r.stratum in ('HGSC_serous','EC_organoid','EC_2D','EC_nonCS') :
        v=f'{r.mean_pct:.1f}% ({r.min_pct:.1f}–{r.max_pct:.1f}%)'
        add(f'equal_depth_{r.stratum}_{r.rule}',v,([f'{r.mean_pct:.1f}% (range {r.min_pct:.1f}', f'{r.mean_pct:.1f}% ({r.min_pct:.1f}'] if r.rule=='strict' and r.stratum!='EC_nonCS' else None),f'source_data/{a.tag}/equal_depth_summary.csv','' if r.rule=='strict' and r.stratum!='EC_nonCS' else f'{ST("sensitivity")} / {SF("equal_depth")} only',supp=('sensitivity','equal_depth_summary','mean_pct',{'stratum':r.stratum,'rule':r.rule},float(r.mean_pct)))
for _,r in SE.iterrows():
    for s_ in ['HGSC_serous','EC_nonCS']:
        add(f'refcomp[{r.reference}]_{s_}_strict',f"{r[s_+'_strict']:.1f}/{int(r[s_+'_n'])}",None,f'source_data/{a.tag}/reference_composition_sensitivity.csv',supp=('sensitivity','reference_composition',s_+'_strict',{'reference':r.reference},float(r[s_+'_strict'])))
# reference-composition claims quoted in text (context: the sensitivity paragraph)
for ref,col,lab,ctx in [('Gyn full + non-gyn 40/project','nonHGSC_loose','nonHGSC_loose','40 tumours/project'),('Gyn ~200/project + non-gyn full','HGSC_serous_strict','HGSC_serous_strict','subsample'),
                        ('Gyn ~200/project + non-gyn full','EC_nonCS_strict','EC_nonCS_strict','subsample'),('Gyn ~200/project + non-gyn full','EC_2D_strict','EC_2D_strict','subsample'),
                        ('Gyn full + non-gyn 40/project','CS_loose','CS_loose','carcinosarcoma loose')]:
    r=SE[SE.reference==ref].iloc[0]; n=int(r[col.rsplit('_',1)[0]+'_n']); v=f'{int(r[col])}/{n}'
    add(f'refcomp[{ref}]_{lab}',v,[v],f'source_data/{a.tag}/reference_composition_sensitivity.csv',ctx=ctx,supp=('sensitivity','reference_composition',col,{'reference':ref},float(r[col])))
# Method B centring reference for non-DepMap models (Sanger organoid tumour types; from the Sanger annotation input)
add('MethodB_Sanger_centring_organoids',str(int(SANGER_TYPES.sum())),[f'{int(SANGER_TYPES.sum())} Sanger organoids'],'input: Sanger data_RNAseq_ALL-organoids.csv + supplementary_table_2_revision.xlsx (data_availability)',ctx='Method B',supp=None)
for k_,v_ in SANGER_TYPES.items():
    add(f'MethodB_Sanger_centring_{k_}',str(int(v_)),[f'{k_.lower()} {int(v_)}'],'input: Sanger data_RNAseq_ALL-organoids.csv + supplementary_table_2_revision.xlsx (data_availability)',ctx='Method B')
loo=K['loo_A']; ng={k:v for k,v in loo.items() if k not in ('OV','UCEC','UCS')}
v=f"{min(ng.values()):.3f}–{max(ng.values()):.3f}"; add('TCGA_LOO_A_nongyn_range',v,[v],f'source_data/{a.tag}/tcga_loo_methodA.csv')
# DepMap 24Q4 sensitivity (same reference/rules/units): matched only on lines that state the 24Q4 sensitivity
p24=R24/'table2_strict_loose.csv'
if not a.tag.startswith('24Q4') and p24.exists():
    T24=pd.read_csv(p24).set_index('stratum')
    for key,cl in [('HGSC/serous (incl. serous NOS)','HGSC_serous'),('Endometrial carcinoma (non-CS)','EC_nonCS'),('  EC 2D lines','EC_2D'),('Uterine carcinosarcoma','CS')]:
        r=T24.loc[key]
        for rule in ['strict','loose']:
            v=f"{int(r[rule+'_A'])}/{int(r.n)}"; add(f'DepMap24Q4_{cl}_{rule}_A',v,[v],'source_data/24Q4/table2_strict_loose.csv','sensitivity (DepMap 24Q4 expression)',ctx=NOT24,ev='DepMap 24Q4')
    # Method B under DepMap 24Q4 (valid run) quoted in the text
    if json.load(open(R24/'key_numbers.json')).get('methodB_qc',{}).get('gene_order_fixed'):
        for key,cl in [('HGSC/serous (incl. serous NOS)','HGSC_serous'),('Endometrial carcinoma (non-CS)','EC_nonCS'),('  EC organoids','EC_organoid')]:
            r=T24.loc[key]; v=f"{int(r.strict_B)}/{int(r.n)}"
            add(f'DepMap24Q4_{cl}_strict_B',v,None,'source_data/24Q4/table2_strict_loose.csv','Method B, DepMap 24Q4 (supplementary only)',ev='DepMap 24Q4',supp=mbloc(cl,'Strict B (DepMap 24Q4)',v))
    if 'chg' in globals():
        v=f"{chg['nearest_A_changed']} of {chg['n']}"; add('release_nearest_A_changed',v,[v],f"tables/{TSTEM('release')}.csv",'; '.join(f"{d['model_unit']} {d['nearest_A_24Q4']}->{d['nearest_A_'+a.tag]}" for d in chg['changed_units']),ev=f'DepMap 24Q4 vs {a.tag}')
# round 4: gene-key sensitivity, Method B organoid/2D, Method B release changes, missingness
def _miss(d):
    m=re.search(r'missing-gene fraction after swap: ([0-9.]+)',open(d/'run_log.txt',encoding='utf-8').read()) if (d/'run_log.txt').exists() else None
    return float(m.group(1)) if m else None
for key,cl in [('  HGSC/serous organoids','HGSC_organoid'),('  HGSC/serous 2D lines','HGSC_2D')]:
    if BV:
        r=T2.loc[key]; v=f"{int(r.strict_B)}/{int(r.n)}"; add(f'{cl}_strict_B',v,None,f'source_data/{a.tag}/table2_strict_loose.csv','Method B (supplementary only)',supp=mbloc(cl,'Strict B',v))
mE=_miss(RES)
if mE is not None: add('DepMap_missing_gene_fraction_entrez',f'{100*mE:.2f}%',[f'{100*mE:.2f}%'],f'source_data/{a.tag}/run_log.txt','DepMap 2D profiles, shared gene space')
if GID=='entrez' and (RSYM/'table2_strict_loose.csv').exists():
    mS=_miss(RSYM)
    if mS is not None: add('DepMap_missing_gene_fraction_symbol',f'{100*mS:.2f}%',[f'{100*mS:.2f}%'],f'source_data/{RSYM.name}/run_log.txt','symbol matching (sensitivity)',ev=EXPR_LABEL+' (symbol)')
    add('feature_genes_absent_symbol_vs_entrez',f"{Ks.get('feature_genes_absent_from_depmap')} vs {K.get('feature_genes_absent_from_depmap')}",[f"36 of the 2,000 feature genes"],f'source_data/{RSYM.name}/key_numbers.json',ev=EXPR_LABEL)
    for key,cl,col in [('Endometrial carcinoma (non-CS)','EC_nonCS','strict_A'),('  EC 2D lines','EC_2D','strict_A'),('HGSC/serous (incl. serous NOS)','HGSC_serous','strict_B'),('Endometrial carcinoma (non-CS)','EC_nonCS','strict_B')]:
        y=Ts.loc[key]; v=f"{int(y[col])}/{int(y.n)}"; add(f'gene_key_symbol_{cl}_{col}',v,[v] if col.endswith('_A') else None,f'source_data/{RSYM.name}/table2_strict_loose.csv','symbol matching (sensitivity)',ctx='gene symbols',ev=EXPR_LABEL+' (symbol)',supp=(('sensitivity','gene_matching_symbol','Strict A (symbol)',{'Stratum':key.strip().replace('OV_HGSC','HGSC')},v) if col.endswith('_A') else mbloc(cl,'Strict B (symbol key)',v)))
    v=f"{ch['nearest_B_changed']}/215"; add('gene_key_symbol_nearest_B_changed',v,None,f'source_data/{a.tag}/gene_matching_sensitivity.json','Method B (supplementary only)',ev=EXPR_LABEL+' (symbol vs Entrez)',supp=('methodB','Unit_changes','Nearest-centroid changes (Method B)',{'Comparison':'Gene key: symbol vs Entrez (same expression)'},v))
    add('gene_key_symbol_nearest_A_changed',', '.join(ch['units_nearest_A_changed']),[ch['units_nearest_A_changed'][0]] if ch['units_nearest_A_changed'] else None,f'source_data/{a.tag}/gene_matching_sensitivity.json',ctx='gene symbols',ev=EXPR_LABEL+' (symbol vs Entrez)')
if 'chg' in globals() and BV and isinstance(chg.get('nearest_B_changed_all_units'),int):
    v=f"{chg['nearest_B_changed_all_units']}/215"; add('release_nearest_B_changed_all_units',v,None,f'source_data/{a.tag}/release_change_summary.json',json.dumps(chg['nearest_B_changed_by_resource'])+'; Method B (supplementary only)',ev=f'DepMap 24Q4 vs {a.tag}',supp=('methodB','Unit_changes','Nearest-centroid changes (Method B)',{'Comparison':f'DepMap expression: 24Q4 vs {a.tag}'},v))
    v24=f"{int((u0.nearest_B=='KIRC').sum())}/215"; add('DepMap24Q4_nearest_B_KIRC',v24,None,'source_data/24Q4/units_215_calls.csv','Method B, DepMap 24Q4 (supplementary only)',ev='DepMap 24Q4',supp=('methodB','KIRC_nearest','Units with KIRC nearest (Method B)',{'Run':'DepMap 24Q4, '+GENE_MATCH},v24))
if (RES/'regression_test_real.json').exists():
    _rt=RT.copy(); cur=_rt[_rt.path=='current']; leg=_rt[_rt.path!='current']
    v=f"{int(cur.changed_nearest_B.min())}–{int(cur.changed_nearest_B.max())}/215"; add('regression_methodB_changed_range',v,[v],f'source_data/{a.tag}/regression_test_real.json','pre-specified C2: <=10/215',supp=('methodB','Stability_test','changed_nearest_B',{'path':'current'},lambda d:f"{int(d.changed_nearest_B.min())}–{int(d.changed_nearest_B.max())}/215"==v))
    v=f"{int((cur.changed_nearest_B>10).sum())} of {len(cur)}"; add('regression_methodB_C2_exceeded',v,[v],f'source_data/{a.tag}/regression_test_real.json','C2 not met')
    v=f"{int(cur.changed_nearest_A.max())}/215"; add('regression_methodA_max_changed',v,[f'at most {v}'],f'source_data/{a.tag}/regression_test_real.json')
    v=f"{leg.cc.min():.3f}–{leg.cc.max():.3f}"; add('regression_legacy_cc_range',v,[v],f'source_data/{a.tag}/regression_test_real.json','original (legacy) code path')
    v=f"{int(leg.changed_nearest_B.min())}–{int(leg.changed_nearest_B.max())}/215"; add('regression_legacy_methodB_changed',v,[v],f'source_data/{a.tag}/regression_test_real.json','original (legacy) code path')
add('identity_threshold_CN',f'{TH_C:.3f}',[f'r > {TH_C:.3f}'],f'source_data/{a.tag}/identity_identity_summary.json',ev='not applicable (DNA)')
if (RES/'histotype_reference_checks.json').exists():
    _hr=json.load(open(RES/'histotype_reference_checks.json'))
    for k_,lab_ in [('GSE295399_LOO','GSE295399 LOO'),('GSE226870_external','GSE226870'),('GSE121103_external','GSE121103')]:
        add(f'histotype_reference_{k_}',_hr[k_],[f'{lab_} {_hr[k_]}'],f'source_data/{a.tag}/histotype_reference_checks.json','03_histotype_reference_checks.py',ev='not applicable (tumour references)')
h=K['histotype']
for col,lab in [('hist_uncentred','uncentred'),('hist_ordinary','ordinary_centring'),('hist_label_balanced','label_balanced')]:
    add(f'histotype_{lab}_nonHGSC',h[col]['nonHGSC'],[h[col]['nonHGSC']],f'source_data/{a.tag}/units_215_calls.csv')
    add(f'histotype_{lab}_HGSC',h[col]['HGSC'],[h[col]['HGSC']],f'source_data/{a.tag}/units_215_calls.csv')
pp=pd.read_csv(RES/'histotype_call_PPV.csv')
for _,r in pp.iterrows(): add(f'histotype_PPV_{r.called_as}',f'{int(r.true_positives)}/{int(r.n_called)}',[f'{int(r.true_positives)}/{int(r.n_called)}'] if r.called_as!='LGSC' else None,f'source_data/{a.tag}/histotype_call_PPV.csv',supp=('histotype','PPV_uncentred','true_positives',{'called_as':r.called_as},lambda d,r=r:int(d.true_positives.iloc[0])==int(r.true_positives) and int(d.n_called.iloc[0])==int(r.n_called)))
mk=K['markers']
add('WT1_low_HGSC_serous',mk['WT1_low_HGSC_serous'],[mk['WT1_low_HGSC_serous']],f'source_data/{a.tag}/units_215_calls.csv')
add('mucinous_GI_like',mk['mucinous_GI_like_n']+' ('+mk['mucinous_GI_like']+')',[mk['mucinous_GI_like_n'].split('/')[0]+'/8 GI-like'],f'source_data/{a.tag}/units_215_calls.csv')
add('paired_RNA_own_vs_other_median',f"{paired['rna_own_median']:.3f} vs {paired['rna_other_median']:.3f}",[f"{paired['rna_own_median']:.3f}"],'source_data/static/rna_paired_purityfilt.csv',ev='HCMI open RNA (not DepMap)')
add('CPTAC_DNA_surrogate',f'{int((Ev.rule_meth==Ev.truth).sum())}/{len(Ev)}',[f'{int((Ev.rule_meth==Ev.truth).sum())}/{len(Ev)}'],'source_data/static/cptac_rule_evaluation.csv',ev='not applicable (DNA)')
add('uterine_DNA_classes',' '.join(f'{c}{n}' for c,n in zip(t3.Class,t3.iloc[:,1])),['MMRd 23'],'tables/Table3_endometrial_classes.xlsx','n=51 incl. DNA-only ACH-003239',ev='not applicable (DNA; MLH1 expression criterion unchanged between releases)')
add('organoid_DNA_classes',' '.join(f'{c}{n}' for c,n in zip(t3.Class,t3.iloc[:,2])),['p53abn 8'],'tables/Table3_endometrial_classes.xlsx','n=19',ev='not applicable (DNA)')
add('TCGA_reference_counts',json.dumps(K['reference_counts']),['GDC Data Release'],'source_data/static/nongyn_full_manifest.csv; gyn_full_manifest.csv',K['gdc_release'],ev='not applicable (TCGA)')
add('independent_units','215 (126 organoid / 89 2D)',['**215**'],'source_data/static/units_215.csv',ev='not applicable')
add(f'{TSTEM("selection")}_ATCC_filled',str(s2_stats['ATCC_filled']),None,f"tables/{TSTEM('selection')}.csv",f"any repository catalogue {s2_stats['any_catalog']}; units with ACH {s2_stats['units_with_ACH']}",ev='not applicable (metadata)',supp=('selection','Selection','ATCC_catalog',{'_nonempty':'ATCC_catalog'},lambda d:len(d)==s2_stats['ATCC_filled']))
# identity / crosswalk claims, recomputed from the package files (no static copies)
XW=pd.read_csv(XWALK/'crosswalk_hcmi_ach_sidm.csv'); XS=open(XWALK/'summary_crosswalk.txt',encoding='utf-8').read()
assert len(XW)==141 and 'UNAVAILABLE' not in XS and set(XW.source_key)=={'DepMap_26Q3_Model.csv'}, 'crosswalk is not the 26Q3-merged version'
off=XW[XW.official_hcmi_from_nextgen_supp.notna()]
add('NextGen_official_HCMI_pairs',str(len(off)),['five official HCMI-linked','Five of 19 NextGen'],'source_data/static/crosswalk_hcmi_ach_sidm.csv (official_hcmi_from_nextgen_supp)','NextGen Nature Supplementary Table 1 Cell Line Name',ev='not applicable (metadata)')
IDc=I.copy(); called_ii=set(IDc.nextgen[IDc.crit_i|IDc.crit_ii])
# criterion (iii) (count above the HCMI between-case maximum and case-specific; Methods 2.3): ACH-003709 (case-specific excess) yes, ACH-003903 (hypermutation, low Jaccard) no
called=called_ii|{'ACH-003709'}
fp=sorted(called-set(off.ACH)); rec=len(set(off.ACH)&called)
add('multiomic_validation',f'{rec}/{len(off)} recovered; {len(fp)} false positives',['5/5'],f'source_data/{a.tag}/identity_nextgen_identity_check.csv',f"by (i)/(ii): {sorted(called_ii)}; ACH-003709 by criterion (iii); ACH-003903 not called",ev=IDS['expression_source'])
for nm in ['ACH-003709','ACH-003903']:
    r_=_i9.loc[nm]; v=f"{int(r_.shared_damaging_genes)} shared; Jaccard {r_.jaccard:.3f}"
    add(f'identity_{nm}_iii',v,[f"{int(r_.shared_damaging_genes)} shared truncating genes" if nm=='ACH-003709' else f"shared {int(r_.shared_damaging_genes)} truncating genes"],f'source_data/{a.tag}/identity_nextgen_identity_check.csv',ev='not applicable (DNA)')
add('identity_thresholds',f'expr {TH_E:.3f}; CN {TH_C:.3f}; mut {TH_M}',[f'r > {TH_E:.3f}'],f'source_data/{a.tag}/identity_identity_summary.json',ev=IDS['expression_source'])
tm=tmp[tmp.cohort=='UCEC']; add('TMP_UCEC_available',f"{int(tm.tmp_model_subtype.notna().sum())}/{len(tm)}",None,'source_data/static/tableS_TMP_vs_DNA_surrogate_consistency.csv',ev='not applicable',supp=('tmp','TMP','tmp_model_subtype',{'cohort':'UCEC'},lambda d:len(d)==len(tm) and int(d.tmp_model_subtype.notna().sum())==int(tm.tmp_model_subtype.notna().sum())))
add('DepMap_26Q3_Model_csv',f'{len(M26)} models; HCMIID non-null {int(M26.HCMIID.notna().sum())}',['2,554 models'],'DepMap Public 26Q3 Model.csv (downloaded 2026-10-07)',ev='not applicable (metadata)')
add('gyn_crosswalk_rows_26Q3',str(len(XW)),[f'{len(XW)} gynaecological DepMap models'],'source_data/static/crosswalk_hcmi_ach_sidm.csv','Ovary/Fallopian/Uterus',ev='not applicable (metadata)')
cl_=XW[XW.ModelType=='Cell Line']; og_=XW[XW.ModelType=='Organoid']
add('gyn_cell_lines_and_organoids_26Q3',f'{len(cl_)} cell lines + {len(og_)} organoids',[f'{len(cl_)} cell lines + {len(og_)} organoids'],'source_data/static/crosswalk_hcmi_ach_sidm.csv',ev='not applicable (metadata)')
v=f'{int(cl_.SIDM.notna().sum())}/{len(cl_)}'; add('gyn_cell_SIDM_26Q3',v,[v],'source_data/static/crosswalk_hcmi_ach_sidm.csv',ev='not applicable (metadata)')
hc=XW[XW.HCMIID.notna()]; add('gyn_HCMIID_nonnull_26Q3',str(len(hc)),['six HCMIID-bearing'],'source_data/static/crosswalk_hcmi_ach_sidm.csv','5 NextGen official + ACH-002419',ev='not applicable (metadata)')
add('ACH_SIDM_HCMIID_triples_26Q3',str(int((XW.SIDM.notna()&XW.HCMIID.notna()).sum())),['triples remain 0'],'source_data/static/crosswalk_hcmi_ach_sidm.csv',ev='not applicable (metadata)')
m5=int(sum(1 for r in off.itertuples() if isinstance(r.HCMIID,str) and str(r.HCMIID)==str(r.official_hcmi_from_nextgen_supp)))
add('NextGen5_in_26Q3_with_matching_HCMIID',f'{m5}/{len(off)}',['all five carry matching'],'source_data/static/crosswalk_hcmi_ach_sidm.csv',ev='not applicable (metadata)')
e903=XW[XW.ACH=='ACH-003903']; add('ACH-003903_in_26Q3_HCMIID_empty',str(bool(len(e903)) and bool(e903.HCMIID.isna().all())),['in 26Q3 with empty HCMIID'],'source_data/static/crosswalk_hcmi_ach_sidm.csv',ev='not applicable (metadata)')
NA=pd.DataFrame(rows); NA.insert(0,'row',range(1,len(NA)+1))
NA=NA[['row','claim','value','manuscript_lines','status','supplement_location','source_file','notes','expression_version']]
NA.to_csv(RES/'number-audit.csv',index=False)
json.dump(dict(s2=s2_stats,paired=paired,methodB_valid=BV),open(RES/'outputs_summary.json','w'),indent=1,default=str)
print('number-audit rows',len(NA),NA.status.value_counts().to_dict(),'NOT FOUND:',int((NA.status=='NOT FOUND').sum()))
print(NA[NA.status=='NOT FOUND'][['claim','value']].to_string(index=False))

# ───────── publish ─────────
# strip directory parts of absolute paths (keep file names) and internal stage names from logs/json
INTERNAL=[(re.compile(r'(?<![\w.])/(?:[\w.\-]+/)+'),''),(re.compile(r'DepMap 2D vs \w+ harmonised'),'DepMap 2D vs harmonised')]
def clean_text(t):
    for rx,rep in INTERNAL: t=rx.sub(rep,t)
    return t
def copy_clean(src,dst,mask_B=False):
    """copy a result file into the package: internal paths removed, internal dataset key relabelled, Method B withheld if invalid"""
    if src.suffix=='.csv':
        df=pd.read_csv(src)
        if 'local_path' in df.columns: df['local_path']=df.local_path.astype(str).map(lambda p: Path(p).name)
        if 'already_local' in df.columns: df=df.drop(columns='already_local')   # download bookkeeping of the study machine, not data (round 5d)
        df=relabel(df)
        if mask_B:
            for c in df.columns:
                if c in ('nearest_B','strict_B','loose_B'): df[c]=PENDING
        df.to_csv(dst,index=False)
    elif src.suffix in ('.json','.txt'):
        t=clean_text(open(src,encoding='utf-8').read())
        if mask_B and src.name=='run_log.txt': t='NOTE: Method B lines in this log are superseded ('+BV_NOTE+')\n'+t
        if mask_B and src.name=='release_change_summary.json':
            d=json.loads(t); d['strict_B_changed']=PENDING; t=json.dumps(d,indent=1)
        open(dst,'w',encoding='utf-8').write(t)
    else: shutil.copy2(src,dst)
def archive_stale(folder,keep_names,arch):
    arch.mkdir(parents=True,exist_ok=True)
    for f_ in folder.glob('*.*'):
        if f_.name not in keep_names: shutil.move(str(f_),str(arch/f_.name))
def internal_audit(folder,base,tag):
    # audits of sensitivity runs check the headline text against sensitivity values; kept internally, not in the public package
    f_=folder/'number-audit.csv'
    if f_.exists():
        d_=base/'archive/internal_audit'/tag; d_.mkdir(parents=True,exist_ok=True); shutil.move(str(f_),str(d_/'number-audit.csv'))
if a.publish:
    assert int((NA.status=='NOT FOUND').sum())==0, 'public number-audit must have zero NOT FOUND rows'
    for base in [b for b in [MS,DELIV] if b]:
        (base/'figures').mkdir(exist_ok=True)
        archive_stale(base/'figures',{p.name for p in FIG.glob('*')},base/'figures/archive')
        for f_ in FIG.glob('*'): shutil.copy2(f_,base/'figures'/f_.name)
        (base/'tables').mkdir(exist_ok=True)
        archive_stale(base/'tables',{p.name for p in TAB.glob('*')},base/'tables/archive')
        for f_ in TAB.glob('*'): shutil.copy2(f_,base/'tables'/f_.name)
        sd=base/'supp/source_data'/a.tag; sd.mkdir(parents=True,exist_ok=True)
        res_files=[f_ for f_ in RES.glob('*.*') if BV or 'methodB' not in f_.name or f_.name=='methodB_validity.json']
        archive_stale(sd,{f_.name for f_ in res_files},base/'supp/source_data/archive'/a.tag)
        for f_ in res_files: copy_clean(f_,sd/f_.name,mask_B=not BV)
        for f_ in IDDIR.glob('*.*'): copy_clean(f_,sd/('identity_'+f_.name))
        if RSYM.exists():   # gene-symbol sensitivity run (source data)
            dS=base/'supp/source_data'/RSYM.name; dS.mkdir(parents=True,exist_ok=True)
            internal_audit(dS,base,RSYM.name)
            archive_stale(dS,{f_.name for f_ in RSYM.glob('*.*') if f_.name!='number-audit.csv'}|{'identity_'+f_.name for f_ in (RSYM/'identity').glob('*.*')},base/'supp/source_data/archive'/RSYM.name)
            for f_ in RSYM.glob('*.*'):
                if f_.name!='number-audit.csv': copy_clean(f_,dS/f_.name)
            for f_ in (RSYM/'identity').glob('*.*'): copy_clean(f_,dS/('identity_'+f_.name))
        st_=base/'supp/source_data/static'; st_.mkdir(parents=True,exist_ok=True)
        for nm_ in ['nextgen_identity_check.csv','nextgen_identity_summary.txt',   # superseded by source_data/<tag>/identity_* (round 4)
                    'uterine_units_plus_nextgen.csv']:   # superseded by uterine_units_classes.csv (05_derive_endometrial_classes.py; round 5c)
            if (st_/nm_).exists(): (base/'supp/source_data/archive/static').mkdir(parents=True,exist_ok=True); shutil.move(str(st_/nm_),str(base/'supp/source_data/archive/static'/nm_))
        # canonical package inputs (crosswalk = DepMap 26Q3-merged version, 141 rows; guarded above)
        for f_ in [PIPE/'inputs/units_215.csv',XWALK/'crosswalk_hcmi_ach_sidm.csv',XWALK/'summary_crosswalk.txt',
                   STATIC/'tableS_TMP_vs_DNA_surrogate_consistency.csv',STATIC/'TMP_blocker_note.txt',STATIC/'hcmi_S12_gyn_own_tumour_rank.csv',
                   cfg_path('cptac_rule_evaluation'),cfg_path('paired_rna'),cfg_path('paired_methylation'),cfg_path('paired_mutation'),
                   cfg_path('uterine_units_classes'),MOL/'mlh1_threshold.json',MOL/'hcmi_sample_molecular_table.csv',MOL/'hcmi_case_labels.csv',MOL/'mlh1_promoter_probes.csv',
                   cfg_path('tcga_gyn_full_manifest'),cfg_path('tcga_nongyn_dir')/'nongyn_full_manifest.csv',cfg_path('tcga_nongyn_dir')/'release_check.txt',
                   cfg_path('tcga_nongyn_dir')/'gyn_manifest_release_check.csv',cfg_path('gdc_status_json')]:
            if f_.exists(): copy_clean(f_,st_/f_.name)
        if a.tag=='26Q3' and a.write_expected:   # expected outputs for the clean-room check (tests/compare_table2.py); maintainers only
            ex_=PIPE/'tests/expected'; ex_.mkdir(parents=True,exist_ok=True)
            for f_ in [RES/'table2_strict_loose.csv',RES/'units_215_calls.csv',RES/'key_numbers.json',TAB/'Table2_lineage_fidelity.xlsx']: shutil.copy2(f_,ex_/f_.name)
        cd=base/'supp/code'; cd.mkdir(parents=True,exist_ok=True)
        for f_ in list(PIPE.glob('*.py'))+[PIPE/'config_paths.example.json',PIPE/'config_paths.cleanroom.json',PIPE/'environment.yml',PIPE/'README_pipeline.md',PIPE/'data_manifest.csv',PIPE/'reproduce_table2.sh',PIPE/'reproduce_all.sh']:
            if f_.exists(): shutil.copy2(f_,cd/f_.name)
        old_=cd/'inputs/intermediate/uterine_units_plus_nextgen.csv'   # legacy file: now only tests/legacy_reference/LEGACY_uterine_units_plus_nextgen.csv (round 5c)
        if old_.exists() and old_.read_bytes()==(PIPE/'tests/legacy_reference/LEGACY_uterine_units_plus_nextgen.csv').read_bytes(): old_.unlink()
        for sub in ['inputs','inputs/crosswalk','inputs/static','inputs/annotation','inputs/manifests','inputs/intermediate','tests','tests/expected','tests/legacy_reference']:
            (cd/sub).mkdir(parents=True,exist_ok=True)
            for f_ in (PIPE/sub).glob('*.*'): shutil.copy2(f_,cd/sub/f_.name)
        r24=R24
        if not a.tag.startswith('24Q4') and r24.exists():   # DepMap 24Q4 sensitivity package (source data + its tables/figures)
            d24=base/'supp/source_data/24Q4'; d24.mkdir(parents=True,exist_ok=True)
            internal_audit(d24,base,'24Q4')
            archive_stale(d24,{f_.name for f_ in r24.glob('*.*') if f_.name!='number-audit.csv'},base/'supp/source_data/archive/24Q4')
            for f_ in r24.glob('*.*'):
                if f_.name!='number-audit.csv': copy_clean(f_,d24/f_.name)
            a24=base/'supp/source_data/archive/24Q4'; a24.mkdir(parents=True,exist_ok=True)   # 24Q4 identity re-run: archive only (Fabel, round 5)
            for f_ in (r24/'identity').glob('*.*'): copy_clean(f_,a24/('identity_'+f_.name))
            for sub in ['tables','figures']:
                dd=base/'supp/sensitivity_DepMap24Q4'/sub; dd.mkdir(parents=True,exist_ok=True)
                keep24={p.name for p in (r24/sub).glob('*') if not (p.name.startswith('Fig1_') or p.name.startswith(TSTEM('identity')))}   # identity panels/table of the 24Q4 run: archive only
                archive_stale(dd,keep24,base/'supp/source_data/archive/24Q4'/sub)
                for f_ in (r24/sub).glob('*'):
                    if f_.name in keep24: shutil.copy2(f_,dd/f_.name)
                    else: (base/'supp/source_data/archive/24Q4'/sub).mkdir(parents=True,exist_ok=True); shutil.copy2(f_,base/'supp/source_data/archive/24Q4'/sub/f_.name)
        shutil.copy2(RES/'number-audit.csv',base/'number-audit.csv')
    print('published to',MS,'and',DELIV)
