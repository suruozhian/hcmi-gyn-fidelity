#!/usr/bin/env python3
"""Supplementary Fig. S7: HCMI uterine models (M) and matched tumours (T), DNA/methylation features underlying the
surrogate classes and the resulting class (not IHC). Input: <molecular_dir>/hcmi_sample_molecular_table.csv
(04_derive_hcmi_molecular.py). Output: <molecular_dir>/hcmi_endometrial_surrogate.png (300 dpi) and .pdf (the PNG
embedded by Pillow, as in the published file; its internal title is the earlier file stem)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns
from matplotlib.patches import Patch
from PIL import Image
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path
OUT=cfg_path('molecular_dir'); STEM='hcmi_endometrial_surrogate'   # copied to figures/ under its supplementary number by make_outputs.py
plt.rcParams.update({'font.size':8,'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False})
S=pd.read_csv(OUT/'hcmi_sample_molecular_table.csv')
order=['OV_serous','OV_mullerian_mixed','OV_borderline','OV_endometrioid','OV_clear_cell','OV_mucinous','UT_endometrioid','UT_mixed','UT_serous','UT_undifferentiated','UT_carcinosarcoma']
u=S[(S.lineage=='uterus')&(S.EC_class!='no data')].copy()
u['h']=pd.Categorical(u.histology,order); u=u.sort_values(['h','case','role'])
feat=pd.DataFrame({'POLE EDM hotspot':u.POLE_hotspot.fillna('').ne('').astype(float),
 'MMR gene truncating':u.MMR_truncating.fillna('').ne('').astype(float),
 'MLH1 promoter β':u.MLH1_prom_beta,
 'TP53 mutation':u.TP53.fillna('').ne('').astype(float),
 'CN fraction altered':u.cn_frac_altered,
 'log10 nonsyn muts / 6':np.log10(u.n_nonsyn.clip(lower=1))/6})
lab=[f"{s} [{'M' if r=='model' else 'T'}] {h.replace('UT_','')}" for s,r,h in zip(u['sample'],u.role,u.histology)]
fig,(ax,ax2)=plt.subplots(1,2,figsize=(7.5,6.5),gridspec_kw={'width_ratios':[6,1.4]})
sns.heatmap(feat.values,ax=ax,cmap='Purples',vmin=0,vmax=1,yticklabels=lab,xticklabels=feat.columns,cbar_kws={'shrink':0.4,'label':'value (0–1 scale)'},mask=feat.isna().values,linewidths=0.3,linecolor='white')
ax.tick_params(axis='y',labelsize=6); ax.tick_params(axis='x',labelsize=6.5,rotation=40)
cmap={'POLEmut (surrogate)':'#1b9e77','MMRd (surrogate)':'#7570b3','p53abn (surrogate)':'#d95f02','NSMP (surrogate)':'#bbbbbb','incomplete (no MAF)':'white'}
for i,c in enumerate(u.EC_class): ax2.add_patch(plt.Rectangle((0,i),1,1,color=cmap[c]))
ax2.set_ylim(len(u),0); ax2.set_xlim(0,1); ax2.axis('off')
ax2.legend(handles=[Patch(color=v,label=k) for k,v in cmap.items() if k in set(u.EC_class)],fontsize=6,loc='upper left',bbox_to_anchor=(0,-0.01),frameon=False)
ax2.set_title('Surrogate class\n(not IHC)',fontsize=7)
fig.suptitle('Uterine carcinoma models (M) and tumours (T): molecular-class surrogates from open GDC data',fontsize=8.5)
plt.tight_layout(); plt.savefig(OUT/f'{STEM}.png',dpi=300,bbox_inches='tight'); plt.close()
Image.open(OUT/f'{STEM}.png').convert('RGB').save(OUT/f'{STEM}.pdf','PDF',resolution=300,title='FigS5_hcmi_endometrial_surrogate')
print('Supplementary Fig. S7:',len(u),'HCMI uterine samples')
