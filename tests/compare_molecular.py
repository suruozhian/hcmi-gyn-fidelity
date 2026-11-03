#!/usr/bin/env python3
"""Compare the molecular intermediates re-derived from public data (04/05/06 scripts, <molecular_dir>) with the copies
shipped in inputs/intermediate and inputs/static (byte level, then cell level). Writes a JSON report; exit 1 on any difference.
The uterine unit classes are compared with tests/legacy_reference/LEGACY_uterine_units_plus_nextgen.csv, the file formerly shipped
(written by an exploratory script with provisional per-patient E-MTAB-10664 units 'Berg-pt<NN>'; see README_legacy.md): rows are
matched on dataset + unit (E-MTAB-10664 via the patient number). The surrogate class must agree for every unit; the one known
histology difference of the legacy file (Berg-pt04 UT_endometrioid; final unit OEC-04-CC UT_clear_cell) is reported and expected."""
import argparse, hashlib, json, re, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from paths import path as cfg_path
PIPE=Path(__file__).resolve().parents[1]; IN=PIPE/'inputs/intermediate'; ST=PIPE/'inputs/static'
ap=argparse.ArgumentParser(); ap.add_argument('--out',default='molecular_rederivation_check.json'); a=ap.parse_args()
D=cfg_path('molecular_dir'); md5=lambda p: hashlib.md5(Path(p).read_bytes()).hexdigest()
def cells(x,y):
    if list(x.columns)!=list(y.columns) or x.shape!=y.shape: return dict(shape_or_columns_differ=True,derived=[list(x.columns),x.shape],shipped=[list(y.columns),y.shape])
    d=[]
    for c in x.columns:
        for i,(u,v) in enumerate(zip(x[c],y[c])):
            if not ((pd.isna(u) and pd.isna(v)) or u==v): d.append(dict(row=i,col=c,derived=str(u),shipped=str(v)))
    return dict(n_cells=int(x.size),n_differ=len(d),differences=d[:50])
R={}
for f in ['rna_paired_purityfilt.csv','paired_mutation_concordance.csv','paired_methylation_concordance.csv','cptac_rule_evaluation.csv']:
    x,y=D/f,IN/f; same=md5(x)==md5(y)
    R[f]=dict(byte_identical=same,md5_derived=md5(x),md5_shipped=md5(y),**({} if same else cells(pd.read_csv(x),pd.read_csv(y))))
# uterine classes vs the shipped legacy table
x=pd.read_csv(D/'uterine_units_classes.csv'); y=pd.read_csv(PIPE/'tests/legacy_reference/LEGACY_uterine_units_plus_nextgen.csv')
def legacy_key(r): return 'Berg-pt'+re.search(r'-(\d+)-',r.model_unit+'-').group(1) if r.dataset.startswith('Berg') else r.model_unit
x['legacy_unit']=[legacy_key(r) for r in x.itertuples()]
m=x.merge(y,left_on=['dataset','legacy_unit'],right_on=['dataset','model_unit'],how='outer',suffixes=('','_shipped'),indicator=True)
KNOWN=[dict(dataset='Berg2021_EMTAB10664',legacy_unit='Berg-pt04',column='histology',derived='UT_clear_cell',shipped='UT_endometrioid')]
KNOWN_SET={tuple(k.values()) for k in KNOWN}
dif=[dict(dataset=r.dataset,model_unit=r.model_unit,legacy_unit=r.legacy_unit,column=c,derived=getattr(r,c),shipped=getattr(r,c+'_shipped'))
     for r in m[m._merge=='both'].itertuples() for c in ['histology','cls'] if getattr(r,c)!=getattr(r,c+'_shipped')]
R['uterine_units_classes.csv']=dict(byte_identical=False,rows_derived=len(x),rows_shipped=len(y),only_derived=int((m._merge=='left_only').sum()),
    only_shipped=int((m._merge=='right_only').sum()),cls_differ=sum(d['column']=='cls' for d in dif),histology_differ=sum(d['column']=='histology' for d in dif),
    differences=dif,class_counts_by_dataset_identical=bool(pd.crosstab(x.dataset,x.cls).equals(pd.crosstab(y.dataset,y.cls))),
    known_legacy_differences=KNOWN,unexpected_differences=[d for d in dif if (d['dataset'],d['legacy_unit'],d['column'],d['derived'],d['shipped']) not in KNOWN_SET],
    note='compared with the legacy file on dataset + unit; row order and E-MTAB-10664 unit names follow the final unit definitions (inputs/units_215.csv)')
for ext in ['png','pdf']:
    x,y=D/f'hcmi_endometrial_surrogate.{ext}',ST/f'FigS5_hcmi_endometrial_surrogate.{ext}'
    bx,by=x.read_bytes(),y.read_bytes()
    if ext=='pdf': bx,by=[re.sub(rb'/(CreationDate|ModDate) \(D:[^)]*\)',b'',b) for b in (bx,by)]
    R[f'FigS7.{ext}']=dict(byte_identical=bx==by,**({'compared':'after removing /CreationDate and /ModDate'} if ext=='pdf' else {}))
ok=all(v['byte_identical'] for k,v in R.items() if k!='uterine_units_classes.csv')
u=R['uterine_units_classes.csv']; uok=u['only_derived']==0 and u['only_shipped']==0 and u['cls_differ']==0 and not u['unexpected_differences'] and u['class_counts_by_dataset_identical']
R['summary']=dict(all_match=bool(ok and uok),note='uterine classes: class identical for every unit; the known legacy histology difference is listed, not counted as a failure')
json.dump(R,open(a.out,'w'),indent=1,default=str)
for k,v in R.items(): print(k,{kk:vv for kk,vv in v.items() if kk not in ('differences',)})
sys.exit(0 if R['summary']['all_match'] else 1)
