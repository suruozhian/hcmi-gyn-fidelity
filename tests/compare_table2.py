#!/usr/bin/env python3
"""Clean-room check: compare a fresh run with the published Table 2, cell for cell.

  python3 tests/compare_table2.py --results results/26Q3 [--expected tests/expected] [--out reproduce_table2_check.json]

Compared (all cells as text, exact):
  1. Table 2 sheet rebuilt from the fresh results (table2.build_table2) vs tests/expected/Table2_lineage_fidelity.xlsx (sheet Table2)
  2. results/<tag>/table2_strict_loose.csv (all strata, Method A and B columns) vs tests/expected/table2_strict_loose.csv
  3. results/<tag>/units_215_calls.csv (215 units x all columns) vs tests/expected/units_215_calls.csv
Exit code 0 only if every compared cell is identical."""
import argparse, json, sys
from pathlib import Path
import pandas as pd
PIPE=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(PIPE))
from table2 import build_table2
ap=argparse.ArgumentParser(); ap.add_argument('--results',required=True); ap.add_argument('--expected',default=str(PIPE/'tests/expected'))
ap.add_argument('--out',default=None); a=ap.parse_args()
R=Path(a.results); E=Path(a.expected)
def cmp(new,old,name):
    new=new.astype(str).reset_index(drop=True); old=old.astype(str).reset_index(drop=True)
    shape_ok=new.shape==old.shape and list(new.columns)==list(old.columns)
    diffs=[]
    if shape_ok:
        for c in old.columns:
            for i in old.index[(old[c]!=new[c])]: diffs.append(dict(row=int(i),column=c,expected=old.at[i,c],got=new.at[i,c]))
    return dict(table=name,shape_expected=list(old.shape),shape_got=list(new.shape),columns_identical=shape_ok,cells_compared=int(old.size) if shape_ok else 0,
                cells_different=len(diffs) if shape_ok else None,examples=diffs[:10])
T2=pd.read_csv(R/'table2_strict_loose.csv').set_index('stratum'); K=json.load(open(R/'key_numbers.json'))
t2_new=build_table2(T2,K)
t2_old=pd.read_excel(E/'Table2_lineage_fidelity.xlsx',sheet_name='Table2',dtype=str).fillna('')
t2_new=t2_new.astype(str).replace('nan','')
res=[cmp(t2_new,t2_old,'Table 2 (sheet Table2)'),
     cmp(pd.read_csv(R/'table2_strict_loose.csv',dtype=str,keep_default_na=False),pd.read_csv(E/'table2_strict_loose.csv',dtype=str,keep_default_na=False),'table2_strict_loose.csv'),
     cmp(pd.read_csv(R/'units_215_calls.csv',dtype=str,keep_default_na=False),pd.read_csv(E/'units_215_calls.csv',dtype=str,keep_default_na=False),'units_215_calls.csv')]
ok=all(r['columns_identical'] and r['cells_different']==0 for r in res)
out=dict(identical=ok,checks=res)
print(json.dumps(out,indent=1))
if a.out: json.dump(out,open(a.out,'w'),indent=1)
sys.exit(0 if ok else 1)
