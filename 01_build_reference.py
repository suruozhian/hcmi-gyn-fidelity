#!/usr/bin/env python3
"""Build the 11-project full TCGA reference (OV/UCEC/UCS full + 8 non-gyn full), log2(TPM+1),
parsing: STAR tpm_unstranded, GENCODE v36 protein-coding genes, log2(TPM+1).
Gene key: --ids entrez (default; Ensembl -> Entrez via HGNC, rules in genemap.py; max over many-to-one)
          --ids symbol (GENCODE v36 gene name; max over duplicate names; sensitivity analysis).
Optional --check-only: verify that every manifest file resolves locally, without parsing."""
import sys
import pandas as pd, numpy as np
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0,str(Path(__file__).resolve().parent))
from paths import path as cfg_path, has as cfg_has, locate_tcga
import genemap
OUT=cfg_path('reference_dir'); OUT.mkdir(exist_ok=True)
ga=pd.read_csv(cfg_path('gene_annotation'),index_col=0); ga['ens']=ga.index.str.split('.').str[0]
pc=ga[ga.gene_type=='protein_coding']; sym=pc.drop_duplicates('ens').set_index('ens').gene_name
GENES=pd.Index(sym.index)   # explicit gene order; every parsed file is reindexed to it (memory-light: arrays, not per-file indexes)
def parse_one(a):
    path,case,project=a
    t=pd.read_csv(path,sep='\t',comment='#',usecols=['gene_id','tpm_unstranded'])
    t=t[t.gene_id.str.startswith('ENSG')].copy(); t['ens']=t.gene_id.str.split('.').str[0]
    t=t[t.ens.isin(sym.index)].drop_duplicates('ens').set_index('ens')
    v=t.tpm_unstranded.reindex(GENES).to_numpy(dtype='float32')
    assert not np.isnan(v).any(), f'{path}: missing protein-coding genes'
    return project,case,v
if __name__=='__main__':
    G=pd.read_csv(cfg_path('tcga_gyn_full_manifest'))
    N=pd.read_csv(cfg_path('tcga_nongyn_dir')/'nongyn_full_manifest.csv')
    assert N.md5_ok.all(), 'non-gyn md5 not all ok'
    man=pd.concat([G[['file_id','project','case','sample','local_path']],N[['file_id','project','case','sample','local_path']]],ignore_index=True)
    assert man.case.is_unique
    args=[(str(locate_tcga(r.local_path)),r.case,r.project.replace('TCGA-','')) for r in man.itertuples()]
    missing=[x[0] for x in args if not Path(x[0]).exists()]
    print('manifest files',len(args),'missing locally',len(missing))
    if '--check-only' in sys.argv: sys.exit(0 if not missing else 1)
    assert not missing, missing[:5]
    with ProcessPoolExecutor(8) as ex: res=list(ex.map(parse_one,args,chunksize=20))
    T=pd.DataFrame(np.column_stack([v for _,_,v in res]),index=GENES,columns=['TCGA|'+c for _,c,_ in res]); del res
    IDS='symbol' if '--ids' in sys.argv and sys.argv[sys.argv.index('--ids')+1]=='symbol' else 'entrez'
    T=np.log2(T+1)
    if IDS=='entrez':
        em=genemap.ens2entrez(); T=genemap.to_entrez(T,em.entrez).astype('float32')
        em.to_csv(OUT/'ensembl_to_entrez_map.csv')
    else:
        T.index=sym.loc[T.index].values; T=T.groupby(level=0).max().astype('float32')
    T.to_parquet(OUT/f'tcga11_full_{IDS}_log2tpm.parquet')
    meta=pd.DataFrame([dict(sample='TCGA|'+r[1],project=r[2],case=r[1]) for r in args])
    meta=meta.merge(man.assign(case=man.case)[['case','sample','file_id']].rename(columns={'sample':'tcga_sample'}),on='case')
    meta.to_csv(OUT/'tcga11_full_samples.csv',index=False)
    # optional check: gyn columns identical to the earlier gynaecological-only parse
    if IDS=='symbol' and cfg_has('tcga_previous_gyn_parse') and cfg_path('tcga_previous_gyn_parse').exists():
        old=pd.read_parquet(cfg_path('tcga_previous_gyn_parse'))
        cols=[c for c in old.columns if c in T.columns]; g=old.index.intersection(T.index)
        d=float(np.nanmax(np.abs(old.loc[g,cols].values-T.loc[g,cols].values)))
    else: cols=[]; d=float('nan')
    s=meta.groupby('project').size().to_dict()
    open(OUT/f'reference_build_summary_{IDS}.txt','w').write(f'gene key {IDS}\nmatrix {T.shape}\nper project {s}\nmax abs diff vs earlier gynaecological parse on {len(cols)} shared columns: {d:.2e}\n')
    print(open(OUT/f'reference_build_summary_{IDS}.txt').read())
