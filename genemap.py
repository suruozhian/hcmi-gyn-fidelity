"""Gene-identifier harmonisation to NCBI Entrez Gene IDs (primary gene key from round 4, 2026-10-08).

Official mapping: HGNC complete set (hgnc_complete_set.txt; source and version recorded in
inputs/annotation/hgnc_provenance.json). Rules (also in Methods 2.2):
  Ensembl gene IDs (GENCODE v36 protein-coding, version suffix removed):
    E1  Ensembl ID -> the HGNC record carrying that ensembl_gene_id -> entrez_id.
        An Ensembl ID carried by >1 HGNC record with different Entrez IDs (one-to-many) is excluded.
    E2  Ensembl ID not in HGNC (or without Entrez): GENCODE v36 gene name = HGNC approved symbol (exact, unique)
        -> entrez_id, only if that Entrez ID is not already assigned by E1.
    E3  otherwise excluded.
  Gene symbols (resources that provide symbols only: Sanger organoids, E-MTAB-10664):
    S1  approved symbol -> entrez_id;  S2  else previous symbol that belongs to exactly one HGNC record -> entrez_id;
    S3  otherwise excluded (alias symbols are not used: ambiguous).
  DepMap / NextGen files: Entrez ID taken from the column header "SYMBOL (Entrez)".
  Many-to-one (several source identifiers -> one Entrez ID, within one resource): maximum per sample
  (same rule as the previous duplicate-symbol rule).
"""
import json, re
from pathlib import Path
import numpy as np, pandas as pd
from paths import path as cfg_path

_H=None
def hgnc():
    global _H
    if _H is None:
        H=pd.read_csv(cfg_path('hgnc_complete_set'),sep='\t',dtype=str,low_memory=False)
        _H=H[H.status=='Approved'].copy()
    return _H

def ensembl_map(gencode_pc):
    """gencode_pc: DataFrame indexed by Ensembl ID (no version) with column gene_name (GENCODE v36 protein-coding).
    Returns DataFrame ens -> entrez (str), rule."""
    H=hgnc().dropna(subset=['entrez_id'])
    he=H.dropna(subset=['ensembl_gene_id'])
    amb=he.groupby('ensembl_gene_id').entrez_id.nunique(); amb=set(amb[amb>1].index)
    e1=he[~he.ensembl_gene_id.isin(amb)].drop_duplicates('ensembl_gene_id').set_index('ensembl_gene_id').entrez_id
    out=pd.DataFrame(index=gencode_pc.index); out['gene_name']=gencode_pc.gene_name
    out['entrez']=out.index.map(e1); out['rule']=np.where(out.entrez.notna(),'E1','')
    out.loc[out.index.isin(amb),'rule']='excluded: one-to-many'
    taken=set(out.entrez.dropna())
    sym=H.groupby('symbol').entrez_id.agg(lambda s: s.iloc[0] if s.nunique()==1 else None)
    need=out.entrez.isna()&(out.rule=='')
    cand=out.loc[need,'gene_name'].map(sym)
    ok=cand.notna()&~cand.isin(taken)
    out.loc[cand.index[ok],'entrez']=cand[ok]; out.loc[cand.index[ok],'rule']='E2'
    out.loc[out.entrez.isna()&(out.rule==''),'rule']='excluded: not in HGNC with Entrez'
    return out

def symbol_map(symbols):
    H=hgnc().dropna(subset=['entrez_id'])
    appr=H.set_index('symbol').entrez_id
    prev={}
    for e,ps in zip(H.entrez_id,H.prev_symbol.fillna('')):
        for s in [x.strip() for x in ps.split('|') if x.strip()]: prev.setdefault(s,set()).add(e)
    rows=[]
    for s in pd.Index(symbols).unique():
        if s in appr.index: rows.append((s,appr[s],'S1'))
        elif s in prev and len(prev[s])==1: rows.append((s,next(iter(prev[s])),'S2'))
        else: rows.append((s,None,'excluded'))
    return pd.DataFrame(rows,columns=['symbol','entrez','rule']).set_index('symbol')

def to_entrez(M,mapping):
    """M: genes x samples indexed by source identifier; mapping: Series source id -> entrez (NaN = drop). Max over many-to-one."""
    m=mapping.dropna(); M=M[M.index.isin(m.index)].copy(); M.index=m.loc[M.index].astype(str).values
    return M.groupby(level=0).max()

def depmap_entrez(cols):
    """'SYMBOL (Entrez)' -> Entrez string (None if no Entrez)."""
    out=[]
    for c in cols:
        mm=re.search(r'\((\d+)\)$',str(c).strip('"')); out.append(mm.group(1) if mm else None)
    return out

_GENCODE=None
def gencode_pc():
    global _GENCODE
    if _GENCODE is None:
        ga=pd.read_csv(cfg_path('gene_annotation'),index_col=0); ga['ens']=ga.index.str.split('.').str[0]
        _GENCODE=ga[ga.gene_type=='protein_coding'].drop_duplicates('ens').set_index('ens')
    return _GENCODE

_EM=None
def ens2entrez():
    global _EM
    if _EM is None: _EM=ensembl_map(gencode_pc())
    return _EM

def entrez2symbol():
    """display symbol for an Entrez ID (HGNC approved symbol)."""
    H=hgnc().dropna(subset=['entrez_id']); return H.drop_duplicates('entrez_id').set_index('entrez_id').symbol

def gencode_symbols_to_entrez(symbols):
    """symbols that came from GENCODE v36 gene_name (e.g. the purity-filtered list): all protein-coding Ensembl IDs with that
    name -> Entrez via E1/E2."""
    g=gencode_pc(); em=ens2entrez()
    ens=g.index[g.gene_name.isin(set(symbols))]
    return sorted(set(em.loc[ens,'entrez'].dropna()))
