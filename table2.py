"""Table 2 (main text): built from results/<tag>/table2_strict_loose.csv and key_numbers.json.
Shared by make_outputs.py (writes tables/Table2_lineage_fidelity.xlsx) and tests/compare_table2.py (clean-room check)."""
import pandas as pd

def build_table2(T2,K):
    """T2: table2_strict_loose.csv indexed by stratum; K: key_numbers.json (dict). Returns the 'Table2' sheet."""
    def t2row(lbl,key):
        r=T2.loc[key]; return {'Stratum':lbl,'Units':int(r.n),'Strict Method A':f'{int(r.strict_A)}/{int(r.n)} ({r.strict_A_pct:.1f}%)','Loose Method A':f'{int(r.loose_A)}/{int(r.n)} ({r.loose_A_pct:.1f}%)',
            'Organoid strict A':f'{int(r.organoid_strict_A)}/{int(r.organoid_n)}' if r.organoid_n else '—','2D strict A':f'{int(r.twoD_strict_A)}/{int(r.twoD_n)}' if r.twoD_n else '—'}
    t2=pd.DataFrame([t2row('HGSC/serous (incl. serous NOS)','HGSC/serous (incl. serous NOS)'),t2row('  High-grade serous (HGSC)','  OV_HGSC'),t2row('  HCMI/Sanger serous','  HCMI/Sanger OV_serous'),
        t2row('  DepMap serous NOS','  DepMap serous NOS'),t2row('Endometrial carcinoma (non-CS)','Endometrial carcinoma (non-CS)'),t2row('  EC organoids','  EC organoids'),t2row('  EC 2D lines','  EC 2D lines'),
        t2row('Uterine carcinosarcoma','Uterine carcinosarcoma')])
    cp=K['cptac']
    t2=pd.concat([t2,pd.DataFrame([{'Stratum':'CPTAC OV tumours (transfer)','Units':cp['OV_n'],'Strict Method A':f"{cp['OV_nearestOV']}/{cp['OV_n']} nearest OV"},
                                   {'Stratum':'CPTAC UCEC tumours (transfer)','Units':cp['UCEC_n'],'Strict Method A':f"{cp['UCEC_nearest_UCEC']}/{cp['UCEC_n']} nearest UCEC",'Loose Method A':f"{cp['UCEC_nearest_UCEC_or_UCS']}/{cp['UCEC_n']} UCEC or UCS"}])]).fillna('—')
    return t2
