"""Extracts from public web APIs (used by fetch_public_data.py; data_manifest.csv rows with retrieval == "api").

Each function writes one CSV deterministically; fetch_public_data.py checks its md5 against data_manifest.csv.
- cBioPortal (https://www.cbioportal.org/api), study ucec_cptac_2020 (CPTAC UCEC; Dou et al., Cell 2020): sample genomic
  subtype (clinical attribute GENOMICS_SUBTYPE), mutations (profile ucec_cptac_2020_mutations) and MLH1 methylation
  (profile ucec_cptac_2020_methylation_epic, Entrez 4292). Open access.
- NCI GDC (https://api.gdc.cancer.gov/cases): primary site and primary diagnoses of the HCMI gynaecological cases
  (project HCMI-CMDC; primary site Ovary, Corpus uteri, Uterus NOS, Cervix uteri). Open access."""
import json, urllib.request
import pandas as pd
UA={'User-Agent':'gyn-model-fidelity-fetch/1.0','Accept':'application/json','Content-Type':'application/json'}
CBP='https://www.cbioportal.org/api'
def _call(url,body=None):
    req=urllib.request.Request(url,data=None if body is None else json.dumps(body).encode(),headers=UA)
    with urllib.request.urlopen(req,timeout=600) as r: return json.load(r)
def cptac_ucec_genomic_subtype(dest):
    d=_call(CBP+'/studies/ucec_cptac_2020/clinical-data?clinicalDataType=SAMPLE&projection=SUMMARY')
    d=pd.DataFrame(d); d=d[d.clinicalAttributeId=='GENOMICS_SUBTYPE']
    d[['sampleId','patientId','value']].rename(columns={'value':'GENOMICS_SUBTYPE'}).sort_values('sampleId').to_csv(dest,index=False)
def cptac_ucec_mutations(dest):
    mu=_call(CBP+'/molecular-profiles/ucec_cptac_2020_mutations/mutations/fetch?projection=DETAILED',{'sampleListId':'ucec_cptac_2020_all'})
    pd.DataFrame([dict(sampleId=x['sampleId'],gene=x['gene']['hugoGeneSymbol'],proteinChange=x.get('proteinChange'),
                       mutationType=x.get('mutationType'),variantType=x.get('variantType')) for x in mu]).to_csv(dest,index=False)
def cptac_ucec_mlh1_methylation(dest):
    m=_call(CBP+'/molecular-profiles/ucec_cptac_2020_methylation_epic/molecular-data?sampleListId=ucec_cptac_2020_all&entrezGeneId=4292&projection=SUMMARY')
    pd.DataFrame([dict(sampleId=x['sampleId'],MLH1_meth=x['value']) for x in m]).to_csv(dest,index=False)
def hcmi_case_diagnoses(dest):
    body={'filters':{'op':'and','content':[{'op':'in','content':{'field':'project.project_id','value':['HCMI-CMDC']}},
          {'op':'in','content':{'field':'primary_site','value':['Ovary','Corpus uteri','Uterus, NOS','Cervix uteri']}}]},
          'fields':'submitter_id,primary_site,diagnoses.primary_diagnosis','size':10000,'format':'json'}
    h=_call('https://api.gdc.cancer.gov/cases',body)['data']['hits']
    pd.DataFrame([dict(case=x['submitter_id'],primary_site=x['primary_site'],
                       primary_diagnosis=';'.join(sorted({d['primary_diagnosis'] for d in x.get('diagnoses',[])}))) for x in h]
                ).sort_values('case').to_csv(dest,index=False)
EXTRACT={f.__name__:f for f in [cptac_ucec_genomic_subtype,cptac_ucec_mutations,cptac_ucec_mlh1_methylation,hcmi_case_diagnoses]}
