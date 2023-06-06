GWAS_XLSX_PATH = "gwas_catalog_cardiovascular_disease.xlsx"
GWAS_TXT_PATH = "temp_files/gwas_catalog.txt"
RISK_ALLELE_PATH = "temp_files/riskAlleles.json"
from apis import EnsemblRestClient
from save import save
import json
import os
import pandas as pd

def clean_gwas_catalog(path_to_gwas):
    df = pd.read_excel(path_to_gwas)
    df['STRONGEST SNP-RISK ALLELE'].to_csv(GWAS_TXT_PATH, index=False, header=False)

def chrToRS(region):
    client = EnsemblRestClient()
    overlap = client.get_overlap(region)
    return overlap[0]["id"] if len(overlap) else None

def get_risk_alleles(path_to_associations=None):
    risk_alleles = {}
    f = open(path_to_associations, "r")
    for line in f:
        a, nt, *_ = line.split('-')
        if ("chr" in a):
            a = chrToRS(a)
        if (a is None):
            continue
        risk_alleles[a] = nt
    f.close()
    return risk_alleles

def save_risk_alleles(risk_alleles, path):
    save(risk_alleles, path)

def create_rsid_df():
    if not os.path.exists(GWAS_TXT_PATH):
        clean_gwas_catalog(GWAS_XLSX_PATH)
    
    if os.path.exists(RISK_ALLELE_PATH):
        f = open(RISK_ALLELE_PATH, "r")
        risk_alleles = json.load(f)
        f.close()
    else:
        risk_alleles = get_risk_alleles(GWAS_TXT_PATH)
        save_risk_alleles(risk_alleles, RISK_ALLELE_PATH)
    
    df = pd.DataFrame()
    df['rs_id'] = [k.strip() for k in risk_alleles.keys()]
    df['nucleotide'] = [v.strip() for v in risk_alleles.values()]
    return df

