from apis import EnsemblRestClient
import json
import os
import pandas as pd
from save import save
#get variant effects for each SNP
VARIANT_INFO = ["transcript_id", "amino_acids", "gene_symbol", "gene_symbol_source", "protein_start", "protein_end", "variant_allele"]
RISK_ALLELE_PATH = "temp_files/riskAlleles.json"
def get_variants(df):
    client = EnsemblRestClient()
    variants = {}
    risk_alleles = dict(zip(df.rs_id, df.nucleotide))
    for snp in list(risk_alleles.keys()):
        effects = client.get_effects("human", snp)
        if not effects or effects[0]["most_severe_consequence"] != "missense_variant":
            continue
        transcripts = []
        for var in effects[0]["transcript_consequences"]:
            if ("missense_variant" in var["consequence_terms"]):
                t = { i : (var[i] if i in var else "?/?") for i in VARIANT_INFO }
                transcripts.append(t)
                #all_transcripts.add(var["transcript_id"])
        variants[snp] = transcripts
        #get all the information on the google doc
    return variants

def get_snps(variants={}, variants_path=None):
    if not variants:
        f = open(variants_path, "r")
        variants = json.load(f)
        f.close()
    snps = {}
    for rs_id_variants in variants.values():
        for variant in rs_id_variants:
            trans_id = variant["transcript_id"]

            #get info for MIND
            snp_wt, snp_var = variant["amino_acids"].split("/")
            snp_site = str(variant["protein_start"])
            snp = f"{snp_wt}_{snp_site}_{snp_var}"
            if (trans_id in snps):
                snps[trans_id].append(snp)
            else:
                snps[trans_id] = [snp]
    return snps

    
def map_snps_to_proteins(df):
    VARIANTS_PATH = "temp_files/variants.json"
    if os.path.exists(VARIANTS_PATH):
        f = open(VARIANTS_PATH, "r")
        variants = json.load(f)
        f.close()
    else:
        variants = get_variants(df)
        save(variants, VARIANTS_PATH)

    df_variants = pd.DataFrame(columns=['transcript_id', 'amino_acids', 'gene_symbol', 'gene_symbol_source', 'protein_start', 'protein_end', 'variant_allele', 'rs_id'])

    for rsid, var_list in variants.items():
        for var in var_list:
            var['rs_id'] = rsid.strip()
            df_variants.loc[len(df_variants)] = var

    df = df.merge(df_variants, on='rs_id', how='right')
    def add_snp(row):
        wt, mut = row['amino_acids'].split('/')
        site = row['protein_start']
        return f"{wt}_{site}_{mut}"
    df['mind_snp'] = df.apply(add_snp, axis=1)

    return df