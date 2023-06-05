from apis import UniProtRestClient
from functools import reduce
import json
import os
import pandas as pd

FASTA_PATH = "mind_files/fastas"

#filter out IDs that are not reviewed
def mapTransIdToAccessId(ids=None):
    if (ids is None):
        f = open("temp_files/ids.json", "r")
        ids = json.load(f)
        f.close()
    filtered_ids = {}
    for transcript in ids["results"]:
        filtered_ids[transcript["to"]["primaryAccession"]] = transcript["from"]
        
    return filtered_ids

def saveFastas(fastas, fasta_dir):
    for t in fastas:
        start = 0
        accession_id = ""
        for i in range(len(t)):
            if(t[i] == "|"):
                if start == 0:
                    start = i + 1
                else:
                    accession_id = t[start:i]
        with open(os.path.join(FASTA_PATH, f"{accession_id}.fasta"), "w") as f:
            f.write(t)


IDS_PATH = "temp_files/ids.json"
def download_fastas_save_ids(df):
    client = UniProtRestClient()
    trans_ids = df['transcript_id'].drop_duplicates().to_list()
    cs_ids = reduce(lambda val, el: val + ',' + el, trans_ids)
    
    if os.path.exists(IDS_PATH):
        f = open(IDS_PATH, 'r')
        ids = json.load(f)
        f.close()
    else:
        fastas, ids = client.id_map(cs_ids)
        saveFastas(fastas, FASTA_PATH)
        ids = mapTransIdToAccessId(ids)

    df_ids = pd.DataFrame(columns=['transcript_id', 'accession_id'])
    df_ids['transcript_id'] = ids.values()
    df_ids['accession_id'] = ids.keys()

    df = df.merge(df_ids, on='transcript_id')
    df = df.drop_duplicates()
    return df