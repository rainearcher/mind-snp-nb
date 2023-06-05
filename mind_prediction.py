import subprocess
import os
def examine_snp_effect(path_to_fasta, path_to_res, snp):
    command = f"""python3 PTMSNP.py --pretrain_name saved_model/MIND_fifteenfold --data_path {path_to_fasta} --res_path {path_to_res} --snp {snp} --n_fold 15"""
    print("running command: " + command)
    subprocess.run([command], shell=True)

RES_PATH = "mind_files/results"
def predict_all_ptms(df):
    for filename in os.listdir("mind_files/fastas"):
        fasta_path = f"mind_files/fastas/{filename}"
        fasta, ext = filename.split(".")
        snps = df.loc[df['accession_id'] == fasta]['mind_snp'].to_list()
        for snp in snps:
            examine_snp_effect(fasta_path, RES_PATH, snp)