import os
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

RES_PATH = "mind_files/results"
SAVE_PATH = "mind_files"

def create_results_df(data_df):
    files = os.listdir(RES_PATH)
    mutant_df = pd.DataFrame(columns=["ProteinID", "Mutation", "PTM", "Mutant_Prob"])
    orig_df = pd.DataFrame(columns = ["ProteinID", "Orig_Prob"])
    for filename in files:
        f = open(os.path.join(RES_PATH, filename), "r")
        dat = json.load(f)
        f.close()
        if "_" in filename:
            name = filename.split("_")
            df = pd.DataFrame.from_dict(dat, orient="index", columns=["Mutant_Prob"])
            df["PTM"] = df.index
            df.reset_index(level=0, inplace=True, drop=True)
            df["ProteinID"] = [ name[0] for i in df.index ]
            df["Mutation"] = [ name[1][:-5] for i in df.index ]
            mutant_df = pd.concat([mutant_df, df])
        else:
            name = filename.split(".")
            df = pd.DataFrame.from_dict(dat, orient="index", columns=["Orig_Prob"])
            df["PTM"] = df.index
            df.reset_index(level=0, inplace=True, drop=True)
            df["ProteinID"] = [ name[0] for i in df.index ]
            df = df[["ProteinID", "PTM", "Orig_Prob"]]
            orig_df = pd.concat([orig_df, df])
    all_df = pd.merge(mutant_df, orig_df, on=["ProteinID", "PTM"])
    all_df = all_df.astype({ "Orig_Prob":"float", "Mutant_Prob":"float" })
    df = all_df

    def rewrite_ptm(row):
        mutation = row['Mutation']
        wt = mutation[0]
        mut = mutation[-1]
        site = mutation[1:-1]
        return f"{wt}_{site}_{mut}"
    df['mind_snp'] = df.apply(rewrite_ptm, axis=1)

    return df.merge(data_df, left_on=['ProteinID', 'mind_snp'], right_on=['accession_id', 'mind_snp'])


def prefilter_results(df):
    def select(row):
        mut = row['Mutant_Prob']
        orig = row['Orig_Prob']
        diff = np.abs(mut - orig)
        if 0.5 < mut and orig < 0.5 and diff >= 0.2:
            return True
        if mut < 0.5 and 0.5 < orig and diff >= 0.2:
            return True
        return False
    filtered_df = df
    filtered_df['Diff'] = filtered_df.apply(lambda row: np.abs(row['Mutant_Prob'] - row['Orig_Prob']), axis=1)
    filtered_df = filtered_df.loc[ lambda row: ((0.5 < row['Mutant_Prob']) & (row['Orig_Prob'] < 0.5) & (0.2 <= row['Diff'])) |
                          ((0.5 > row['Mutant_Prob']) & (row['Orig_Prob'] > 0.5) & (0.2 <= row['Diff'])) ]
    filtered_df = filtered_df.sort_values(["Orig_Prob", "Mutant_Prob"], ascending=[False, True])
    filtered_df = filtered_df.reset_index(level=0, drop=True)
    return filtered_df

def filter_impacted_ptm_snps(df):
    df = df.loc[lambda row: ((row['Orig_Prob'] >= 0.6) & (row['Mutant_Prob'] <= 0.4))
                | ((row['Mutant_Prob'] >= 0.6) & (row['Orig_Prob'] <= 0.4)) ]
    return df.copy()

def save_results(results_df, fname="all_results.tsv"):
    results_df.to_csv(os.path.join(SAVE_PATH, fname), sep="\t")


def plot_pairs(df):
    df = df.copy()
    df = df.loc[ ( (0.5 < df["Mutant_Prob"]) & (df["Orig_Prob"] < 0.5) ) | ( (df["Mutant_Prob"] < 0.5) & (0.5 < df["Orig_Prob"]) ) ]
    df.sort_values(["Orig_Prob", "Mutant_Prob"], ascending=[False, True], inplace=True)
    df = df[["ProteinID", "Mutation", "PTM", "Orig_Prob", "Mutant_Prob"]]
    df.reset_index(level=0, drop=True, inplace=True)

    df['index'] = df.index
    df.rename(columns={'Orig_Prob': 'Original', 'Mutant_Prob': 'Mutant'}, inplace=True)
    cols = ['index', 'Original', 'Mutant']
    df.drop([c for c in df.columns if c not in cols], axis=1, inplace=True)
    dfm = df.melt('index', var_name='Protein', value_name='vals')

    plt.figure()
    sns.set_theme(style='dark')
    sns.color_palette("muted")
    ax = sns.scatterplot(dfm, x='index', y='vals', hue='Protein', linewidth=0, alpha=0.65)
    #plt.scatter(X, orig_prob, color="b",label="Original")
    #plt.scatter(X, mutant_prob, color="r",label="Mutant")
    ax.set(title="Impact of SNP on PTM Presence")
    ax.set(ylabel="probability")
    ax.set(xticklabels=[])
    ax.set(xlabel="SNP/PTM Pair")
    ax.tick_params(bottom=False)
    ax.axhline(0.5, linestyle='--', alpha=0.5)
    ax.axhline(0.6)
    ax.axhline(0.4)
    ax.axhspan(ymin=0.4, ymax=0.6, alpha=0.4, color='gray')
    #plt.title("Impact of SNP on PTM Presence")
    #plt.legend()
    plt.savefig("snp_ptm_probability.png", bbox_inches='tight', transparent=True)

CATALOG_PATH = 'Impacted SNPs catalog.xlsx'
def load_gwas_catalog(path=CATALOG_PATH):
    df = pd.read_excel(path, dtype='str')
    return df

def plot_disease_pie_chart(df):
    pie_data = df[['MAPPED_TRAIT', 'SNPS']].groupby(by='MAPPED_TRAIT').nunique().sort_values(by='SNPS', ascending=False)
    N = 6
    labels = list(pie_data.index[0:N]) + ['_nolegend_' for i in range(len(pie_data.index)-N)]
    import matplotlib.pyplot as plt
    #plt.figure(figsize=(20, 20))
    plot = pie_data.plot.pie(y='SNPS', labels=labels, labeldistance=None, figsize=(8, 10))
    plot.set(title='Disease Distribution')
