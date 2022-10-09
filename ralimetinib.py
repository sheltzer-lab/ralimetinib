#%%
import os

import numpy as np
import pandas as pd
import wget

#%%


data_dir = 'raw-data'
outdir = 'outdir'
def download_file(url, name, ddir=data_dir):
    path = os.path.join(ddir, name)
    wget.download(url, path)


def download_data():
     os.makedirs(data_dir, exist_ok=True)
     download_file('https://ndownloader.figshare.com/files/34990036', 'depmap_crispr_gene_effect.csv')
     download_file('https://ndownloader.figshare.com/files/17741420', 'primary_screen_lfc.csv')
     download_file('https://ndownloader.figshare.com/files/20237712', 'primary_screen_treatment_info.csv')
     download_file('https://storage.googleapis.com/public-smith-sheltzer-cancer-analysis/ralimetinib/depmap_rnai.csv', 'depmap_rnai.csv')

def load_data(ddir=data_dir):
    drugs = pd.read_csv(os.path.join(ddir, 'primary_screen_treatment_info.csv'), index_col=1)
    crispr = pd.read_csv(os.path.join(ddir, 'depmap_crispr_gene_effect.csv'), index_col=0)
    rnai = pd.read_csv(os.path.join(ddir, 'depmap_rnai.csv'), index_col=0)
    rnai = rnai.loc[:, rnai.isna().sum() < 250]

    ps = pd.read_csv(os.path.join(ddir, 'primary_screen_lfc.csv'), index_col=0)
    ps.columns = ps.columns.str.split('::', expand=True).get_level_values(0)


    return drugs, crispr, rnai, ps


def corrs(broad, psc, df, modality, name):
    os.makedirs(outdir, exist_ok=True)
    ral_lfc = psc[broad]
    df_corrs = df.corrwith(ral_lfc).sort_values(ascending=False)
    df_corrs.to_csv(os.path.join(outdir, name + '_' + modality + '_corrs.csv'))
    return df_corrs

download_data()
drugs, crispr, rnai, ps = load_data()

ralimetinib = 'LY2228820'
ralimetinib_broad = 'BRD-K59227464-334-02-4'
ralimetinib_metadata = drugs[drugs.name == ralimetinib]

corrs(ralimetinib_broad, ps, crispr, 'crispr', 'ralimetinib')
corrs(ralimetinib_broad, ps, rnai, 'rnai', 'ralimetinib')

gefitinib_broad = 'BRD-K64052750-001-17-5'
corrs(gefitinib_broad, ps, crispr, 'crispr', 'gefitinib')
corrs(gefitinib_broad, ps, rnai, 'rnai', 'gefitinib')

vemurafenib_broad = 'BRD-K56343971-001-10-6'
corrs(vemurafenib_broad, ps, crispr, 'crispr', 'vemurafenib')
corrs(vemurafenib_broad, ps, rnai, 'rnai', 'vemurafenib')

alpelisib_broad = 'BRD-K54997624-001-06-0'
corrs(alpelisib_broad, ps, crispr, 'crispr', 'alpelsib')
corrs(alpelisib_broad, ps, rnai, 'rnai', 'alpelsib')

