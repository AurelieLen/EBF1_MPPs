import scanpy as sc
import pandas as pd
import re
from . import proc

def load_data(data_keys, path='./data/hpc_processed/', ygenes=None, testrun=False):
    adatas = {key: sc.read_10x_h5(path + key + "_filtered_feature_bc_matrix.h5") for key in data_keys}
    if testrun:
        for i in adatas:
            sc.pp.subsample(adatas[i], n_obs=500)

    for x in adatas.items():
        key = x[0]
        val = x[1]

        print("Analysing sample: " + key)
        print(val)
        val.var['symbol'] = val.var.index
        val.var_names_make_unique()

        # QC
        proc.do_qc(val, min_genes=1200, mt_frac=0.1)

        # Assigning male/female to each sample
        val_logn = val.copy()
        proc.lognorm(val_logn)
        proc.assign_sex(val_logn, ygenes=ygenes, use_raw=False)

        #Removing male/female doublets and unassigned
        val.obs = val_logn.obs.copy()
        val = val[val.obs.index[val.obs.sex.isin(['female', 'male'])],:]

    df = pd.DataFrame({'sample' : [x for x in adatas],
                       'median_genes' : [x.obs.n_genes.median() for x in adatas.values() ],
                       'median_counts': [x.obs.n_counts.median() for x in adatas.values()]})
    print(df)

    comb = adatas[data_keys[0]]
    comb = comb.concatenate([adatas[i] for i in data_keys[1:]], batch_categories=data_keys)
    # Setting up necessary factors
    comb.obs['condition'] = [re.sub('(WT|KO)(LK|LSK)(.*)', '\\1', i) for i in comb.obs.batch]
    comb.obs['gate'] = [re.sub('(WT|KO)(LK|LSK)(.*)', '\\2', i) for i in comb.obs.batch]
    comb.obs['experiment'] = [re.sub('(WT|KO)(LK|LSK)_(.*)', '\\3', i) for i in comb.obs.batch]
    comb.obs['batch_sex'] = comb.obs.batch.str.cat(comb.obs.sex, sep='_')
    comb.obs['conditiongate_sex'] = comb.obs.condition + comb.obs.gate + '_' + comb.obs.sex

    # Extimating doublets with scrublet
    doublets = proc.doubletperrun(comb, 'batch')
    doublets['predicted_doublets_tr025'] = doublets.doublet_score > 0.25
    print('No of cells with doublet score above 0.25: ' + str(sum(doublets.predicted_doublets_tr025)))
    comb.obs = pd.concat((comb.obs, doublets.iloc[:,1:4]), axis=1)

    # Removing unassigned sex cells and female/male doublets
    comb = comb[comb.obs.index[comb.obs.sex.isin(['female', 'male'])],:].copy()

    return(comb)

def summarise_cellnos(adata):
    print('Final stats:')
    batch_sum = adata.obs.groupby(['batch'])[['n_genes', 'n_counts']].median()
    batch_sum['cellno'] = adata.obs.batch.value_counts()
    batch_sex_sum = adata.obs.groupby(['batch_sex'])[['n_genes', 'n_counts']].median()
    batch_sex_sum['cellno'] = adata.obs.batch_sex.value_counts()
    print(batch_sum)
    print(batch_sex_sum)
   

def plot_basicinfo(adata, savename='comb'):
#Umaps with basic information
    sc.pl.umap(adata, color=['sex', 'batch', 'experiment', 'batch_sex', 'n_genes', 'leiden', 'gate', 'phase', 'doublet_score'],
               alpha=0.6, wspace=0.6, save='_' + savename + '_info.pdf')
    sc.pl.umap(adata, color=['leiden'], alpha=0.6, save='_' + savename + '_leiden.pdf', legend_loc='on data')
    sc.pl.umap(adata, color=['Procr', 'Klf1', 'Pf4', 'Dntt', 'Irf8', 'Elane', 'Hoxb5', 'Mpo'],
               save='_' + savename + '_markers1.pdf')
    sc.pl.umap(adata, color=['Ly6d', 'Ebf1', 'Pax5', 'Cd19', 'Cd79a', 'Cd79b', 'Prtn3', 'Gata3', 'Cd3e'],
               save='_' + savename + '_markers2.pdf')
    sc.pl.umap(adata, color=['Ms4a2', 'Cma1', 'Gzmb', 'Prss34', 'Mcpt8', 'Prg2', 'Prg3'],
               save='_' + savename + '_markers3.pdf')
    sc.pl.umap(adata, color=['Mpo', 'Irf8', 'Irf4', 'Sirpa', 'Icam1', 'Siglech', 'Adgre1', 'Bst2',
                             'Csf1r', 'Ly6c1', 'Itgam', 'Spn', 'Sell', 'Itgax', 'Fcgr2b'],
               save='_' + savename + '_markers4.pdf')
    sc.pl.paga_compare(adata, threshold=0.1, save='_' + savename + '_paga.pdf')
