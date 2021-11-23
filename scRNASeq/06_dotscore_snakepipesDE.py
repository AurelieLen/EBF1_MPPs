# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
#Singularity container used:
# !echo $SINGULARITY_CONTAINER
# !hostname

# %% [markdown]
# # Lenaerts et al. 2021 analysis - DoT score analysis (bulk RNA-Seq data)

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import tables as tb
import matplotlib as mpl
from anndata import AnnData
import dotscore

#Matplotlib config
mpl.rcParams['pdf.fonttype'] = 42 #Ensures readable fonts in illustrator
mpl.rcParams['ps.fonttype'] = 42
plt.rc('axes', axisbelow=True) #Ensure that gridlines are behind points

# %%
mpl.rcParams['figure.figsize'] = (4, 4) #1:1 plot ratio
mpl.rc('xtick', labelsize=14) 
mpl.rc('ytick', labelsize=14) 
mpl.rc('axes', labelsize=16) 
mpl.rc('axes', labelsize=16) 

# %%
sc.settings.figdir = './figures/06script/'
base_figures = './figures/06script/'
base_procdata = './procdata/06script/'
import pathlib
for i in [sc.settings.figdir, base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# ## Loading data

# %%
LKadata = sc.read('./data/LKdata/LKLSKdata_QC_processed_small_annotated.h5ad')
LKadata.obsm['X_umap'] = LKadata.obsm['X_umap_noG']
sc.pl.umap(LKadata, color = 'leiden_noG', legend_loc = 'on data', save = '_LKdata_clusters.pdf')

# Subsetting for running
# sc.pp.subsample(LKadata, n_obs = 4000)

# %% [markdown]
# ## Calculating DoTscore - setting origins scales

# %%
scaling_clusters = ['3', '19']
means_cluster = LKadata[LKadata.obs.index[LKadata.obs.leiden_noG.isin(scaling_clusters)],:].X.mean(axis = 0)
dotscore.custom_scale(LKadata, mean = means_cluster)

# %% [markdown]
# ### MPP3 population - Ebf1 KO vs WT comparison

# %%
mpp3 = pd.read_csv('./data/DE_snakepipes/DEseq_results_MPP3_WTvsKO.txt', delimiter='\t')
mpp3['log2FoldChange'] = -mpp3['log2FoldChange'] # Need the KO vs WT comparison

mpp3 = mpp3.loc[~mpp3.padj.isna(),:]
allgenes = mpp3.gene.values
mpp3DE = mpp3.loc[mpp3.padj < 0.1,:]
allfolds = mpp3DE.log2FoldChange.values

# %%
score = dotscore.get_DoTscore(LKadata, 
                              de = mpp3DE, 
                              allgenes = allgenes, 
                              allfolds = allfolds,
                              zscore = True, 
                              id_col = 'gene',
                              simno = 2000)
score = dotscore.qfilt(score, 0.995, 0.005)
LKadata.obs['MPP3_DoT'] = score
c1 = dotscore.cmap_RdBu(LKadata.obs['MPP3_DoT'])
sc.pl.umap(LKadata, color=['MPP3_DoT'], alpha = 0.7, color_map = c1, frameon = False, save = '_LKdata_MPP3_dotscore_snakeDE.pdf')

# %%
leiden_scores = dotscore.get_genescore_pergroup(LKadata,
                                                mpp3DE,
                                                group = 'leiden_noG',
                                                sortby = '3',
                                                gene_symbols = 'symbol',
                                                id_col = 'gene',
                                         weight_col = 'log2FoldChange')
leiden_scores = leiden_scores.loc[leiden_scores.sum(axis = 1) != 0, :]
leiden_scores.to_csv(base_procdata + 'MPP3_EbfKO_leiden_scores.csv')

# %%
x = leiden_scores['8'].sort_values()
plt.scatter(range(len(x)), np.cumsum(x))

# %% [markdown]
# ### MPP4 population - Ebf1 KO vs WT comparison

# %%
mpp4 = pd.read_csv('./data/DE_snakepipes/DEseq_results_MPP4_WTvsKO.txt', delimiter='\t')
mpp4['log2FoldChange'] = -mpp4['log2FoldChange'] # Need the KO vs WT comparison

mpp4 = mpp4.loc[~mpp4.padj.isna(),:]
allgenes = mpp4.gene.values
mpp4DE = mpp4.loc[mpp4.padj < 0.1,:]
allfolds = mpp4DE.log2FoldChange.values

# %%
score = dotscore.get_DoTscore(LKadata, 
                              de = mpp4DE, 
                              allgenes = allgenes, 
                              allfolds = allfolds,
                              zscore = True, 
                              id_col = 'gene',
                              simno = 2000)
score = dotscore.qfilt(score, 0.995, 0.005)
LKadata.obs['MPP4_DoT'] = score
c1 = dotscore.cmap_RdBu(LKadata.obs['MPP4_DoT'])
sc.pl.umap(LKadata, color=['MPP4_DoT'], alpha = 0.7, color_map = c1, frameon = False, save = '_LKdata_MPP4_dotscore_snakeDE.pdf')

# %%
leiden_scores = dotscore.get_genescore_pergroup(LKadata,
                                                mpp4DE,
                                                group = 'leiden_noG',
                                                sortby = '3',
                                                gene_symbols = 'symbol',
                                                id_col = 'gene',
                                         weight_col = 'log2FoldChange')
leiden_scores = leiden_scores.loc[leiden_scores.sum(axis = 1) != 0, :]
leiden_scores.to_csv(base_procdata + 'MPP4_EbfKO_leiden_scores.csv')

# %%
x = leiden_scores['8'].sort_values()
plt.scatter(range(len(x)), np.cumsum(x))
