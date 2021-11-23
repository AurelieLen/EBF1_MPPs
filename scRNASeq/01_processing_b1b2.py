# -*- coding: utf-8 -*-
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
# !echo $SINGULARITY_CONTAINER
# !hostname

# %% [markdown]
# # Lenaerts et al. 2021 scRNA-Seq analysis

# %% [markdown]
# Data consists of 8 samples:
# - WT LK - wild type cells sorted from the LK gate
# - WT LSK - widl type cells sorted from the LSK gate
# - KO LK - Ebf1 KO cells sorted from the LK gate 
# - KO LSK - Ebf1 KO cells sorted from the LSK gate 
# In each case the sample is a 1:1 mixture of female+male cells

# %% [markdown]
# ## Setup

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import re
from utils import proc
from utils import load
import cellproject as cp
import pathlib

sc.set_figure_params(color_map='viridis', dpi_save=350)
sc.settings.verbosity = 3
# Setting up target directories
sc.settings.figdir = './figures/01script/'
base_figures = './figures/01script/'
base_procdata = './procdata/01script/'
for i in [sc.settings.figdir, base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# %%
# %load_ext autoreload
# %autoreload 2

# %%
# Matplotlib config
mpl.rcParams['figure.figsize'] = (4, 4)
mpl.rcParams['pdf.fonttype'] = 42 # Ensures readable fonts in illustrator
mpl.rcParams['ps.fonttype'] = 42
plt.rc('axes', axisbelow=True) # Ensure that gridlines are behind points
mplparams = mpl.rcParams.copy()


# %% [markdown]
# ## Loading data

# %% [markdown]
# ### Various gene-sets

# %%
# Loading a single dataset (needed for genenames etc)
adata = sc.read_10x_h5("./data/hpc_processed/" + 'WTLK_b1' + "_filtered_feature_bc_matrix.h5")
adata.var['symbol'] = adata.var.index
adata.var_names_make_unique()

# Getting y chromosome genes
ygenes = pd.read_csv('./data/Y_genes.csv')
realY_filter = [not bool(re.search("predicted gene", i)) for i in ygenes['Gene description']]
ygenes = ygenes.loc[realY_filter, :]
ygenes = ygenes.loc[ygenes['Gene name'].isin(adata.var.symbol),:]
ygenes = ygenes['Gene name']

# Getting cell cycle genes, which we typically remove from variable genes
ccgenes_toremove = pd.read_csv('./data/CC_genes_forremoval.txt', header=None)[0]
toremove = ccgenes_toremove.copy()
# Adding the Y chromosome genes and Xist (to avoid splitting by sex)
toremove = toremove.append(pd.Series(['Xist', 'Ebf1']))
toremove = toremove.append(ygenes)

# Getting cell cycle genes per phase
ccgenes = pd.read_csv('./data/macosko2015_hum_ccgenes.csv') # These are human cell cycle genes with just names swapped for mouse
ccgenes = ccgenes.loc[ccgenes['gene'].isin(adata.var.symbol),:]
Sgenes = ccgenes.loc[ccgenes['phase'] == 'S','gene']
G2Mgenes = ccgenes.loc[ccgenes['phase'] == 'G2_M', 'gene']

del adata

# %% [markdown]
# ### Loading counts, QC, pre-processing (each dataset)

# %% tags=[]
comb = load.load_data(data_keys=['WTLK_b1', 'WTLSK_b1', 'KOLK_b1', 'KOLSK_b1',
                                 'WTLK_b2', 'WTLSK_b2', 'KOLK_b2', 'KOLSK_b2'], ygenes=ygenes, testrun=False)
load.summarise_cellnos(comb)

# %% [markdown]
# ### Computing the landscape

# %%
sc.pp.filter_genes(comb, min_counts=1)
comb.layers['counts'] = comb.X.copy()

# %%
# Estimating variable genes
proc.process(comb,
             compute_hivar=True,
             n_pcs=50,
             n_neighbors=8,
             leiden_resolution=0.8,
             n_variable=7000,
             lognorm_data=True,  # !!!!!!CHECK!!!!!
             remove_genes=toremove,
             Sgenes=Sgenes,
             G2Mgenes=G2Mgenes,
             regress_cc=True,
             batch_correct='harmony',
             batch_key='experiment')

# %% [markdown]
# #### UMAP

# %% tags=[]
umapref = cp.quick_umap(comb, n_neighbors=8, use_rep='X_pca_harmony')
sc.tl.paga(comb)
load.plot_basicinfo(comb, savename='comb')
sc.pl.umap(comb, color=['Flt3', 'Satb1'], save='_markers_Flt3_Satb1.pdf')

comb.obsm['X_umap_2d'] = comb.obsm['X_umap'].copy()
umapref3d = cp.quick_umap(comb, n_neighbors=8,
                          n_components=3, use_rep='X_pca_harmony')
comb.obsm['X_umap_3d'] = comb.obsm['X_umap'].copy()
comb.obsm['X_umap'] = comb.obsm['X_umap_2d'].copy()

# %%
# Saving (substituting cell cycle regresed values for logn values to save space
comb.X = comb.raw[:, comb.var.index].X.copy()
del comb.raw
comb.write(base_procdata + 'aur4_comb_proc.h5ad', compression='lzf')
comb.obs.to_csv(base_procdata + 'aur4_comb_obs.csv')
