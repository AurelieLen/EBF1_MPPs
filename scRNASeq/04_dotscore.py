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
# # Lenaerts et al. 2021 analysis - DoT score

# %% [markdown]
# ## Setup

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import tables as tb
import matplotlib as mpl
import re
import cellproject as cp
import dotscore as dt


sc.set_figure_params(color_map = 'viridis', dpi_save = 350)
sc.settings.verbosity = 3
# Setting up target directories
sc.settings.figdir = './figures/04script/'
base_figures = './figures/04script/'
base_procdata = './procdata/04script/'
import pathlib
for i in [sc.settings.figdir, base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)
    
#Matplotlib config
mpl.rcParams['figure.figsize'] = (4, 4)
mpl.rcParams['pdf.fonttype'] = 42 #Ensures readable fonts in illustrator
mpl.rcParams['ps.fonttype'] = 42
plt.rc('axes', axisbelow=True) #Ensure that gridlines are behind points
mplparams = mpl.rcParams.copy()

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# ## Loading data

# %%
aur4 = sc.read('./procdata/01script/aur4_comb_proc.h5ad')
aur4 = aur4[aur4.obs.index[aur4.obs.condition == 'WT'], :].copy()
# sc.pp.subsample(aur4, n_obs = 1000)
sc.pl.umap(aur4, color='leiden', legend_loc='on data')

# %%
DE = pd.read_csv('./procdata/02script/clusterDEpb_b1b2/DEallclusters_sig.csv')

# %%
cl4DE = DE.loc[DE.cluster=='cluster4',:]

# %% [markdown]
# ### Scaling reference data

# %%
means_cluster = aur4[aur4.obs.index[aur4.obs.leiden == '4'],:].X.mean(axis = 0)
dt.custom_scale(aur4, mean = means_cluster)

# %% [markdown]
# ### DoTscore

# %%
allgenes = aur4.var.index.values.copy()
allfolds = cl4DE.logFC.values.copy()

score = dt.get_DoTscore(aur4, 
                              de = cl4DE, 
                              allgenes = allgenes, 
                              allfolds = allfolds,
                              zscore = True,
                              id_col = 'genesymbol',
                              weight_col='logFC',
                              simno = 1000)
aur4.obs['cl4_Ebf1_dtscore'] = dt.qfilt(score, 0.995, 0.005)

# %%
c1 = dt.cmap_RdBu(aur4.obs['cl4_Ebf1_dtscore'])

sc.pl.umap(aur4, color = 'cl4_Ebf1_dtscore', cmap = c1, save = 'cl4_Ebf1KO_dotscore.pdf')

# %%
leiden_scores = dt.get_genescore_pergroup(aur4,
                                                cl4DE,
                                                group = 'leiden',
                                                sortby = '3',
                                                gene_symbols = 'symbol',
                                                id_col = 'genesymbol',
                                         weight_col = 'logFC')
leiden_scores = leiden_scores.loc[leiden_scores.sum(axis = 1) != 0, :]
leiden_scores.to_csv(base_procdata + 'cl4_EbfKO_leiden_scores.csv')


# %% [markdown]
# Looking at cummulative weights in example clusters

# %%
x = leiden_scores['3'].sort_values()
plt.scatter(range(len(x)), np.cumsum(x))

# %%
x = leiden_scores['8'].sort_values()
plt.scatter(range(len(x)), np.cumsum(x))

# %%
x = leiden_scores['11'].sort_values()
plt.scatter(range(len(x)), np.cumsum(x))

# %%
