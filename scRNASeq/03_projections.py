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

# %% [markdown]
# # Lenaerts et al. 2021 analysis - cell projections

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
from utils.plot import umap_subgroups
import seaborn as sns
import cellproject as cp

sc.set_figure_params(color_map='viridis', dpi_save=350)
sc.settings.verbosity = 3
# Setting up target directories
sc.settings.figdir = './figures/03script/'
base_figures = './figures/03script/'
base_procdata = './procdata/03script/'
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
sc.pl.umap(aur4, color='leiden', legend_loc='on data')

# %%
sfdata = sc.read('./data/sfdata/sfdata_nlog.h5ad')

# %% [markdown]
# Creating transient reference data with logn values

# %% [markdown]
# ## Projecting Nestorowa et al. 2016 data onto the landscape

# %%
#Unifying genes
ref = aur4.copy()

ref = ref[:, ref.var.index[ref.var.index.isin(sfdata.var.index)]].copy()
ref.X = ref[:,ref.var.index].X.copy()
sfdata = sfdata[:, ref.var.index].copy()

comb = ref.concatenate(sfdata, batch_key='data_type', batch_categories=['10x', 'SS2'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-10x']

#Subsetting just for variable genes to speed up
comb = comb[:, comb.var.index[comb.var.highly_variable]].copy()

sfdata = sfdata[:,comb.var.index].copy()
ref = ref[:,comb.var.index].copy()

# %%
#Running the batch correction with Seurat
comb_cor = cp.run_SeuratCCA(comb)

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref, n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
sfdata.X = comb_cor[sfdata.obs.index, sfdata.var.index].X.copy()

#Projecting into the ref PC space and identifying neighbors
cp.project_cells(sfdata, ref,
                 obs_columns=['leiden'],
                 fit_pca=True,
                 scale_data=True)

# %%
toplot = []
for i in sfdata.obs.celltype.unique():
    x = sfdata[sfdata.obs.index[sfdata.obs.celltype == i],:].copy()
    nn = sfdata.uns['cross_nn'][sfdata.obs.celltype == i,:]
    ascore = nn.sum(axis=0).A1 / nn.shape[0]
    ref.obs[i] = ascore
    toplot.append(i)
    
from matplotlib.colors import LinearSegmentedColormap
cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, '#DDD3D3'),
                                                    (0.001, 'tan'),
                                                    (1, 'blue')])
sc.pl.umap(ref, color = toplot, cmap=cmap2, save='sfdata_projected_nns.pdf')

# %% [markdown]
# Bringing back the harmony PCs

# %%
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()

cp.nnregress(sfdata,
             ref,
             regress=['pca'],
             weighted=True)

# %%
umapref = cp.quick_umap(ref, use_rep='X_pca', n_neighbors=8)

# %%
sfdata.obsm['X_umap'] = cp.quick_umap_proj(sfdata, umap_ref=umapref, rep_n_components=50)

# %%
sc.pl.umap(sfdata, color='celltype', save='umap_sfdataCP_celltype.pdf')
umap_subgroups(sfdata, key='celltype', toplot=sfdata.obs.celltype.unique(), file='umapsubgroups_sfdataCP_celltype.pdf')

# %%
sfdata.write(base_procdata + 'sfdata_aur4fit.h5ad', compression='lzf')

# %% [markdown]
# ## Plotting gene expression

# %%
sfdata_orig = sc.read('./data/sfdata/sfdata_nlog.h5ad')
sfdata_orig.obsm['X_umap'] = sfdata.obsm['X_umap'].copy()

# %%
from matplotlib.colors import LinearSegmentedColormap
cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, '#e8e8e8'),
                                                     (1, '#cc0404')])

# %%
sc.pl.umap(sfdata_orig, color='Ebf1', save='_sfdata_Ebf1.pdf', cmap=cmap2)

# %%
sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'WT'],:],
           color='Ebf1', save='_Ebf1.pdf', cmap=cmap2, size=5)

# %%
# Plotting lymphoid genes
lygenes = ['Igkc', 'Iglc1', 'Jchain', 'Ighm']
sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'WT'],:],
           color=lygenes, save='_Lygenes_WT.pdf', cmap=cmap2, size=5)
sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'KO'],:],
           color=lygenes, save='_Lygenes_KO.pdf', cmap=cmap2, size=5)
sc.pl.umap(sfdata_orig, color=lygenes, save='_Lygenes.pdf', cmap=cmap2, size=50)
