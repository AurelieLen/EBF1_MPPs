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
# # Lenaerts et al. 2021 analysis - hFL data analysis

# %% [markdown]
# ## Setup

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import tables as tb
import matplotlib as mpl
from utils import proc


sc.set_figure_params(color_map='viridis', dpi_save=350)
sc.settings.verbosity = 3
# Setting up target directories
sc.settings.figdir = './figures/07script/'
base_figures = './figures/07script/'
base_procdata = './procdata/07script/'
import pathlib
for i in [sc.settings.figdir, base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)
    
#Matplotlib config
mpl.rcParams['figure.figsize'] = (4, 4)
mpl.rcParams['pdf.fonttype'] = 42 #Ensures readable fonts in illustrator
mpl.rcParams['ps.fonttype'] = 42
plt.rc('axes', axisbelow=True) #Ensure that gridlines are behind points
mplparams = mpl.rcParams.copy()


from matplotlib.colors import LinearSegmentedColormap
cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, '#e8e8e8'),
                                                     (1, '#cc0404')])

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# ## Loading data

# %%
# hfl_BL_subset.h5ad contains all the blood cells + endothelium from the Popescu et al. 2019 publication
# All blood cells are in loaded from file below, which is subsequently subset to remove T cells/NK cells,
# Endothelium and late erythroid and myeloid
# hfl = sc.read('./data/ext_data/hFL_BL_QC_processed.h5ad')

# sc.tl.leiden(hfl, resolution=0.8)
# sc.pl.umap(hfl, color='leiden', legend_loc='on data')
# keep = ['0', '8', '18', '12', '23', '6', '10', '14', '19']
# hfl_small = hfl[hfl.obs['leiden'].isin(keep), :].copy()
# sc.pl.umap(hfl_small)
# hfl_small.write('./data/ext_data/hfl_BL_subset.h5ad', compression='lzf')

# %%
hfl_small = sc.read('./data/ext_data/hfl_BL_subset.h5ad')
hfl_small = hfl[hfl.obs['leiden'].isin(keep), :].copy()
proc.process(hfl_small,
             compute_hivar=True,
             n_pcs=50,
             n_neighbors=15,
             leiden_resolution=0.8,
             n_variable=5000,
             lognorm_data=False,
             batch_correct='harmony',
             batch_key='batch')
sc.pl.umap(hfl_small, color=['leiden', 'batch', 'Cell.Labels'], wspace = 1)

# %%
sc.pl.umap(hfl_small, color='leiden', legend_loc='on data', save='_hfl_subset_leiden.pdf')

# %%
#markers
sc.pl.umap(hfl_small,
           color=['PROCR', 'CD34', 'CD38', 'THY1'],
           cmap=cmap2, save='_hfl_subset_markers1.pdf')


# %% tags=[]
#markers
sc.pl.umap(hfl_small,
           color=['PROCR', 'CD34', 'CD38'],
           cmap=cmap2, save='_hfl_subset_markers1.pdf')

#Ly markers
sc.pl.umap(hfl_small,
           color=['FLT3', 'IL7R', 'DNTT', #Ly progenitos
                  'VPREB1', 'DNTT', 'SOX4', 'PAX5', #Pre-pro B cells
                  'CD19', #B cells
                  'MZB1'],
           cmap=cmap2,
          save='hfl_subset_Ly_markers.pdf') #margnal zone B cells

#Ery/Meg/Bas/MC markers
sc.pl.umap(hfl_small,
           color=['PF4', 'VWF', #Meg
                             'KLF1', 'GATA1', #Ery
                             'CMA1', 'TPSB2', #MC 
                             'PRG2', #Eos
                             'CLC'], #Bas
           cmap=cmap2,
          save='_hfl_subset_MegEryBasMC_markers.pdf')
           
#Myeloid markers
sc.pl.umap(hfl_small,
           color=['CD14', 'VCAN', 'CTSS', 'FCN1', 'CD163', #Class monocytes
                  'FCGR3A', 'LST1', 'MS4A7', 'SAT1', #non-class mono
                  'IRF8', 'IL3RA', 'SHD', 'CCDC50', 'DNTT', 'JCHAIN', #pDC progenitors
                  'IRF8', 'IL3RA', 'UGCG', 'CCDC50', 'CD4', 'TCF4', 'DERL3', #pDC
                  'S100A10', 'FCER1A', 'HAVCR2', 'CD2', 'C12orf75', #cDC1
                  'S100A10', 'FCER1A', 'CD1C', 'HLA-DQB1'],#cDC2
           cmap=cmap2,
          save='_hfl_subset_myomarkers.pdf')

# %%
anno_man = {'0' : 'DC',
            '1' : 'Ery prog',
            '2' : 'Pre/pro B cells',
            '3' : 'Bas/MC/Ery prog', 
            '4' : 'Ery prog',
            '5' : 'Mono class. ',
            '6' : 'HSC/MPP',
            '7' : 'Meg prog',
            '8' : 'Mono class.',
            '9' : 'Meg/Ery/Bas/MC prog',
            '10' : 'B cells',
            '11' : 'MC',
            '12' : 'Mono/DC prog',
            '13' : 'Meg/Ery/Bas/MC prog',
            '14' : 'Neu/Mono/DC prog',
            '15' : 'Meg prog',
           '16' : 'Ly/B prog',
           '17' : 'Mono class II',
           '18' : 'pDC',
           '19' : 'Bas',
           '20' : 'Pot. doublet',
           '21' : 'Mono non-class.'}
hfl_small.obs['anno_man'] = [anno_man[i] for i in hfl_small.obs.leiden]
sc.pl.umap(hfl_small, color='anno_man', save='_hfl_subset_annoman.pdf')

# %%
toplot = ['leiden', 'EBF1', 'CEBPA', 'PAX5', 'CD34', 'THY1', 'IKZF1', 'SPI1', 'CD19']
sc.pl.umap(hfl_small, color=toplot, save='_hfl_subset_chosengenes.pdf', cmap=cmap2)

# %%
toplot = ['leiden', 'EBF1', 'CEBPA', 'PAX5', 'CD34', 'THY1', 'IKZF1', 'SPI1', 'CD19']
sc.pl.umap(hfl_small[hfl_small.obs.leiden.isin(['6', '14', '16']), :], color=toplot, save='_hfl_subset_chosengenes_zoomin.pdf', cmap=cmap2)
