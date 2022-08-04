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
# # Lenaerts et al. 2021 analysis - various figures

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
import seaborn as sns
import cellproject as cp
import dotscore as dt
import plotnine as pn
from utils.HSCscore import calculate_HSCscore
from anndata import AnnData

sc.set_figure_params(color_map='viridis', dpi_save=350)
sc.settings.verbosity = 3
# Setting up target directories
sc.settings.figdir = './figures/05script/'
base_figures = './figures/05script/'
base_procdata = './procdata/05script/'
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
sfdata = sc.read('procdata/03script/sfdata_aur4fit.h5ad')

# %%
sc.pl.umap(aur4, color='leiden', legend_loc='on data', save='_aur4_leiden.pdf')
sc.pl.umap(aur4, color='condition', save='_aur4_condition.pdf')

# %% [markdown]
# ## Plotting projected cells from Nestorowa et al. 2016 paper

# %%
sfdata.obs['celltype_fin'] = sfdata.obs['celltype'].astype(str).copy()
sfdata.obs.loc[sfdata.obs.celltype_fin == 'HSC', 'celltype_fin'] = 'unassigned'
sfdata.obs.loc[sfdata.obs.celltype_e == 'ESLAM', 'celltype_fin'] = 'ESLAM'

# %%
fig, ax = plt.subplots()

fig.set_size_inches(5, 4)
sc.pl.umap(aur4, alpha =0.25, ax=ax, show=False)
x = sfdata[sfdata.obs.index[~sfdata.obs.celltype_fin.isin(['unassigned', 'MPP2', 'CMP', 'MEP', 'STHSC', 'MPP1'])],:]
x.uns['celltype_fin_colors'] = ["#85a100",
"#dbc400",
"#1272ba",
"#e87b2a"]
# from scanpy.plotting.palettes import default_20
# default_20[0:len(x.obs.celltype_fin.unique())]

sc.pl.umap(x,
           color='celltype_fin', ax=ax, show=False, size = 38, frameon=False)
# fig.tight_layout()

ax.set_box_aspect(1)

plt.savefig(base_figures + 'sfdata_aur4_proj_4cats.pdf')

# %%
fig, ax = plt.subplots()
fig.set_size_inches(5, 4)
sc.pl.umap(aur4, alpha =0.25, ax=ax, show=False)
x = sfdata[sfdata.obs.index[sfdata.obs.celltype_fin.isin(['ESLAM', 'STHSC'])],:]
x.uns['celltype_fin_colors'] = ["#85a100",
"#c453c9"]

sc.pl.umap(x,
           color='celltype_fin', ax=ax, show=False, size = 38, frameon=False)
ax.set_box_aspect(1)
plt.savefig(base_figures + 'sfdata_aur4_proj_HSCs.pdf')

# %% [markdown]
# ## HSC score

# %%
HSCscore_genes = pd.Series(np.genfromtxt('data/HSCscore/model_molo_genes.txt', dtype='str'))
missing = HSCscore_genes[~HSCscore_genes.isin(aur4.var.index)]

tempX = np.concatenate((aur4.layers['counts'].toarray(), np.zeros((aur4.n_obs, len(missing)))),
               axis = 1)
temp = AnnData(X = tempX,
              obs = aur4.obs,
              var = pd.DataFrame(index = np.concatenate((aur4.var.index.values, missing.values))))

aur4.obs['HSCscore'] = calculate_HSCscore(temp)
sc.pl.umap(aur4, color = 'HSCscore', save = '_HSCscore.pdf')

# %%
# Simple smoothing with nearest neighbors
x = aur4.obsp['connectivities']
inds = x.tolil().rows #Simple way to get the indices for nearest neighbors. lil is a lis of list format.
y = aur4.obs.HSCscore.values
nn_HSCscores = [y[i + [n]].mean() for n, i in enumerate(inds)] #indices are missing the cell itself so adding it back
aur4.obs['nn_HSCscore'] = nn_HSCscores

#Choosing the root
root = np.where(aur4.obs.nn_HSCscore == aur4.obs.nn_HSCscore.max())[0][0]

aur4.uns['iroot'] = root
aur4.obs['isroot'] = False
aur4.obs.isroot[root] = True
aur4.obs['isroot'] = pd.Categorical(aur4.obs.isroot)
sc.pl.umap(aur4, color = 'isroot')
from utils.plot import umap3d
umap3d(aur4, color = 'isroot', key = 'X_umap_3d', save='rootcheck.html')

#Computing pseudotime
sc.tl.dpt(aur4, n_dcs = 15)
sc.pl.umap(aur4, color='dpt_pseudotime', cmap='turbo', vmax=0.25, save='_dpt_vmax025.pdf')

# %% [markdown]
# ## Gene signature scors and Myc scores

# %%
toplot =  ['Odc1', 'Cdc20', 'Cyc1', 'Npm1', 'Cdk4', 'Eif3b', 'Myc']
vmaxes = aur4[:, toplot].X.todense().max(axis=0).A1
print(vmaxes)
vmaxes = list(np.round(1.05*vmaxes, 1))
print(vmaxes)

sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'WT'],:],
           color=toplot, save='_cc_growth_WT.pdf', cmap='viridis', size=7, vmax=vmaxes)
sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'KO'],:],
           color=toplot, save='_cc_growth_KO.pdf', cmap='viridis', size=7, vmax=vmaxes)

# %%
sc.pl.violin(aur4[aur4.obs.leiden == '4', :], keys=toplot, groupby='condition', save='_cc_growth.pdf')

# %%
sc.pl.umap(aur4, color=['Satb1', 'Flt3', 'Dntt'])

# %% [markdown]
# ## Ebf1 expression

# %%
sfdata_orig = sc.read('./data/sfdata/sfdata_nlog.h5ad')
sfdata_orig.obsm['X_umap'] = sfdata.obsm['X_umap'].copy()

# %%
from matplotlib.colors import LinearSegmentedColormap
cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, '#e8e8e8'),
                                                     (1, '#cc0404')])
sc.pl.umap(sfdata_orig, color='Ebf1', save='_sfdata_Ebf1.pdf', cmap=cmap2)

sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'WT'],:],
           color='Ebf1', save='_Ebf1.pdf', cmap=cmap2, size=5)

# %% [markdown]
# ## Lymphoid gene expression

# %%
lygenes = ['Igkc', 'Iglc1', 'Jchain', 'Ighm', 'Iglc1', 'Iglc2', 'Iglc3', 'Ebf1', 'Il7r', 'Cebpa']
aur4.obs['mouseid'] = aur4.obs.condition.astype(str) + '_' + aur4.obs.sex.astype(str) + '_' + aur4.obs.experiment.astype(str)

# %%
df = aur4.obs[['mouseid', 'leiden', 'condition']].join(pd.DataFrame(aur4[:,lygenes].X.toarray(), columns=lygenes, index=aur4.obs.index))
df['leiden'] = df.leiden.astype(str)
df = df.loc[df.leiden.isin(['0', '1', '2', '14', '4', '7']),:]
df = df.melt(id_vars=['mouseid', 'leiden', 'condition'], var_name='gene')
df = df.groupby(by=['leiden', 'mouseid', 'condition', 'gene']).mean()
df.reset_index(inplace=True)

df = df.melt(id_vars=['mouseid', 'leiden', 'condition', 'gene'], value_name='mean_expr')
df = df.dropna()

dfmean = df.groupby(['leiden', 'condition', 'gene']).mean()
dfmean.reset_index(inplace=True)
# dfsem = df.groupby(['leiden', 'condition', 'gene']).sem()
# dfsem.reset_index(inplace=True)

# %% tags=[]
df['leiden'] = pd.Categorical(df.leiden, categories=['0', '1', '2', '4', '7', '14'])

h1 = (pn.ggplot(df, pn.aes('leiden', 'mean_expr', color='condition'))
 + pn.geom_point(data=dfmean, mapping=pn.aes('leiden', 'mean_expr', fill='condition'),
                 position=pn.positions.position_dodge(width=0.7),
                shape='_', size=4)
 + pn.geom_point(stat='identity', position=pn.positions.position_dodge(width=0.7), alpha=0.5, size=1.7)
 + pn.theme_bw()
 + pn.facet_wrap('~gene', scales='free')
 + pn.ggtitle('Lymphoid gene expression by condition')
 +pn.theme(subplots_adjust={'wspace':0.3})
)
h1.save(base_figures + 'Lygenes_dotplots.pdf', width=16, height=6)
h1.draw()

# %% tags=[]
sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'WT'],:],
           color=lygenes, save='_Lygenes_WT.pdf', cmap=cmap2, size=7, vmax=7)
sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'KO'],:],
           color=lygenes, save='_Lygenes_KO.pdf', cmap=cmap2, size=7, vmax=7)

# %%
sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'WT'],:],
           color='Cebpa', save='_Lygenes_Cebpa_WT.pdf', cmap=cmap2, size=7, vmax=2.8)
sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'KO'],:],
           color='Cebpa', save='_Lygenes_Cebpa_KO.pdf', cmap=cmap2, size=7, vmax=2.8)

sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'WT'],:],
           color='Il7r', save='_Lygenes_Il7r_WT.pdf', cmap=cmap2, size=7, vmax=4.3)
sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'KO'],:],
           color='Il7r', save='_Lygenes_Il7r_KO.pdf', cmap=cmap2, size=7, vmax=4.3)

# %%
sc.pl.umap(aur4[aur4.obs.index[aur4.obs.condition == 'WT'],:],
           color=['Ly6d', 'Rag1', 'Dntt', 'Satb1'], save='_Lygenes2_WT.pdf', cmap=cmap2, size=5)

# %% [markdown]
# ## Common DE gene expression

# %%
genes = ['AY036118', 'Thg1l', 'mt-Nd3', 'H2-Ab1', 'Rpl36al', 'H2-D1', 'mt-Atp8', 'mt-Nd4l', 'Ly6c2', 'F13a1', 'Tap1', 'Nme2']

# %%
df = aur4.obs[['mouseid', 'leiden', 'condition']].join(pd.DataFrame(aur4[:,genes].X.toarray(), columns=genes, index=aur4.obs.index))
df['leiden'] = df.leiden.astype(str)
df = df.loc[df.leiden.isin(['0', '1', '2', '14', '4', '7']),:]
df = df.melt(id_vars=['mouseid', 'leiden', 'condition'], var_name='gene')
df = df.groupby(by=['leiden', 'mouseid', 'condition', 'gene']).mean()
df.reset_index(inplace=True)

df = df.melt(id_vars=['mouseid', 'leiden', 'condition', 'gene'], value_name='mean_expr')
df = df.dropna()

dfmean = df.groupby(['leiden', 'condition', 'gene']).mean()
dfmean.reset_index(inplace=True)
# dfsem = df.groupby(['leiden', 'condition', 'gene']).sem()
# dfsem.reset_index(inplace=True)

# %% tags=[]

h1 = (pn.ggplot(df, pn.aes('leiden', 'mean_expr', color='condition'))
 + pn.geom_point(data=dfmean, mapping=pn.aes('leiden', 'mean_expr', fill='condition'),
                 position=pn.positions.position_dodge(width=0.7),
                shape='_', size=4)
 + pn.geom_point(stat='identity', position=pn.positions.position_dodge(width=0.7), alpha=0.3)
 + pn.theme_bw()
 + pn.facet_wrap('~gene', scales='free')
 + pn.ggtitle('CommonDE gene expression by condition')
 +pn.theme(subplots_adjust={'wspace':0.3})
)
h1.save(base_figures + 'commonDE_genes_dotplots.pdf', width=16, height=6)
h1.draw()

# %% [markdown]
# ## Number of DE genes per cluster

# %%
deno = pd.read_csv('./procdata/02script/clusterDEpb_b1b2/DEsum.csv')

deno['clusterno'] = [int(x[7:]) for x in deno.cluster]
deno = deno.loc[deno.clusterno <=15,:]

deno = deno.melt(id_vars=['DEno', 'experiment', 'clusterno', 'cluster'])
deno['clusterno'] = pd.Categorical(deno['clusterno'], categories = np.sort(deno['clusterno'].unique()))

# %% tags=[]
h1 = (pn.ggplot(deno, pn.aes('clusterno', 'value'))
 + pn.geom_bar(stat='identity')
 + pn.theme_bw()
 + pn.facet_wrap('~variable', scales='free')
 + pn.ggtitle('Number of DE genes')
 +pn.theme(subplots_adjust={'wspace':0.3})
)
h1.save(base_figures + 'DEno.pdf', width=10, height=2)
h1.draw()

# %% [markdown]
# ## Enrichr and GSEA

# %%
cl4de = pd.read_csv('./procdata/02script/clusterDEpb_b1b2/DEb1b2_DEsig_cluster4.csv')

# %%
cl4_rnk = pd.DataFrame({0 : cl4de.genesymbol.str.upper(), 1 : -cl4de.logFC* np.log2(cl4de.pvalue)})

# %%
import gseapy as gp
from gseapy.plot import barplot, dotplot

# %%
cl4de_UP = cl4de.loc[cl4de.logFC > 0,'genesymbol'].to_list()
enr = gp.enrichr(gene_list=cl4de_UP,
                 gene_sets=['MSigDB_Hallmark_2020'],
                 organism='Mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                 description='test_name',
                 outdir='test/enrichr_msigH',
                 cutoff=0.1)
dotplot(enr.res2d, title='MSigDB_Hallmark_2020_UP',cmap='viridis_r')
dotplot(enr.res2d, title='MSigDB_Hallmark_2020_UP',cmap='viridis_r', ofname=base_figures + 'enrichr_msigDBH_cl4DEUP.pdf')

# %%
cl4de_DOWN = cl4de.loc[cl4de.logFC < 0,'genesymbol'].to_list()
enr = gp.enrichr(gene_list=cl4de_DOWN,
                 gene_sets=['MSigDB_Hallmark_2020'],
                 organism='Mouse', # don't forget to set organism to the one you desired! e.g. Yeast
                 description='test_name',
                 outdir='test/enrichr_msigH',
                 cutoff=0.1)
dotplot(enr.res2d, title='MSigDB_Hallmark_2020_UP',cmap='viridis_r')
dotplot(enr.res2d, title='MSigDB_Hallmark_2020_DOWN',cmap='viridis_r', ofname=base_figures + 'enrichr_msigDBH_cl4DEDOWN.pdf')

# %%
pre_res = gp.prerank(rnk=cl4_rnk, gene_sets='MSigDB_Hallmark_2020',
                     processes=4,
                     permutation_num=100, # reduce number to speed up testing
                     outdir='test/prerank_report_kegg', format='png', seed=6)


# %% [markdown]
# https://maayanlab.cloud/Enrichr/#libraries

# %%
pre_res = gp.prerank(rnk=cl4_rnk, gene_sets='GO_Biological_Process_2021',
                     processes=4,
                     permutation_num=100, # reduce number to speed up testing
                     outdir='test/prerank_report_kegg', format='png', seed=6)
