import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from .proc import cmap_RdBu2


def umap_subgroups(adata, key, toplot, subsample=None, alpha=0.5, file=None):

    n = len(toplot)
    nrows = -(-n // 4) # hack for ceil division
    print(nrows)
    fig, axes = plt.subplots(nrows=nrows, ncols=4)
    fig.set_size_inches(20, nrows * 5)
    axes = axes.ravel()
    
    for n, i in enumerate(toplot):
        sc.pl.umap(adata, show=False, ax=axes[n], size=16, alpha=alpha)
        # In the old version of scanpy there is no na_color and the palette
        # argument does not work for me, As a workaround copying the object
        # and plotting in default blue
        # In newer just set the na_color to a given value
        temp = adata[adata.obs.index[adata.obs[key] == i],:].copy()
        if subsample is not None:
            sc.pp.subsample(temp, n_obs=subsample)
        temp.obs['temp'] = 'a'
        sc.pl.umap(temp, color='temp', show=False, ax=axes[n], size=16,
                   legend_loc=None, title=i)

        axes[n].margins(0.05)

    fig.tight_layout()
    fig.show()
    if file is not None:
        plt.savefig(file)