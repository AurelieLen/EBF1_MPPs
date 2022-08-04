import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.sparse import issparse
import plotly as plotly
import plotly.graph_objs as go
from .proc import cmap_RdBu2
from pandas.api.types import is_numeric_dtype


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


######## 3d plots ########
def umap3d(adata,
           color,
           key='X_umap',
           colorscale='Viridis',
           save=None,
           cmid=None,
           cmax=None,
           cmin=None,
           use_raw=True):
    """Create colourcoded 3d UMAP plots with plotly.

    Outputs 3d plots (using plotly framework) of the UMAp representation. Cells
    are colourcoded according to either gene expression or a column in the adata.obs
    slot.

    Parameters
    ----------
    adata : AnnData object
        AnnData objects with gene expression values and UMAP coordinates
    color : str
        Name of a column in the adata.obs (will be treated as categorical) slot
        or name of the gene (in adata.var) to be used for colourcoding
    key : str
        Name of the dimensionality reduction coordinates (key in .obsm attribute)
        to be used, defaults to 'X_umap'
    colorscale
        Name of the matplotlib colorscale to use
    save : str
        If not none, string indicating file destination (uses sc.settings.figdir)
        otherwise the plot will be lanuched inside of jupyter notebook
    cmid
        value for the middle of colorscale
    cmax
        value for the maximum of colorscale
    cmin
        value for the minimum of colorscale
    use_raw
        bool, whether to use the .raw field

    Returns
    -------
    Either returns the figure (if no file is provided) or saves to the indicated file
    """

    if key == "X_diffmap":
        x, y, z = adata.obsm[key][:, 1], adata.obsm[key][:, 2], adata.obsm[key][:, 3]
    else:
        x, y, z = adata.obsm[key][:, 0], adata.obsm[key][:, 1], adata.obsm[key][:, 2]

    traces = []
    if color in adata.obs.columns:
        if hasattr(adata.obs[color], 'cat'):
            cats = adata.obs[color].cat.categories
            for n, i in enumerate(cats):

                xI = x[adata.obs[color] == i]
                yI = y[adata.obs[color] == i]
                zI = z[adata.obs[color] == i]

                if color + '_colors' in adata.uns.keys():
                    col = adata.uns[color + '_colors'][n]
                else:
                    if len(cats) <= 20:
                        col = default_20[n]
                    else:
                        col = godsnot_102[n]

                trace1 = go.Scatter3d(
                    x=xI,
                    y=yI,
                    z=zI,
                    mode='markers',
                    marker=dict(
                        size=5,
                        color=col
                    ),
                    text=adata.obs.index[adata.obs[color] == i],
                    name=str(i))

                traces.append(trace1)
                data = traces

        elif is_numeric_dtype(adata.obs[color]):
            col = adata.obs[color]
            traces = go.Scatter3d(
                x=x,
                y=y,
                z=z,
                text=col, #[str(x) for x in color],
                mode='markers',
                marker=dict(
                    size=5,
                    color=col,  # set color to an array/list of desired values
                    colorscale=colorscale)  # choose a colorscale
            )
            # assembling the traces into data
            data = [traces]

        else:
            raise Exception("The column in the .obs slot needs to be either category",
                  "int or float32/float64")

    elif (color in adata.var.index) or (color in adata.raw.var.index):

        if use_raw:
            col = adata.raw.var.index.get_loc(color)
            col = adata.raw.X[:, col]
        else:
            col = adata.var.index.get_loc(color)
            col = adata.X[:, col]

        if issparse(col):
            col = col.todense()
            col = col.A1

        traces = go.Scatter3d(
            x=x,
            y=y,
            z=z,
            text=col,
            mode='markers',
            marker=dict(
                size=5,
                color=col,  # set color to an array/list of desired values
                colorscale=colorscale,  # choose a colorscale
                )
            )
        # assembling the traces into data
        data = [traces]

    else:
        raise Exception("Color not found in adata.obs or in the adata.var")

    # Adjusting the colorscale
    # setting the center of the continuous color scale (for instance with negative and
    # positive values one often wants this at 0)
    if cmid is not None:
        traces.marker['cmid'] = cmid
    if cmin is not None:
        traces.marker['cmin'] = cmin
    if cmax is not None:
        traces.marker['cmax'] = cmax

    fig = go.Figure(data=data)  # layout=layout)

    if save is not None:
        plotly.offline.plot(fig, filename=str(sc.settings.figdir / ('umap3d' + save)), auto_open=False)
    else:
        return fig
