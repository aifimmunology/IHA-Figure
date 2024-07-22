import pandas as pd
import numpy as np
import scanpy as sc


def hex_to_rgb(hex):
    """
    Converts a hexadecimal color code to an RGB tuple.

    Parameters:
    hex (str): A string representing a hexadecimal color code, prefixed with a `#`.

    Returns:
    tuple: A tuple containing three float values representing the RGB components of the color, each in the range [0, 1].

    """
    hex = hex.lstrip('#')
    return tuple(int(hex[i:i+2], 16)/255.0 for i in (0, 2, 4))

def create_cmap_from_hex(hex_list, cmap_name='custom_cmap'):
    """
    Creates a matplotlib colormap from a list of hexadecimal color codes.

    Parameters:
    hex_list (list): A list of strings, each representing a hexadecimal color code.
    cmap_name (str, optional): The name to assign to the colormap. Default is 'custom_cmap'.

    Returns:
    LinearSegmentedColormap: A colormap object that can be used in matplotlib for plotting.

    """
    rgb_list = [hex_to_rgb(hex_code) for hex_code in hex_list]
    cmap = mcolors.LinearSegmentedColormap.from_list(cmap_name, rgb_list)
    return cmap

def plot_nhood_graph(adata,alpha: float = 0.1,min_logFC: float = 0,cmap='RdBu_r',min_size: int = 10,plot_edges: bool = False,title: str = "DA log-Fold Change"):
    nhood_adata = adata.uns["nhood_adata"].copy()

    if "Nhood_size" not in nhood_adata.obs.columns:
        raise KeyError(
            'Cannot find "Nhood_size" column in adata.uns["nhood_adata"].obs -- \
                please run milopy.utils.build_nhood_graph(adata)'
        )

    nhood_adata.obs["graph_color"] = nhood_adata.obs["logFC"]
    nhood_adata.obs.loc[nhood_adata.obs["SpatialFDR"]
                        > alpha, "graph_color"] = np.nan
    nhood_adata.obs["abs_logFC"] = abs(nhood_adata.obs["logFC"])
    nhood_adata.obs.loc[nhood_adata.obs["abs_logFC"]
                        < min_logFC, "graph_color"] = np.nan

    # Plotting order - extreme logFC on top
    nhood_adata.obs.loc[nhood_adata.obs["graph_color"].isna(),
                        "abs_logFC"] = np.nan
    ordered = nhood_adata.obs.sort_values(
        'abs_logFC', na_position='first').index
    nhood_adata = nhood_adata[ordered]

    vmax = np.max([nhood_adata.obs["graph_color"].max(),
                  abs(nhood_adata.obs["graph_color"].min())])
    vmin = - vmax

    sc.pl.embedding(nhood_adata, "X_milo_graph",
                    color="graph_color", cmap=cmap,
                    size=adata.uns["nhood_adata"].obs["Nhood_size"]*min_size,
                    edges=plot_edges, neighbors_key="nhood",
                    # edge_width =
                    sort_order=False,
                    frameon=False,
                    vmax=vmax, vmin=vmin,
                    title=title)

    
def grouped_obs_sum_raw(adata_filt, group_key, layer=None, gene_symbols=None):
    """
    Compute the sum of raw counts for each group defined by `group_key` in the given AnnData object.

    Parameters:
    -----------
    adata_filt : AnnData
        Filtered AnnData object containing the single-cell data.
    group_key : str
        Key in `adata_filt.obs` to group the observations by.
    layer : str, optional
        If specified, the layer in `adata_filt` to use for calculations. If None, use `adata_filt.X`.
    gene_symbols : list of str, optional
        List of gene symbols to include in the sum. If None, include all genes.

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the summed observations for each group. Rows correspond to genes and columns to groups.
    """
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X

    if gene_symbols is not None:
        idx = adata_filt.var_names.isin(gene_symbols)
        new_idx = adata_filt.var_names[idx]
    else:
        new_idx = adata_filt.var_names

    grouped = adata_filt.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((len(new_idx), len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=new_idx
    )

    for group, idx in grouped.indices.items():
        X = getX(adata_filt[idx])
        out[group] = np.ravel(X.sum(axis=0, dtype=np.float64))

    return out

def grouped_obs_mean(adata_filt, group_key, layer=None, gene_symbols=None):
    """
    Compute the mean of log-normalized expression for each group defined by `group_key` in the given AnnData object.

    Parameters:
    -----------
    adata_filt : AnnData
        Filtered AnnData object containing the single-cell data.
    group_key : str
        Key in `adata_filt.obs` to group the observations by.
    layer : str, optional
        If specified, the layer in `adata_filt` to use for calculations. If None, use `adata_filt.X`.
    gene_symbols : list of str, optional
        List of gene symbols to include in the mean calculation. If None, include all genes.

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the mean observations for each group. Rows correspond to genes and columns to groups.
    """
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X

    if gene_symbols is not None:
        idx = adata_filt.var_names.isin(gene_symbols)
        new_idx = adata_filt.var_names[idx]
    else:
        new_idx = adata_filt.var_names

    grouped = adata_filt.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((len(new_idx), len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=new_idx
    )

    for group, idx in grouped.indices.items():
        X = getX(adata_filt[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))

    return out