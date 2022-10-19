#!/usr/bin/env python
# 🐛

# Script with modified functions

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #
# Tweaked code for scglue get_metacells to return metadata for metacells  #
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * #

import os
from collections import defaultdict
from itertools import chain
from typing import Callable, List, Mapping, Optional

import anndata
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
import scipy.stats
import sklearn.cluster
import sklearn.decomposition
import sklearn.feature_extraction.text
import sklearn.linear_model
import sklearn.neighbors
import sklearn.preprocessing
import sklearn.utils.extmath
from anndata import AnnData
from networkx.algorithms.bipartite import biadjacency_matrix
from sklearn.preprocessing import normalize
from sparse import COO

from scglue import genomics, num
from scglue.typehint import Kws
from scglue.utils import logged, smart_tqdm

# Modified code to output a dictionary including metacell metadata
# See original at https://github.com/gao-lab/GLUE/blob/60e0655b14503de1741f9554c1de229c6e3f140d/scglue/data.py
def get_metacells(
        *adatas: AnnData, use_rep: str = None, n_meta: int = None,
        common: bool = True, seed: int = 0,
        agg_kwargs: Optional[List[Kws]] = None
) -> List[AnnData]:
    r"""
    Aggregate datasets into metacells
    Parameters
    ----------
    *adatas
        Datasets to be correlated
    use_rep
        Data representation based on which to cluster meta-cells
    n_meta
        Number of metacells to use
    common
        Whether to return only metacells common to all datasets
    seed
        Random seed for k-Means clustering
    agg_kwargs
        Keyword arguments per dataset passed to :func:`aggregate_obs`
    Returns
    -------
    adatas
        A list of AnnData objects containing the metacells
    Note
    ----
    When a single dataset is provided, the metacells are clustered
    with the dataset itself.
    When multiple datasets are provided, the metacells are clustered
    jointly with all datasets.
    """
    if use_rep is None:
        raise ValueError("Missing required argument `use_rep`!")
    if n_meta is None:
        raise ValueError("Missing required argument `n_meta`!")
    adatas = [
        AnnData(
            X=adata.X,
            obs=adata.obs.set_index(adata.obs_names + f"-{i}"), var=adata.var,
            obsm=adata.obsm, varm=adata.varm, layers=adata.layers
        ) for i, adata in enumerate(adatas)
    ]  # Avoid unwanted updates to the input objects

    print("Clustering metacells...")
    combined = anndata.concat(adatas)
    try:
        import faiss
        kmeans = faiss.Kmeans(
            combined.obsm[use_rep].shape[1], n_meta,
            gpu=False, seed=seed
        )
        kmeans.train(combined.obsm[use_rep])
        _, combined.obs["metacell"] = kmeans.index.search(combined.obsm[use_rep], 1)
    except ImportError:
        print(
            "`faiss` is not installed, using `sklearn` instead... "
            "This might be slow with a large number of cells. "
            "Consider installing `faiss` following the guide from "
            "https://github.com/facebookresearch/faiss/blob/main/INSTALL.md"
        )
        kmeans = sklearn.cluster.KMeans(n_clusters=n_meta, random_state=seed)
        combined.obs["metacell"] = kmeans.fit_predict(combined.obsm[use_rep])
    for adata in adatas:
        adata.obs["metacell"] = combined[adata.obs_names].obs["metacell"]

    print("Aggregating metacells...")
    agg_kwargs = agg_kwargs or [{}] * len(adatas)
    if not len(agg_kwargs) == len(adatas):
        raise ValueError("Length of `agg_kwargs` must match the number of datasets!")
    adatas = [
        scglue.data.aggregate_obs(adata, "metacell", **kwargs)
        for adata, kwargs in zip(adatas, agg_kwargs)
    ]
    if common:
        common_metacells = list(set.intersection(*(
            set(adata.obs_names) for adata in adatas
        )))
        if len(common_metacells) == 0:
            raise RuntimeError("No common metacells found!")
        return {'meta':combined.obs,'embedding':combined.obsm[use_rep],'adatas':[adata[common_metacells].copy() for adata in adatas]}
    return {'meta':combined.obs,'embedding':combined.obsm[use_rep],'adatas':adatas}
