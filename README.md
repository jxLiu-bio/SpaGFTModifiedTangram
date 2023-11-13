Tangram is a Python package, written in PyTorch and based on scanpy, for mapping single-cell (or single-nucleus) gene expression data onto spatial gene expression data. The single-cell dataset and the spatial dataset should be collected from the same anatomical region/tissue type, ideally from a biological replicate, and need to share a set of genes. Tangram_SpaGFT, is an improved version for orignal Tangram by adding the frequency constraints.

## How to install Tangram_SpaGFT?
SpaGCN_SpaGFT is based on Tangram.  To install  it, make sure you have [PyTorch](https://pytorch.org/) and [scanpy](https://scanpy.readthedocs.io/en/stable/) installed.

Create a conda environment
```
conda create -n tangram_spagft_env python==3.8.5
```
and activate the environment
```
conda activate tangram_spagft_env
```

Next, install Tangram_SpaGFT
```
python setup.py install
```

## Run Tangram_SpaGFT.
To begin with, import following modules in python environment.
```
import scanpy as sc
import Tangram_SpaGFT as tg
```

Load your spatial data and your single cell data (which should be in AnnData format), and pre-process them using tg.pp_adatas:
```
ad_sp = sc.read_h5ad(path)
ad_sc = sc.read_h5ad(path)
tg.pp_adatas(ad_sc, ad_sp, genes=None)
```

The function `pp_adatas` finds the common genes between `adata_sc`, `adata_sp`, and saves them in two `adatas.uns` for mapping and analysis later. Also, it subsets the intersected genes to a set of training genes passed by genes. If genes=None, Tangram maps using all genes shared by the two datasets. Once the datasets are pre-processed we can map:
```
ad_map = tg.map_cells_to_space_gft(ad_sc, ad_sp)
```

The returned AnnData,`ad_map`, is a cell-by-voxel structure where `ad_map.X[i, j]` gives the probability for cell `i` to be in voxel `j`. This structure can be used to project gene expression from the single cell data to space, which is achieved via `tg.project_genes`. 
```
ad_ge = tg.project_genes(ad_map, ad_sc)
```
The returned `ad_ge` is a voxel-by-gene AnnData, similar to spatial data `ad_sp`, but where gene expression has been projected from the single cells. This allows to extend gene throughput, or correct for dropouts, if the single cells have higher quality (or more genes) than spatial data. It can also be used to transfer cell types onto space. 

### How to run Tangram_SpaGFT at cluster level
Prepare the input data as the same you would do for cell level Tangram mapping. Then map using following code:

```
    ad_map = tg.map_cells_to_space_gft(
                   ad_sc, 
                   ad_sp,         
                   mode='clusters',
                   cluster_label='subclass_label')
```

Provided cluster_label must belong to ad_sc.obs. Above example code is to map at 'subclass_label' level, and the 'subclass_label' is in ad_sc.obs.

To project gene expression to space, use `tg.project_genes` and be sure to set the `cluster_label` argument to the same cluster label in mapping.

```
    ad_ge = tg.project_genes(
                  ad_map, 
                  ad_sc,
                  cluster_label='subclass_label')
```