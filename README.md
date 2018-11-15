

https://github.com/CWTSLeiden/networkanalysis

https://github.com/vtraag/leidenalg

## Install

This package requires the 'leidenalg' and 'igraph' modules for python (2) to be installed on your system. For example:

``pip install leidenalg igraph``

If you do not have root access, you can use `pip install --user` or `pip install --prefix` to install these in your user directory (which you have write permissions for) and ensure that this directory is in your PATH so that Python can find it.

The 'devtools' package will be used to install 'leiden' and the dependancies (igraph and reticulate):

```R
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("TomKellyGenetics/leiden")
```

## Usage

This package provides a function to perform clustering with the Leiden algorithm:

```R
membership <- leiden(adjacency_matrix)
```

For an igraph object 'graph' in R:

```R
adjacency_matrix <- igraph::as_adjacency_matrix(graph)
membership <- leiden(adjacency_matrix)
```

To use Leiden with the Seurat pipeline for a Seurat Object `object` that has an SNN computed (for example with `Seurat::FindClusters` with `save.SNN = TRUE`). This will compute the Leiden clusters and add them to the Seurat Object Class.

```R
adjacency_matrix <- as.matrix(object@snn)
membership <- leiden(adjacency_matrix)
object@ident <- membership
object@meta.data$ident <- membership
```


### Citation

Traag, V.A., Waltman. L., Van Eck, N.-J. (2018). From Louvain to
       Leiden: guaranteeing well-connected communities.
       `arXiv:1810.08473 <https://arxiv.org/abs/1810.08473>`_
