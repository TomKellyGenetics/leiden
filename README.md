# Leiden Algorithm

## leiden version 0.2.3

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/leiden)](https://cran.r-project.org/package=leiden)
[![Travis Build Status](https://travis-ci.org/TomKellyGenetics/leiden.svg?branch=master)](https://travis-ci.org/TomKellyGenetics/leiden)
[![CircleCI](https://circleci.com/gh/TomKellyGenetics/leiden.svg?style=svg)](https://circleci.com/gh/TomKellyGenetics/leiden)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/TomKellyGenetics/leiden?branch=master&svg=true)](https://ci.appveyor.com/project/TomKellyGenetics/leiden)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/TomKellyGenetics/leiden/branch/master/graph/badge.svg)](https://codecov.io/gh/TomKellyGenetics/leiden)

## Clustering with the Leiden Algorithm in R

This package allows calling the Leiden algorithm for clustering on an igraph object from R. See the Python and Java implementations for more details: 

https://github.com/CWTSLeiden/networkanalysis

https://github.com/vtraag/leidenalg

## Install

This package requires the 'leidenalg' and 'igraph' modules for python (2) to be installed on your system. For example:

```
pip install leidenalg numpy igraph
```

If you do not have root access, you can use `pip install --user` or `pip install --prefix` to install these in your user directory (which you have write permissions for) and ensure that this directory is in your PATH so that Python can find it.

Dependancies can also be installed from a conda repository. This is recommended fro Windows users:

```
conda -c vtraag python-igraph leidenalg
```

The stable 'leiden' package and the dependancies can be installed from CRAN:

```R
install.packages("leiden")
```

### Development version

The 'devtools' package can also be used to install development version of 'leiden' and the dependancies (igraph and reticulate) from GitHub:

```R
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("TomKellyGenetics/leiden")
```

## Usage

This package provides a function to perform clustering with the Leiden algorithm:

```R
partition <- leiden(adjacency_matrix)
```

For an igraph object 'graph' in R:

```R
adjacency_matrix <- igraph::as_adjacency_matrix(graph)
partition <- leiden(adjacency_matrix)
```

To use Leiden with the Seurat pipeline for a Seurat Object `object` that has an SNN computed (for example with `Seurat::FindClusters` with `save.SNN = TRUE`). This will compute the Leiden clusters and add them to the Seurat Object Class.

```R
adjacency_matrix <- as.matrix(object@snn)
partition <- leiden(adjacency_matrix)
object@ident <- as.factor(partition)
names(test@ident) <- rownames(test@meta.data)
object@meta.data$ident <- as.factor(partition)
```

Note that this code is designed for Seurat version 2 releases. For Seurat version 3 objects, the Leiden algorithm will be implemented in the Seurat version 3 package with `Seurat::FindClusters` and `algorithm = "leiden"`).  

These clusters can then be plotted with:

```R
library("RColorBrewer")
colourPal <- function(groups) colorRampPalette(brewer.pal(min(length(names(table(groups))), 11), "Set3"))(length(names(table(groups))))

object <- RunPCA(object = object, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = object, colors.use = colourPal(object@ident), group.by = "ident")

object <- RunTSNE(object = object, dims.use = 1:20, do.fast = TRUE, dim.embed = 2)
TSNEPlot(object = object, colors.use = colourPal(object@ident), group.by = "ident")

object <- RunUMAP(object = object, dims.use = 1:20, metric = "correlation", max.dim = 2)
DimPlot(object, reduction.use = "umap", colors.use = colourPal(object@ident), group.by = "ident")
```

### Example

```
#generate example data
adjacency_matrix <- rbind(cbind(matrix(round(rbinom(4000, 1, 0.8)), 20, 20), matrix(round(rbinom(4000, 1, 0.3)), 20, 20), matrix(round(rbinom(400, 1, 0.1)), 20, 20)),
##'                           cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20), matrix(round(rbinom(400, 1, 0.8)), 20, 20), matrix(round(rbinom(4000, 1, 0.2)), 20, 20)),
##'                           cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20), matrix(round(rbinom(4000, 1, 0.1)), 20, 20), matrix(round(rbinom(4000, 1, 0.9)), 20, 20)))
library("igraph")
rownames(adjacency_matrix) <- 1:60
colnames(adjacency_matrix) <- 1:60
graph_object <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
#plot graph structure
library("devtools")
install_github("TomKellyGenetics/igraph.extensions")
library("plot.igraph")
plot_directed(graph_object, cex.arrow = 0.3, col.arrow = "grey50")
#generate partitions
partition <- leiden(adjacency_matrix)
table(partition)
#plot results
library("RColorBrewer")
node.cols <- brewer.pal(max(partition),"Pastel1")[partition]
plot_directed(graph_object, cex.arrow = 0.3, col.arrow = "grey50", fill.node = node.cols)
```

<img src="https://github.com/TomKellyGenetics/leiden/blob/master/images/example_plot.png?raw=true" alt="A graph plot of results showing distinct clusters" width="600px">
</img>

### Vignette

For more details see the follow vignette:

[https://github.com/TomKellyGenetics/leiden/blob/master/vignettes/run_leiden.html](https://github.com/TomKellyGenetics/leiden/blob/master/vignettes/run_leiden.html)


### Citation

Please cite this implementation R in if you use it:

```
To cite the leiden package in publications use:

  S. Thomas Kelly (2019). leiden: R implementation of the Leiden algorithm. R
  package version 0.2.3 https://github.com/TomKellyGenetics/leiden

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {leiden: R implementation of the Leiden algorithm},
    author = {S. Thomas Kelly},
    year = {2018},
    note = {R package version 0.2.3},
    url = {https://github.com/TomKellyGenetics/leiden},
  }
 ```

Please also cite the original publication of this algorithm.

```
Traag, V.A., Waltman. L., Van Eck, N.-J. (2018). From Louvain to
       Leiden: guaranteeing well-connected communities.
       `arXiv:1810.08473 <https://arxiv.org/abs/1810.08473>`
```
