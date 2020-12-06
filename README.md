# Leiden Algorithm

## leiden development version 0.4.0

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/leiden)](https://cran.r-project.org/package=leiden)
[![Travis Build Status](https://travis-ci.org/TomKellyGenetics/leiden.svg?branch=master)](https://travis-ci.org/TomKellyGenetics/leiden)
[![CircleCI](https://circleci.com/gh/TomKellyGenetics/leiden.svg?style=svg)](https://circleci.com/gh/TomKellyGenetics/leiden)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/TomKellyGenetics/leiden?branch=master&svg=true)](https://ci.appveyor.com/project/TomKellyGenetics/leiden)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![codecov](https://codecov.io/gh/TomKellyGenetics/leiden/branch/master/graph/badge.svg)](https://codecov.io/gh/TomKellyGenetics/leiden)
[![Downloads](https://cranlogs.r-pkg.org/badges/leiden)](https://CRAN.R-project.org/package=leiden)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/leiden?color=orange)](https://CRAN.R-project.org/package=leiden)
[![GitHub Views](http://hits.dwyl.com/tomkellygenetics/leiden.svg)](http://hits.dwyl.com/tomkellygenetics/leiden)

## Clustering with the Leiden Algorithm in R

This package allows calling the Leiden algorithm for clustering on an igraph object from R. See the Python and Java implementations for more details: 

[https://github.com/CWTSLeiden/networkanalysis](https://github.com/CWTSLeiden/networkanalysis)

[https://github.com/vtraag/leidenalg](https://github.com/vtraag/leidenalg)

## Install

Note: this is the _development_ version of the `leiden` R package. This version
has remote dependencies on the development version of the R 
[igraph](https://github.com/igraph/rigraph/tree/dev) package. This must be
cloned and compiled from source. It depends on functions not on CRAN (yet).

This development version is for testing an upcoming release. It is _not_
recommended to use this unless you require features not supported in previous
releases.

### Dependancies

This package requires the 'leidenalg' and 'igraph' modules for python (2) to be installed on your system. For example:

```
pip install leidenalg numpy python-igraph
```

Note you may need to uninstall the igraph 0.1.11 (now deprecated to jgraph) and install python-igraph or igraph-0.7.0:

```
pip uninstall igraph
pip install leidenalg python-igraph
```

The python version can be installed with pip or conda:
  
```
pip uninstall -y igraph
pip install -U -q leidenalg python-igraph
```

```
conda install -c vtraag leidenalg
```

It is also possible to install the python dependencies with reticulate in R.

```{r, eval = FALSE}
library("reticulate")
py_install("python-igraph")
py_install("leidenalg", forge = TRUE)
```

If you do not have root access, you can use `pip install --user` or `pip install --prefix` to install these in your user directory (which you have write permissions for) and ensure that this directory is in your PATH so that Python can find it.

Dependancies can also be installed from a conda repository. This is recommended for Windows users:

```
conda -c vtraag python-igraph leidenalg
```

### Stable release

The stable 'leiden' package and the dependancies can be installed from CRAN:

```R
install.packages("leiden")
```

### Development version

The 'devtools' package can also be used to install development version of 'leiden' and the dependancies (igraph and reticulate) from GitHub:

```R
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("TomKellyGenetics/leiden", ref = "master")
```

### Development version

To use or test the development version, install the "dev" branch from GitHub.

```R
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("TomKellyGenetics/leiden", ref = "dev")
```

Please submit pull requests to the "dev" branch. This can be downloaded to your system with:

```
git clone --branch dev git@github.com:TomKellyGenetics/leiden.git
```

## Usage

This package provides a function to perform clustering with the Leiden algorithm:

```R
partition <- leiden(adjacency_matrix)
```

### Use with iGraph

For an igraph object 'graph' in R:

```R
adjacency_matrix <- igraph::as_adjacency_matrix(graph)
partition <- leiden(adjacency_matrix)
```

Calling leiden directly on a graph object is also available:

```R
partition <- leiden(graph_object)
```

See the benchmarking vignette on details of performance.

### Computing partitions on data matrices or dimension reductions

To generate an adjacency matrix from a dataset, we can compute the shared nearest neighbours (SNN) from the data. For example, for a dataset `data_mat` with `n` features (rows) by `m` samples or cells (columns), we generate an adjacency matrix of nearest neighbours between samples.

```R
library(RANN)
snn <- RANN::nn2(t(data_mat), k=30)$nn.idx
adjacency_matrix <- matrix(0L, ncol(data_mat), ncol(data_mat))
rownames(adjacency_matrix) <- colnames(adjacency_matrix) <- colnames(data_mat)
for(ii in 1:ncol(data_mat)) {
    adjacency_matrix[i,colnames(data_mat)[snn[ii,]]] <- 1L
}
#check that rows add to k
sum(adjacency_matrix[1,]) == 30
table(apply(adjacency_matrix, 1, sum))
```

For a dimension reduction `embedding` of `m` samples (rows) by `n` dimensions (columns):

```R
library(RANN)
snn <- RANN::nn2(embedding, k=30)$nn.idx
adjacency_matrix <- matrix(0L, nrow(embedding), nrow(embedding))
rownames(adjacency_matrix) <- colnames(adjacency_matrix) <- colnames(data_mat)
for(ii in 1:nrow(embedding)) {
    adjacency_matrix[ii,rownames(data_mat)[snn[ii,]]] <- 1L
}
#check that rows add to k
sum(adjacency_matrix[1,]) == 30
table(apply(adjacency_matrix, 1, sum))
```

This is compatible with PCA, tSNE, or UMAP results.

### Use with Seurat

#### Seurat version 2

To use Leiden with the Seurat pipeline for a Seurat Object `object` that has an SNN computed (for example with `Seurat::FindClusters` with `save.SNN = TRUE`). This will compute the Leiden clusters and add them to the Seurat Object Class.

```R
library("Seurat")
FindClusters(pbmc_small)
adjacency_matrix <- as.matrix(pbmc_small@snn)
partition <- leiden(adjacency_matrix)
pbmc_small@ident <- as.factor(partition)
names(test@ident) <- rownames(test@meta.data)
pbmc_small@meta.data$ident <- as.factor(partition)
```

Seurat objects contain an SNN graph that can be passed directly to the igraph method. For example

```R
library("Seurat")
FindClusters(pbmc_small)
membership <- leiden(pbmc_small@snn)
table(membership)
pbmc_small@ident <- as.factor(membership)
names(pbmc_small@ident) <- rownames(pbmc_small@meta.data)
pbmc_small@meta.data$ident <- as.factor(membership)
```

```R
library("RColorBrewer")
colourPal <- function(groups) colorRampPalette(brewer.pal(min(length(names(table(groups))), 11), "Set3"))(length(names(table(groups))))

pbmc_small <- RunPCA(object = pbmc_small, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = pbmc_small, colors.use = colourPal(pbmc_small@ident), group.by = "ident")

pbmc_small <- RunTSNE(object = pbmc_small, dims.use = 1:20, do.fast = TRUE, dim.embed = 2)
TSNEPlot(object = pbmc_small, colors.use = colourPal(pbmc_small@ident), group.by = "ident")

pbmc_small <- RunUMAP(object = pbmc_small, dims.use = 1:20, metric = "correlation", max.dim = 2)
DimPlot(pbmc_small, reduction.use = "umap", colors.use = colourPal(pbmc_small@ident), group.by = "ident")
```

#### Seurat version 3 (or higher)

Note that this code is designed for Seurat version 2 releases. For Seurat version 3 objects, the Leiden algorithm will be implemented in the Seurat version 3 package with `Seurat::FindClusters` and `algorithm = "leiden"`).  

```R
library("Seurat")
FindClusters(pbmc_small, algorithm = 4)
```

These clusters can then be plotted with:

```R
library("RColorBrewer")
colourPal <- function(groups) colorRampPalette(brewer.pal(min(length(names(table(groups))), 11), "Set3"))(length(names(table(groups))))

PCAPlot(object = pbmc_small, colors.use = colourPal(pbmc_small@active.ident), group.by = "ident")

TSNEPlot(object = pbmc_small, colors.use = colourPal(pbmc_small@active.ident), group.by = "ident")

pbmc_small <- RunUMAP(object = pbmc_small, reduction.use = "pca", dims.use = 1:20, metric = "correlation", max.dim = 2)
DimPlot(pbmc_small, reduction.use = "umap", colors.use = colourPal(pbmc_small@active.ident), group.by = "ident")
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

For more details see the follow vignettes:

* running leiden on an adjacency matrix

[https://github.com/TomKellyGenetics/leiden/blob/master/vignettes/run_leiden.html](https://github.com/TomKellyGenetics/leiden/blob/master/vignettes/run_leiden.html)

* running leiden on an igraph object

[https://github.com/TomKellyGenetics/leiden/blob/master/vignettes/run_igraph.html](https://github.com/TomKellyGenetics/leiden/blob/master/vignettes/run_igraph.html)

* comparing running leiden in Python to various methods in R

[https://github.com/TomKellyGenetics/leiden/blob/master/vignettes/benchmarking.html](https://github.com/TomKellyGenetics/leiden/blob/master/vignettes/run_benchmarking.html)


### Citation

Please cite this implementation R in if you use it:

```
To cite the leiden package in publications use:

  S. Thomas Kelly (2020). leiden: R implementation of the Leiden algorithm. R
  package version 0.4.0 https://github.com/TomKellyGenetics/leiden

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {leiden: R implementation of the Leiden algorithm},
    author = {S. Thomas Kelly},
    year = {2020},
    note = {R package version 0.4.0},
    url = {https://github.com/TomKellyGenetics/leiden},
  }
 ```

Please also cite the original publication of this algorithm.

```
Traag, V.A., Waltman. L., Van Eck, N.-J. (2019). From Louvain to
       Leiden: guaranteeing well-connected communities.
       Sci Rep 9, 5233 <https://doi.org/10.1038/s41598-019-41695-z>
```
