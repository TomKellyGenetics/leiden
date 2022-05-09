# leiden 0.4.2

* migrates bug fixes to conda environment and limit on total cells to refactored package

* allows calling igraph::community_leiden for supported parameters (requires igraph v1.2.7 or later)

* automatically calls native R version of leiden rather than Python to improve performance

* updates vignettes and unit tests to ensure consistent results with past versions

# leiden 0.4.1

* migrates changes to retain on CRAN to leiden 0.4.0 (alpha)

# leiden 0.4.0

* migrate to calling community_leiden in igraph

* updates the benchmarking vignette to compare performance to legacy versions using reticulate

# leiden 0.3.10

* removes limitation on number of cells (disables scientific notation within function call): resolves #12

* resolves conflict between base and r-reticulate conda environments on loading: resolves #20

* updates conda environment in interactive sessions only for compliance to CRAN checks

* resolves formatting error in Rmarkdown vignettes (https://github.com/yihui/knitr/issues/2057)

* update testing for bipartite graphs for compatibility with newer version

# leiden 0.3.9

Updates maintainer contact details.

# leiden 0.3.8

* bug fixes for vignettes to retain on CRAN

# leiden 0.3.7

* remove cairo graphics from vignette to retain on CRAN

# leiden 0.3.6

* add methods for multiplex community detection from a list of graphs (requires leidenalg 0.7.1)

* add support for maximum community size (depends on leidenalg 0.8.2), if available

# leiden 0.3.5

* background changes to build vignettes on CRAN

# leiden 0.3.4

* add support for bipartite graphs (requires leidenalg 0.6.1 or later)

* changes to install python leidenalg from vtraag channel on Windows

* bug fixes to install documentation

* improve automated conda configuration in background on loading library

* optionally derive edge weights from the Laplacian matrix

# leiden 0.3.3

* bug fixes for documentation in response to changes to R

see development version: https://bugs.r-project.org/bugzilla/show_bug.cgi?id=16223

# leiden 0.3.2

* added support for passing weighted igraph objects

* improved handling of sparse matrices

* bug fixes to ensure same results from matrix and igraph methods

# leiden 0.3.1

* method for sparse matrices that passes to igraph without casting to dense

* added seed and n_iterations to find_partitions

# leiden 0.3.0

* Implements calling leiden directly on an igraph object

* Separate methods for igraph objects and adjacency matrices

* Support for sparse matrices

* Benchmarking added to vignettes

# leiden 0.2.3

* Adds passing weighted adjacency matrices to derive weight parameters

# leiden 0.2.2

* Changes to ensure compatibility with CRAN. Updates to vignettes and documentation.

# leiden 0.2.1

* Enable passing arguments to Python functions: initial_membership, weights, and
node_sizes.

# leiden 0.2.0

* Adds passing arguments to the Python implementation: the partition_type and resolution_parameter. Runs the RBConfigurationVertexPartition by default (which is equivalent to ModularityVertexPartition with a resolution_parameter of 1).

# leiden 0.1.1

* Removes dependancy on igraph R package and avoids writing to disk (compatible with CRAN). Passes adjacency matrix directly to python as a numpy array.

# leiden 0.1.0

* Implements the Leiden algorithm in R by calling the leidenalg python library. Runs ModularityVertexPartition with defaults.
