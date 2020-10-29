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
