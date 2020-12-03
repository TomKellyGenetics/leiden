##' internal function to compute bipartite paritions
##' @param py_graph an igraph object in python (reticulate)
##' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition, ModularityVertexPartition.Bipartite, CPMVertexPartition.Bipartite (see the Leiden python module documentation for more details)
##' @param initial_membership,weights,node_sizes Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None). Weights are derived from weighted igraph objects and non-zero integer values of adjacency matrices.
##' @param resolution_parameter_1,resolution_parameter_0,resolution_parameter_01 A parameter controlling the coarseness of the clusters
##' @param seed Seed for the random number generator. By default uses a random seed if nothing is specified.
##' @param n_iterations Number of iterations to run the Leiden algorithm. By default, 2 iterations are run. If the number of iterations is negative, the Leiden algorithm is run until an iteration in which there was no improvement.
##' @param max_comm_size (non-negative int) – Maximal total size of nodes in a community. If zero (the default), then communities can be of any size.
##' @param degree_as_node_size (defaults to FALSE). If True use degree as node size instead of 1, to mimic modularity for Bipartite graphs.
##' @param types (deftaults to "type"). Defines edge attribute for bipartite graphs.
##' @noRd
##' @description internal function to compute partitions for bipartite graphs
##' @keywords internal
run_bipartite_partitioning <- function(py_graph,
                                       initial_membership = NULL,
                                       weights = NULL,
                                       resolution_parameter_01 = 1,
                                       resolution_parameter_0 = 0,
                                       resolution_parameter_1 = 0,
                                       degree_as_node_size = FALSE,
                                       types = "type",
                                       seed = NULL,
                                       n_iterations = 2L,
                                       max_comm_size = 0L){
  #import python modules with reticulate
  numpy <- import("numpy", delay_load = TRUE)
  leidenalg <- import("leidenalg", delay_load = TRUE)
  ig <- import("igraph", delay_load = TRUE)

  self.optimiser = leidenalg$Optimiser()
  if(py_has_attr(self.optimiser, "max_comm_size")){
    py_set_attr(self.optimiser, "max_comm_size", r_to_py(as.integer(max_comm_size)))
  }

  if(!is.null(seed)){
    self.optimiser$set_rng_seed(r_to_py(as.integer(seed)))
  }
  for(ii in 1:n_iterations){
    bipartite_layers <- leidenalg$CPMVertexPartition$Bipartite(py_graph,
                                                               initial_membership = initial_membership,
                                                               weights = weights,
                                                               resolution_parameter_01 = resolution_parameter_01,
                                                               resolution_parameter_0 = resolution_parameter_0,
                                                               resolution_parameter_1 = resolution_parameter_1,
                                                               degree_as_node_size = degree_as_node_size,
                                                               types = r_to_py(types))
  }
  bipartite_layers <- r_to_py(bipartite_layers)
  self.optimiser$optimise_partition_multiplex(
    bipartite_layers,
    layer_weights=r_to_py(c(1L, -1L, -1L)))
  part <- py_to_r(bipartite_layers[[1]])
  part
}

##' call leiden on a python graph object with reticulate
##' @param py_graph an igraph object
##' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition, ModularityVertexPartition.Bipartite, CPMVertexPartition.Bipartite (see the Leiden python module documentation for more details)
##' @param initial_membership,weights,node_sizes Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None). Weights are derived from weighted igraph objects and non-zero integer values of adjacency matrices.
##' @param resolution_parameter A parameter controlling the coarseness of the clusters
##' @param seed Seed for the random number generator. By default uses a random seed if nothing is specified.
##' @param n_iterations Number of iterations to run the Leiden algorithm. By default, 2 iterations are run. If the number of iterations is negative, the Leiden algorithm is run until an iteration in which there was no improvement.
##' @param max_comm_size (non-negative int) – Maximal total size of nodes in a community. If zero (the default), then communities can be of any size.
##' @param degree_as_node_size (defaults to FALSE). If True use degree as node size instead of 1, to mimic modularity for Bipartite graphs.
##' @noRd
##' @description internal function to compute partitions by calling Python with reticulate
##' @keywords internal
find_partition <- function(py_graph, partition_type = c(
  'RBConfigurationVertexPartition',
  'ModularityVertexPartition',
  'RBERVertexPartition',
  'CPMVertexPartition',
  'MutableVertexPartition',
  'SignificanceVertexPartition',
  'SurpriseVertexPartition',
  'ModularityVertexPartition.Bipartite',
  'CPMVertexPartition.Bipartite'
),
initial_membership = NULL,
weights = NULL,
node_sizes = NULL,
resolution_parameter = 1,
seed = NULL,
n_iterations = 2L,
max_comm_size = 0L,
degree_as_node_size = TRUE,
legacy = FALSE
) {
  #import python modules with reticulate
  numpy <- import("numpy", delay_load = TRUE)
  leidenalg <- import("leidenalg", delay_load = TRUE)
  ig <- import("igraph", delay_load = TRUE)

  partition_type <- match.arg(partition_type)
  if(partition_type == 'ModularityVertexPartition.Bipartite') degree_as_node_size <- TRUE
  if(!is.null(seed)) seed <- as.integer(seed)
  if (is.integer(n_iterations)) n_iterations <- as.integer(n_iterations)

  self.optimiser = leidenalg$Optimiser()
  if(py_has_attr(self.optimiser, "max_comm_size")){
    py_set_attr(self.optimiser, "max_comm_size", r_to_py(as.integer(max_comm_size)))
    part <- switch(
      EXPR = partition_type,
      'RBConfigurationVertexPartition' = leidenalg$find_partition(
        py_graph,
        leidenalg$RBConfigurationVertexPartition,
        initial_membership = initial_membership, weights = weights,
        seed = seed,
        n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size)),
        resolution_parameter = resolution_parameter
      ),
      'ModularityVertexPartition' = leidenalg$find_partition(
        py_graph,
        leidenalg$ModularityVertexPartition,
        initial_membership = initial_membership, weights = weights,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size))
      ),
      'RBERVertexPartition' = leidenalg$find_partition(
        py_graph,
        leidenalg$RBERVertexPartition,
        initial_membership = initial_membership, weights = weights,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size)),
        node_sizes = node_sizes,
        resolution_parameter = resolution_parameter
      ),
      'CPMVertexPartition' = leidenalg$find_partition(
        py_graph,
        leidenalg$CPMVertexPartition,
        initial_membership = initial_membership, weights = weights,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size)),
        node_sizes = node_sizes,
        resolution_parameter = resolution_parameter
      ),
      'MutableVertexPartition' = leidenalg$find_partition(
        py_graph,
        leidenalg$MutableVertexPartition,
        initial_membership = initial_membership,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size))
      ),
      'SignificanceVertexPartition' = leidenalg$find_partition(
        py_graph,
        leidenalg$SignificanceVertexPartition,
        initial_membership = initial_membership,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size)),
        node_sizes = node_sizes,
        resolution_parameter = resolution_parameter
      ),
      'SurpriseVertexPartition' = leidenalg$find_partition(
        py_graph,
        leidenalg$SurpriseVertexPartition,
        initial_membership = initial_membership, weights = weights,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size)),
        node_sizes = node_sizes
      ),
      'ModularityVertexPartition.Bipartite' = run_bipartite_partitioning(py_graph,
                                                                         initial_membership = initial_membership, weights = weights,
                                                                         resolution_parameter_01 = resolution_parameter,
                                                                         resolution_parameter_0 = 0,
                                                                         resolution_parameter_1 = 0,
                                                                         degree_as_node_size = TRUE,
                                                                         types = "type",
                                                                         seed = seed,
                                                                         n_iterations = n_iterations),
      'CPMVertexPartition.Bipartite' = run_bipartite_partitioning(py_graph,
                                                                  initial_membership = initial_membership, weights = weights,
                                                                  resolution_parameter_01 = resolution_parameter,
                                                                  resolution_parameter_0 = 0,
                                                                  resolution_parameter_1 = 0,
                                                                  degree_as_node_size = degree_as_node_size,
                                                                  types = "type",
                                                                  seed = seed,
                                                                  n_iterations = n_iterations),
      stop("please specify a partition type as a string out of those documented")
    )
  } else {

  part <- switch(
    EXPR = partition_type,
    'RBConfigurationVertexPartition' = leidenalg$find_partition(
      py_graph,
      leidenalg$RBConfigurationVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed,
      n_iterations = n_iterations,
      resolution_parameter = resolution_parameter
    ),
    'ModularityVertexPartition' = leidenalg$find_partition(
      py_graph,
      leidenalg$ModularityVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations
    ),
    'RBERVertexPartition' = leidenalg$find_partition(
      py_graph,
      leidenalg$RBERVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'CPMVertexPartition' = leidenalg$find_partition(
      py_graph,
      leidenalg$CPMVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'MutableVertexPartition' = leidenalg$find_partition(
      py_graph,
      leidenalg$MutableVertexPartition,
      initial_membership = initial_membership,
      seed = seed, n_iterations = n_iterations
    ),
    'SignificanceVertexPartition' = leidenalg$find_partition(
      py_graph,
      leidenalg$SignificanceVertexPartition,
      initial_membership = initial_membership,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'SurpriseVertexPartition' = leidenalg$find_partition(
      py_graph,
      leidenalg$SurpriseVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes
    ),
    'ModularityVertexPartition.Bipartite' = run_bipartite_partitioning(py_graph,
                                                                initial_membership = initial_membership, weights = weights,
                                                                resolution_parameter_01 = resolution_parameter,
                                                                resolution_parameter_0 = 0,
                                                                resolution_parameter_1 = 0,
                                                                degree_as_node_size = TRUE,
                                                                types = "type",
                                                                seed = seed,
                                                                n_iterations = n_iterations,
                                                                max_comm_size = max_comm_size),
    'CPMVertexPartition.Bipartite' = run_bipartite_partitioning(py_graph,
                                                                initial_membership = initial_membership, weights = weights,
                                                                resolution_parameter_01 = resolution_parameter,
                                                                resolution_parameter_0 = 0,
                                                                resolution_parameter_1 = 0,
                                                                degree_as_node_size = degree_as_node_size,
                                                                types = "type",
                                                                seed = seed,
                                                                n_iterations = n_iterations,
                                                                max_comm_size = max_comm_size),
    stop("please specify a partition type as a string out of those documented")
  )
}
  partition <- part$membership+1
  partition
}

##' call leiden on a python graph_list object with reticulate
##' @param graph_list list of igraph objects
##' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition, ModularityVertexPartition.Bipartite, CPMVertexPartition.Bipartite (see the Leiden python module documentation for more details)
##' @param initial_membership,weights,node_sizes Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None). Weights are derived from weighted igraph objects and non-zero integer values of adjacency matrices.
##' @param resolution_parameter A parameter controlling the coarseness of the clusters
##' @param seed Seed for the random number generator. By default uses a random seed if nothing is specified.
##' @param n_iterations Number of iterations to run the Leiden algorithm. By default, 2 iterations are run. If the number of iterations is negative, the Leiden algorithm is run until an iteration in which there was no improvement.
##' @param max_comm_size (non-negative int) – Maximal total size of nodes in a community. If zero (the default), then communities can be of any size.
##' @param degree_as_node_size (defaults to FALSE). If True use degree as node size instead of 1, to mimic modularity for Bipartite graphs.
##' @noRd
##' @description internal function to compute partitions by calling Python with reticulate
##' @keywords internal
find_partition_multiplex <- function(graph_list, partition_type = c(
  'RBConfigurationVertexPartition',
  'ModularityVertexPartition',
  'RBERVertexPartition',
  'CPMVertexPartition',
  'MutableVertexPartition',
  'SignificanceVertexPartition',
  'SurpriseVertexPartition'
),
initial_membership = NULL,
weights = NULL,
node_sizes = NULL,
resolution_parameter = 1,
seed = NULL,
n_iterations = 2L,
max_comm_size = 0L,
degree_as_node_size = TRUE,
legacy = FALSE
) {
  #import python modules with reticulate
  numpy <- import("numpy", delay_load = TRUE)
  leidenalg <- import("leidenalg", delay_load = TRUE)
  ig <- import("igraph", delay_load = TRUE)

  partition_type <- match.arg(partition_type)
  if(!is.null(seed)) seed <- as.integer(seed)
  if (is.integer(n_iterations)) n_iterations <- as.integer(n_iterations)

  self.optimiser = leidenalg$Optimiser()
  if(py_has_attr(self.optimiser, "max_comm_size")){
    py_set_attr(self.optimiser, "max_comm_size", r_to_py(as.integer(max_comm_size)))
    part <- switch(
      EXPR = partition_type,
      'RBConfigurationVertexPartition' = leidenalg$find_partition_multiplex(
        graph_list,
        leidenalg$RBConfigurationVertexPartition,
        initial_membership = initial_membership, weights = weights,
        seed = seed,
        n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size)),
        resolution_parameter = resolution_parameter
      ),
      'ModularityVertexPartition' = leidenalg$find_partition_multiplex(
        graph_list,
        leidenalg$ModularityVertexPartition,
        initial_membership = initial_membership, weights = weights,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size))
      ),
      'RBERVertexPartition' = leidenalg$find_partition_multiplex(
        graph_list,
        leidenalg$RBERVertexPartition,
        initial_membership = initial_membership, weights = weights,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size)),
        node_sizes = node_sizes,
        resolution_parameter = resolution_parameter
      ),
      'CPMVertexPartition' = leidenalg$find_partition_multiplex(
        graph_list,
        leidenalg$CPMVertexPartition,
        initial_membership = initial_membership, weights = weights,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size)),
        node_sizes = node_sizes,
        resolution_parameter = resolution_parameter
      ),
      'MutableVertexPartition' = leidenalg$find_partition_multiplex(
        graph_list,
        leidenalg$MutableVertexPartition,
        initial_membership = initial_membership,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size))
      ),
      'SignificanceVertexPartition' = leidenalg$find_partition_multiplex(
        graph_list,
        leidenalg$SignificanceVertexPartition,
        initial_membership = initial_membership,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size)),
        node_sizes = node_sizes,
        resolution_parameter = resolution_parameter
      ),
      'SurpriseVertexPartition' = leidenalg$find_partition_multiplex(
        graph_list,
        leidenalg$SurpriseVertexPartition,
        initial_membership = initial_membership, weights = weights,
        seed = seed, n_iterations = n_iterations,
        max_comm_size = r_to_py(as.integer(max_comm_size)),
        node_sizes = node_sizes
      ),
      stop("please specify a partition type as a string out of those documented")
    )
  } else {
  part <- switch(
    EXPR = partition_type,
    'RBConfigurationVertexPartition' = leidenalg$find_partition_multiplex(
      graph_list,
      leidenalg$RBConfigurationVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed,
      n_iterations = n_iterations,
      resolution_parameter = resolution_parameter
    ),
    'ModularityVertexPartition' = leidenalg$find_partition_multiplex(
      graph_list,
      leidenalg$ModularityVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations
    ),
    'RBERVertexPartition' = leidenalg$find_partition_multiplex(
      graph_list,
      leidenalg$RBERVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'CPMVertexPartition' = leidenalg$find_partition_multiplex(
      graph_list,
      leidenalg$CPMVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'MutableVertexPartition' = leidenalg$find_partition_multiplex(
      graph_list,
      leidenalg$MutableVertexPartition,
      initial_membership = initial_membership,
      seed = seed, n_iterations = n_iterations
    ),
    'SignificanceVertexPartition' = leidenalg$find_partition_multiplex(
      graph_list,
      leidenalg$SignificanceVertexPartition,
      initial_membership = initial_membership,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'SurpriseVertexPartition' = leidenalg$find_partition_multiplex(
      graph_list,
      leidenalg$SurpriseVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes
    ),
    stop("please specify a partition type as a string out of those documented")
  )
}
  names(part) <- c("membership", "improv")
  partition <- part$membership+1
  partition
}
