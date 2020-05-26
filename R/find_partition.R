#call leiden on a python snn_graph object with reticulate
find_partition <- function(snn_graph, partition_type = c(
  'RBConfigurationVertexPartition',
  'ModularityVertexPartition',
  'RBERVertexPartition',
  'CPMVertexPartition',
  'MutableVertexPartition',
  'SignificanceVertexPartition',
  'SurpriseVertexPartition',
  'CPMVertexPartition.Bipartite'
),
initial_membership = NULL,
weights = NULL,
node_sizes = NULL,
resolution_parameter = 1,
seed = NULL,
n_iterations = 2L
) {
  partition_type <- match.arg(partition_type)
  if(!is.null(seed)) seed <- as.integer(seed)
  if (is.integer(n_iterations)) n_iterations <- as.integer(n_iterations)
  part <- switch(
    EXPR = partition_type,
    'RBConfigurationVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBConfigurationVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed,
      n_iterations = n_iterations,
      resolution_parameter = resolution_parameter
    ),
    'ModularityVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$ModularityVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations
    ),
    'RBERVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBERVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'CPMVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$CPMVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'MutableVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$MutableVertexPartition,
      initial_membership = initial_membership,
      seed = seed, n_iterations = n_iterations
    ),
    'SignificanceVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SignificanceVertexPartition,
      initial_membership = initial_membership,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'SurpriseVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SurpriseVertexPartition,
      initial_membership = initial_membership, weights = weights,
      seed = seed, n_iterations = n_iterations, node_sizes = node_sizes
    ),
    'CPMVertexPartition.Bipartite' =
      optimiser <- leidenalg$Optimiser()
      bipartite_layers <- leidenalg$CPMVertexPartition$Bipartite(snn_graph,
                                             resolution_parameter_01 = 0.1,
                                             resolution_parameter_0 = 0,
                                             resolution_parameter_1 = 0,
                                             degree_as_node_size = FALSE,
                                             types = "type")
      optimiser$optimise_partition_multiplex(r_to_py(bipartite_partition), layer_weights = r_to_py(c(1L, -1L, -1L)))
   ,
    stop("please specify a partition type as a string out of those documented")
  )
  partition <- part$membership+1
  partition
}
