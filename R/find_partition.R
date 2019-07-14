#call leiden on a python snn_graph object with reticulate
find_partition <- function(snn_graph, partition_type = c(
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
resolution_parameter = 1
) {
  partition_type <- match.arg(partition_type)
  part <- switch(
    EXPR = partition_type,
    'RBConfigurationVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBConfigurationVertexPartition,
      initial_membership = initial_membership, weights = weights,
      resolution_parameter = resolution_parameter
    ),
    'ModularityVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$ModularityVertexPartition,
      initial_membership = initial_membership, weights = weights
    ),
    'RBERVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$RBERVertexPartition,
      initial_membership = initial_membership, weights = weights, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'CPMVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$CPMVertexPartition,
      initial_membership = initial_membership, weights = weights, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'MutableVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$MutableVertexPartition,
      initial_membership = initial_membership
    ),
    'SignificanceVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SignificanceVertexPartition,
      initial_membership = initial_membership, node_sizes = node_sizes,
      resolution_parameter = resolution_parameter
    ),
    'SurpriseVertexPartition' = leidenalg$find_partition(
      snn_graph,
      leidenalg$SurpriseVertexPartition,
      initial_membership = initial_membership, weights = weights, node_sizes = node_sizes
    ),
    stop("please specify a partition type as a string out of those documented")
  )
  partition <- part$membership+1
  partition
}
