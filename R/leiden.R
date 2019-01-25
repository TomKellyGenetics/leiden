##' @name Leiden State Matrix
##' @rdname leiden
##'
##' @title Run Leiden clustering algorithm
##'
##' @description Implements the Leiden clustering algorithm in R using reticulate to run the Python version. Requires the python "leidenalg" and "igraph" modules to be installed. Returns a vector of partition indices.
##'
##' @param adj_mat An adjacency matrix compatible with \code{\link[igraph]{igraph}} object.
##' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition (see the Leiden python module documentation for more details)
##' @param resolution_parameter A parameter controlling the coarseness of the clusters. Higher values lead to more clusters. (defaults to 1.0 for partition types that accept a resolution parameter)
##' @param ... Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None).
##'
##' @keywords graph network igraph mvtnorm simulation
##' @import reticulate
##' @export

leiden <- function(adj_mat, partition_type = c('RBConfigurationVertexPartition', 'ModularityVertexPartition', 'RBERVertexPartition', 'CPMVertexPartition', 'MutableVertexPartition', 'SignificanceVertexPartition', 'SurpriseVertexPartition'), resolution_parameter = 1, ...){
    #import python modules with reticulate
    leidenalg <- reticulate::import("leidenalg")
    ig <- reticulate::import("igraph")

    #convert matrix input
    adj_mat <- as.matrix(ceiling(adj_mat))

    ##convert to python numpy.ndarray, then a list
    adj_mat_py <- reticulate::r_to_py(adj_mat)
    adj_mat_py <- adj_mat_py$tolist()

    #convert graph structure to a Python compatible object
    snn_graph <- ig$Graph$Adjacency(adj_mat_py)

    #compute partitions
    partition_type <- partition_type[1]
    if(partition_type == "RBConfigurationVertexPartition"){
        part <- leidenalg$find_partition(snn_graph, leidenalg$RBConfigurationVertexPartition, resolution_parameter = resolution_parameter, ...)
    } else if(partition_type == "ModularityVertexPartition"){
        part <- leidenalg$find_partition(snn_graph, leidenalg$ModularityVertexPartition, ...)
    } else if(partition_type == "RBERVertexPartition"){
        part <- leidenalg$find_partition(snn_graph, leidenalg$RBERVertexPartition, resolution_parameter = resolution_parameter, ...)
    } else if(partition_type == "CPMVertexPartition"){
        part <- leidenalg$find_partition(snn_graph, leidenalg$CPMVertexPartition, resolution_parameter = resolution_parameter, ...)
    } else if(partition_type == "MutableVertexPartition"){
        part <- leidenalg$find_partition(snn_graph, leidenalg$MutableVertexPartition, ...)
    } else if(partition_type == "SignificanceVertexPartition"){
        part <- leidenalg$find_partition(snn_graph, leidenalg$SignificanceVertexPartition, resolution_parameter = resolution_parameter, ...)
    } else if(partition_type == "SurpriseVertexPartition"){
        part <- leidenalg$find_partition(snn_graph, leidenalg$SurpriseVertexPartition, resolution_parameter = resolution_parameter, ...)
    } else {
        stop("please specify a partition type as a string out of those documented")
    }
    return(part$membership+1)
}
