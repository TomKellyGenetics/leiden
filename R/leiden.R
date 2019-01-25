##' @name Leiden State Matrix
##' @rdname leiden
##'
##' @title Run Leiden clustering algorithm
##'
##' @description Implements the Leiden clustering algorithm in R using reticulate to run the Python version. Requires the python "leidenalg" and "igraph" modules to be installed. Returns a vector of partition indices.
##'
##' @param adj_mat An adjacency matrix compatible with \code{\link[igraph]{igraph}} object.
##' @param ... Parameters to pass to the Python leidenalg function.
##'
##' @keywords graph network igraph mvtnorm simulation
##' @import reticulate
##' @export

leiden <- function(adj_mat, ...){
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
    part <- leidenalg$find_partition(snn_graph, leidenalg$ModularityVertexPartition, ...)
    return(part$membership+1)
}
