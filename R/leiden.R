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
##' @importFrom igraph graph_from_adjacency_matrix write.graph
##' @export

leiden <- function(adj_mat, ...){
    #import python modules with reticulate
    leidenalg <- reticulate::import("leidenalg")
    ig <- reticulate::import("igraph")

    #export graph structure to call from Python
    snn_graph <- igraph::graph_from_adjacency_matrix(ceiling(as.matrix(adj_mat)))
    igraph::write.graph(snn_graph, file = ".graph.edges", format = "graphml")

    #import graph structure as a Python compatible object
    snn_graph <-  ig$Graph$Read_GraphML(".graph.edges")

    #compute partitions
    part <- leidenalg$find_partition(snn_graph, leidenalg$ModularityVertexPartition, ...)
    return(part$membership+1)
}
