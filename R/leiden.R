##' Run Leiden clustering algorithm
##'
##' @description Implements the Leiden clustering algorithm in R using reticulate to run the Python version. Requires the python "leidenalg" and "igraph" modules to be installed. Returns a vector of partition indices.
##' @param adj_mat An adjacency matrix compatible with \code{\link[igraph]{igraph}} object.
##' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition (see the Leiden python module documentation for more details)
##' @param initial_membership,weights,node_sizes Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None).
##' @param resolution_parameter A parameter controlling the coarseness of the clusters
##' @return A parition of clusters as a vector of integers
##' @examples
##' #check if python is availble
##' modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
##' if(modules){
##' #generate example data
##' adjacency_matrix <- rbind(cbind(matrix(round(rbinom(4000, 1, 0.8)), 20, 20),
##'                                 matrix(round(rbinom(4000, 1, 0.3)), 20, 20),
##'                                 matrix(round(rbinom(400, 1, 0.1)), 20, 20)),
##'                           cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20),
##'                                 matrix(round(rbinom(400, 1, 0.8)), 20, 20),
##'                                 matrix(round(rbinom(4000, 1, 0.2)), 20, 20)),
##'                           cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20),
##'                                 matrix(round(rbinom(4000, 1, 0.1)), 20, 20),
##'                                 matrix(round(rbinom(4000, 1, 0.9)), 20, 20)))
##' rownames(adjacency_matrix) <- 1:60
##' colnames(adjacency_matrix) <- 1:60
##' #generate partitions
##' partition <- leiden(adjacency_matrix)
##' table(partition)
##'
##' #generate partitions at a lower resolution
##' partition <- leiden(adjacency_matrix, resolution_parameter = 0.5)
##' table(partition)
##'
##' #generate example weights
##' weights <- sample(1:10, sum(adjacency_matrix!=0), replace=TRUE)
##' partition <- leiden(adjacency_matrix, weights = weights)
##' table(partition)
##'
##' #generate example weighted matrix
##' adjacency_matrix[adjacency_matrix == 1] <- weights
##' partition <- leiden(adjacency_matrix)
##' table(partition)
##' }
##'
##'
##' @keywords graph network igraph mvtnorm simulation
##' @importFrom reticulate import r_to_py
##' @export
##'
#' @export
#' @usage NULL
leiden <- function(x, ...) {
    UseMethod("leiden")
}

leiden.data.frame <- function(...) {
    leiden.matrix(...)
}

leiden.Matrix <- function(...) {
    leiden.matrix(...)
}


leiden.matrix <- function(adj_mat,
                          partition_type = c(
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
    #import python modules with reticulate
    leidenalg <- import("leidenalg", delay_load = TRUE)
    ig <- import("igraph", delay_load = TRUE)

    #convert matrix input (corrects for sparse matrix input)
    adj_mat <- as.matrix(adj_mat)

    #compute weights if non-binary adjacency matrix given
    is_pure_adj <- all(as.logical(adj_mat) == adj_mat)
    if (is.null(weights) && !is_pure_adj) {
        #assign weights to edges (without dependancy on igraph)
        weights <- t(adj_mat)[t(adj_mat)!=0]
        #remove zeroes from rows of matrix and return vector of length edges
    }

    ##convert to python numpy.ndarray, then a list
    adj_mat_py <- r_to_py(adj_mat)
    adj_mat_py <- adj_mat_py$tolist()

    #convert graph structure to a Python compatible object
    GraphClass <- if (!is.null(weights) && !is_pure_adj){
        ig$Graph$Weighted_Adjacency
    } else {
        ig$Graph$Adjacency
    }
    snn_graph <- GraphClass(adj_mat_py)

    #compute partitions
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

leiden.igraph <- function(snn_graph,
                          partition_type = c(
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
    #import python modules with reticulate
    leidenalg <- import("leidenalg", delay_load = TRUE)
    ig <- import("igraph", delay_load = TRUE)

    ##convert to python numpy.ndarray, then a list
    vertices <- as.list(names(V(snn_graph)))
    edges <- as_edgelist(snn_graph)
    dim(edges)
    edgelist <- list(rep(NA, nrow(edges)))
    for(ii in 1:nrow(edges)){
        edgelist[[ii]] <- edges[ii,]
    }

    snn_graph_py <- ig$Graph()
    snn_graph_py$add_vertices(r_to_py(vertices))
    snn_graph_py$add_edges(r_to_py(edgelist))

    #compute weights if weighted graph given
    if (is.weighted(snn_graph)) {
        #assign weights to edges (without dependancy on igraph)
        weights <- r_to_py(weights(snn_graph))
        snn_graph_py$es['weights'] <- weights
    }

    snn_graph <- snn_graph_py

    #compute partitions
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


# global reference to python modules (will be initialized in .onLoad)
leidenalg <- NULL
ig <- NULL

.onLoad = function(libname, pkgname) {
    if(reticulate::py_available()){
        install_python_modules <- function(method = "auto", conda = "auto") {
            reticulate::py_install("leidenalg", method = method, conda = conda)
            reticulate::py_install("igraph", method = method, conda = conda)
        }
    }
    if (suppressWarnings(suppressMessages(requireNamespace("reticulate")))) {
        modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
        if (modules) {
            ## assignment in parent environment!
            leidenalg <<- reticulate::import("leidenalg", delay_load = TRUE)
            ig <<- reticulate::import("igraph", delay_load = TRUE)
        }
    }
}
