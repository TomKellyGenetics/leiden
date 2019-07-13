##' Run Leiden clustering algorithm
##'
##' @description Implements the Leiden clustering algorithm in R using reticulate to run the Python version. Requires the python "leidenalg" and "igraph" modules to be installed. Returns a vector of partition indices.
##' @param object An adjacency matrix compatible with \code{\link[igraph]{igraph}} object or an input graph as an \code{\link[igraph]{igraph}} object (e.g., shared nearest neighbours).
##' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition (see the Leiden python module documentation for more details)
##' @param initial_membership,weights,node_sizes Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None).
##' @param resolution_parameter A parameter controlling the coarseness of the clusters
##' @param ... 	Arguments to be passed to methods
##' @return A partition of clusters as a vector of integers
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
##'
##'
##' # generate (unweighted) igraph object in R
##' library("igraph")
##' adjacency_matrix[adjacency_matrix > 1] <- 1
##' snn_graph <- graph_from_adjacency_matrix(adjacency_matrix)
##' partition <- leiden(snn_graph)
##' table(partition)
##'
##' # generate (weighted) igraph object in R
##' library("igraph")
##' adjacency_matrix[adjacency_matrix >= 1] <- weights
##' snn_graph <- graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE)
##' partition <- leiden(snn_graph)
##' table(partition)
##'
##' # pass weights to python leidenalg
##' adjacency_matrix[adjacency_matrix >= 1 ] <- 1
##' snn_graph <- graph_from_adjacency_matrix(adjacency_matrix, weighted = NULL)
##' weights <- sample(1:10, sum(adjacency_matrix!=0), replace=TRUE)
##' partition <- leiden(snn_graph, weights = weights)
##' table(partition)
##'
##' # run only if python is available (for testing)
##' }
##'
##' @keywords graph network igraph mvtnorm simulation
##' @importFrom reticulate import r_to_py
##' @export
##'
#' @export
#' @rdname leiden
leiden <- function(object,
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
                   resolution_parameter = 1) {
    UseMethod("leiden", object)
}

#' @export
leiden.matrix <- function(object,
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
    if(is.matrix(object) || is(object, "Matrix")){
        adj_mat <- object
    } else{
        adj_mat <- as.matrix(object)
    }

    #compute weights if non-binary adjacency matrix given
    is_pure_adj <- all(as.logical(adj_mat) == adj_mat)
    if (is.null(weights) && !is_pure_adj) {
        if(!is.matrix(object)) adj_mat <- as.matrix(adj_mat)
        #assign weights to edges (without dependancy on igraph)
        t_mat <- t(adj_mat)
        weights <- t_mat[t_mat!=0]
        #remove zeroes from rows of matrix and return vector of length edges
    }

    ##convert to python numpy.ndarray, then a list
    adj_mat_py <- r_to_py(adj_mat)
    if(is(object, "Matrix")){
        adj_mat_py <- adj_mat_py$todense()
    }
    adj_mat_py <- adj_mat_py$tolist()

    #convert graph structure to a Python compatible object
    GraphClass <- if (!is.null(weights) && !is_pure_adj){
        ig$Graph$Weighted_Adjacency
    } else {
        ig$Graph$Adjacency
    }
    snn_graph <- GraphClass(adj_mat_py)
    # test performance for computing matrix to graph in R
    # other option is to passing snn_graph to Python

    #compute partitions
    #compute partitions
    partition <- find_partition(snn_graph, partition_type = partition_type,
                                initial_membership = initial_membership ,
                                weights = weights,
                                node_sizes = node_sizes,
                                resolution_parameter = resolution_parameter
    )
    partition
}

#' @export
leiden.data.frame <- leiden.matrix
#' @export
#' @importClassesFrom Matrix dgCMatrix
leiden.Matrix <- leiden.matrix
#' @export
leiden.default <- leiden.matrix

#' @importFrom igraph V as_edgelist is.weighted is.named
#' @export
leiden.igraph <- function(object,
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
    if(!is.named(object)){
        vertices <- as.list(as.character(V(object)))
    } else {
        vertices <- as.list(names(V(object)))
    }

    edges <- as_edgelist(object)
    dim(edges)
    edgelist <- list(rep(NA, nrow(edges)))
    for(ii in 1:nrow(edges)){
        edgelist[[ii]] <- as.character(edges[ii,])
    }

    snn_graph <- ig$Graph()
    snn_graph$add_vertices(r_to_py(vertices))
    snn_graph$add_edges(r_to_py(edgelist))

    #compute weights if weighted graph given
    if (is.weighted(object)) {
        #assign weights to edges (without dependancy on igraph)
        weights <- r_to_py(weights(object))
        snn_graph$es$set_attribute_values('weight', weights)
    }

    # from here is the same as method for matrix
    # would be better to refactor to call from matrix methof

    #compute partitions
    partition <- find_partition(snn_graph, partition_type = partition_type,
                                initial_membership = initial_membership ,
                                weights = weights,
                                node_sizes = node_sizes,
                                resolution_parameter = resolution_parameter
    )
    partition
}

leiden.default <- leiden.matrix

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
