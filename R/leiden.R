##' Run Leiden clustering algorithm
##'
##' @description Implements the Leiden clustering algorithm in R using reticulate to run the Python version. Requires the python "leidenalg" and "igraph" modules to be installed. Returns a vector of partition indices.
##' @param adj_mat An adjacency matrix compatible with \code{\link[igraph]{igraph}} object.
##' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition (see the Leiden python module documentation for more details)
##' @param initial_membership,weights,node_sizes Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None).
##' @param resolution_parameter A parameter controlling the coarseness of the clusters
##' @return A parition of clusters as a vector of integers
##' @examples
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
#' \dontrun{
#' library("igraph")
#' graph_object <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
#' }
##' #plot graph structure
#' \dontrun{
#' library("devtools")
#' install_github("TomKellyGenetics/igraph.extensions")
#' library("plot.igraph")
#' plot_directed(graph_object, cex.arrow = 0.3, col.arrow = "grey50")
#' }
##' #generate partitions
##' partition <- leiden(adjacency_matrix)
##' table(partition)
##' #plot results
#' \dontrun{
#' library("RColorBrewer")
#' node.cols <- brewer.pal(max(partition),"Pastel1")[partition]
#' plot_directed(graph_object, cex.arrow = 0.3, col.arrow = "grey50", fill.node = node.cols)
#' }
##'
##' @keywords graph network igraph mvtnorm simulation
##' @importFrom reticulate import r_to_py
##' @export

leiden <- function(adj_mat,
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

    #convert matrix input
    adj_mat <- as.matrix(ceiling(adj_mat))

    ##convert to python numpy.ndarray, then a list
    adj_mat_py <- r_to_py(adj_mat)
    adj_mat_py <- adj_mat_py$tolist()

    #convert graph structure to a Python compatible object
    snn_graph <- ig$Graph$Adjacency(adj_mat_py)

    #compute partitions
    partition_type <- partition_type[1]
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
    return(part$membership+1)
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
