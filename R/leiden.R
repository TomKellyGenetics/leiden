##' Run Leiden clustering algorithm
##'
##' @description Implements the Leiden clustering algorithm in R using reticulate to run the Python version. Requires the python "leidenalg" and "igraph" modules to be installed. Returns a vector of partition indices.
##' Windows users can still this with devtools::install_github("rstudio/reticulate", ref = "86ebb56"); reticulate::use_condaenv("r-reticulate"); reticulate::conda_install("r-reticulate", "leidenalg", channel = "vtraag")
##' @param object An adjacency matrix compatible with \code{\link[igraph]{igraph}} object or an input graph as an \code{\link[igraph]{igraph}} object (e.g., shared nearest neighbours).
##' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition, ModularityVertexPartition.Bipartite, CPMVertexPartition.Bipartite (see the Leiden python module documentation for more details)
##' @param initial_membership,weights,node_sizes Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None). Weights are derived from weighted igraph objects and non-zero integer values of adjacency matrices.
##' @param resolution_parameter A parameter controlling the coarseness of the clusters
##' @param seed Seed for the random number generator. By default uses a random seed if nothing is specified.
##' @param n_iterations Number of iterations to run the Leiden algorithm. By default, 2 iterations are run. If the number of iterations is negative, the Leiden algorithm is run until an iteration in which there was no improvement.
##' @param degree_as_node_size (defaults to FALSE). If True use degree as node size instead of 1, to mimic modularity for Bipartite graphs.
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
##' @rdname leiden
##' @export
leiden <- function(object,
                   partition_type = c(
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
                   degree_as_node_size = FALSE,
                   laplacian = FALSE) {
    UseMethod("leiden", object)
}

##' @export
##' @importFrom methods is
leiden.matrix <- function(object,
                          partition_type = c(
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
                          degree_as_node_size = FALSE,
                          laplacian = FALSE
) {
    #import python modules with reticulate
    leidenalg <- import("leidenalg", delay_load = TRUE)
    ig <- import("igraph", delay_load = TRUE)

    #convert matrix input (corrects for sparse matrix input)
    if(is.matrix(object) || is(object, "dgCMatrix")){
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
    adj_mat_py <- r_to_py(adj_mat, convert = T)
    if(is(object, "dgCMatrix")){
        adj_mat_py <- adj_mat_py$toarray()
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
                                resolution_parameter = resolution_parameter,
                                seed = seed,
                                n_iterations = n_iterations,
                                degree_as_node_size = degree_as_node_size
    )
    partition
}

##' @export
leiden.data.frame <- leiden.matrix

##' @importFrom igraph graph_from_adjacency_matrix edge.attributes set.edge.attribute E
##' @importFrom methods as
##' @importClassesFrom Matrix dgCMatrix dgeMatrix
##' @export
leiden.Matrix <- function(object,
                          partition_type = c(
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
                          degree_as_node_size = FALSE,
                          laplacian = FALSE
) {
    #cast to sparse matrix
    adj_mat <- as(object, "dgCMatrix")
    #run as igraph object (passes to reticulate)
    if(is.null(weights)){
        object <- graph_from_adjacency_matrix(adjmatrix = adj_mat, weighted = TRUE)
        weights <- edge.attributes(object)$weight
    } else {
        object <- graph_from_adjacency_matrix(adjmatrix = adj_mat, weighted = TRUE)
        object <- set.edge.attribute(object, "weight", index=E(object), weights)
    }

    leiden.igraph(object,
        partition_type = partition_type,
        weights = weights,
        node_sizes = node_sizes,
        resolution_parameter = resolution_parameter,
        seed = seed,
        n_iterations = n_iterations,
        degree_as_node_size = degree_as_node_size,
        laplacian = laplacian
    )
}

##' @export
leiden.default <- leiden.matrix

##' @importFrom igraph V as_edgelist is.weighted is.named edge.attributes as_adjacency_matrix laplacian_matrix get.vertex.attribute is.bipartite
##' @export
leiden.igraph <- function(object,
                          partition_type = c(
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
                          degree_as_node_size = FALSE,
                          laplacian = FALSE
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

    #derive Laplacian
    if(laplacian == TRUE){
        laplacian <- laplacian_matrix(object)
    }

    #compute weights if weighted graph given
    if (is.weighted(object)) {
        #assign weights to edges (without dependancy on igraph)
        weights <- r_to_py(edge.attributes(object)$weight)
        snn_graph$es$set_attribute_values('weight', weights)
    }

    if(!is.null(get.vertex.attribute(object, "type")) || is.bipartite(object)){
        type <- as.integer(V(object)$type)
        snn_graph$vs$set_attribute_values('type', r_to_py(as.integer(type)))
    }

    # from here is the same as method for matrix
    # would be better to refactor to call from matrix methof

    #compute partitions
    partition <- find_partition(snn_graph, partition_type = partition_type,
                                initial_membership = initial_membership ,
                                weights = weights,
                                node_sizes = node_sizes,
                                resolution_parameter = resolution_parameter,
                                seed = seed,
                                n_iterations = n_iterations,
                                degree_as_node_size = degree_as_node_size
    )
    partition
}


# global reference to python modules (will be initialized in .onLoad)
leidenalg <- NULL
ig <- NULL

.onLoad = function(libname, pkgname) {
    if(!reticulate::py_available()){
        tryCatch({
            if(!("r-reticulate" %in% reticulate::conda_list()$name)){
                reticulate::conda_create(envname = "r-reticulate")
                reticulate::conda_install(envname = "r-reticulate", packages = "conda")
            }
            reticulate::use_python(reticulate::conda_python())
            reticulate::use_condaenv("r-reticulate")
        }, error = function(e){
            print("Unable to set up conda environment r-reticulate")
            print("run in terminal:")
            print("conda init")
            print("conda create -n r-reticulate")
        },
        finally = print("conda environment r-reticulate installed"))
    }
    tryCatch({
        if(reticulate::py_available() || sum("r-reticulate" == reticulate::conda_list()$name) == 1){
            install_python_modules <- function(method = "auto", conda = "auto") {
                if(!is.null(reticulate::conda_binary())){
                    reticulate::use_python(reticulate::conda_python())
                    if(!("r-reticulate" %in% reticulate::conda_list()$name)){
                        reticulate::conda_create(envname = "r-reticulate", )
                        reticulate::conda_install(envname = "r-reticulate", packages = "conda")
                    }
                    reticulate::use_condaenv("r-reticulate")
                    if(.Platform$OS.type == "windows"){
                        install.packages("devtools",  quiet = TRUE)
                        devtools::install_github("rstudio/reticulate", ref = "86ebb56",  quiet = TRUE)
                        reticulate::conda_install(envname = "r-reticulate", packages = "python-igraph")
                        reticulate::conda_install(envname = "r-reticulate", packages = "mkl", channel = "intel")
                        reticulate::conda_install(envname = "r-reticulate", packages = "leidenalg", channel = "conda-forge")
                        install.packages("reticulate",  quiet = TRUE)
                        reticulate::conda_install(envname = "r-reticulate", packages = "leidenalg") #, channel = "conda-forge")
                        utils::install.packages("reticulate",  quiet = TRUE)
                    } else {
                        reticulate::conda_install("r-reticulate", "python-igraph")
                        reticulate::conda_install("r-reticulate", "leidenalg", forge = TRUE)
                    }
                } else {
                    # shell <- strsplit(Sys.getenv("SHELL"), "/")[[1]]
                    # shell <- shell[length(shell)]
                    # eval(parse(text = paste0(c('system("conda init ', shell, '")'), collapse = "")))
                    # eval(parse(text = paste0(c('system("source ~/.', shell, 'rc")'), collapse = "")))
                    # shell <- as.list(system("echo $0"))
                    # if(shell == sh) shell <- "bash"
                    # system("conda init")
                    # eval(parse(text = paste0(c('system("source ~/.', shell, '_profile")'), collapse = "")))
                    # system("conda init")
                    # system("conda activate r-reticulate")
                    reticulate::py_install("python-igraph", method = method, conda = conda)
                    reticulate::py_install("leidenalg", method = method, conda = conda, forge = TRUE)
                }
            }
        }
    }, error = function(e){
        print("Unable to install python modules igraph and leidenalg")
        print("run in terminal:")
        print("conda install -c conda-forge vtraag python-igraph")
    },
    finally = print("python modules igraph and leidenalg installed"))
    if (suppressWarnings(suppressMessages(requireNamespace("reticulate")))) {
        modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
        if (modules) {
            ## assignment in parent environment!
            leidenalg <<- reticulate::import("leidenalg", delay_load = TRUE)
            ig <<- reticulate::import("igraph", delay_load = TRUE)
        }
    }
}
