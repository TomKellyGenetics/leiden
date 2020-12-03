#' @include find_partition.R
#'
NULL

##' Run Leiden clustering algorithm
##'
##' @description Implements the Leiden clustering algorithm in R using reticulate to run the Python version. Requires the python "leidenalg" and "igraph" modules to be installed. Returns a vector of partition indices.
##' Windows users can still this with devtools::install_github("rstudio/reticulate", ref = "86ebb56"); reticulate::use_condaenv("r-reticulate"); reticulate::conda_install("r-reticulate", "leidenalg", channel = "vtraag")
##' @param object An adjacency matrix compatible with \code{\link[igraph]{igraph}} object or an input graph as an \code{\link[igraph]{igraph}} object (e.g., shared nearest neighbours). A list of multiple graph objects can be passed for multiplex community detection.
##' @param partition_type Type of partition to use. Defaults to RBConfigurationVertexPartition. Options include: ModularityVertexPartition, RBERVertexPartition, CPMVertexPartition, MutableVertexPartition, SignificanceVertexPartition, SurpriseVertexPartition, ModularityVertexPartition.Bipartite, CPMVertexPartition.Bipartite (see the Leiden python module documentation for more details)
##' @param initial_membership,weights,node_sizes Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None). Weights are derived from weighted igraph objects and non-zero integer values of adjacency matrices.
##' @param resolution_parameter A parameter controlling the coarseness of the clusters
##' @param seed Seed for the random number generator. By default uses a random seed if nothing is specified.
##' @param n_iterations Number of iterations to run the Leiden algorithm. By default, 2 iterations are run. If the number of iterations is negative, the Leiden algorithm is run until an iteration in which there was no improvement.
##' @param max_comm_size (non-negative int) – Maximal total size of nodes in a community. If zero (the default), then communities can be of any size.
##' @param degree_as_node_size (defaults to FALSE). If True use degree as node size instead of 1, to mimic modularity for Bipartite graphs.
##' @param laplacian (defaults to FALSE). Derive edge weights from the Laplacian matrix.
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
##' my_graph <- graph_from_adjacency_matrix(adjacency_matrix)
##' partition <- leiden(my_graph)
##' table(partition)
##'
##' # generate (weighted) igraph object in R
##' library("igraph")
##' adjacency_matrix[adjacency_matrix >= 1] <- weights
##' my_graph <- graph_from_adjacency_matrix(adjacency_matrix, weighted = TRUE)
##' partition <- leiden(my_graph)
##' table(partition)
##'
##' # pass weights to python leidenalg
##' adjacency_matrix[adjacency_matrix >= 1 ] <- 1
##' my_graph <- graph_from_adjacency_matrix(adjacency_matrix, weighted = NULL)
##' weights <- sample(1:10, sum(adjacency_matrix!=0), replace=TRUE)
##' partition <- leiden(my_graph, weights = weights)
##' table(partition)
##'
##' # run only if python is available (for testing)
##' }
##'
##' @keywords graph network igraph mvtnorm simulation
##' @importFrom reticulate import py_to_r r_to_py py_has_attr py_get_attr py_set_attr
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
                   max_comm_size = 0L,
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
                          max_comm_size = 0L,
                          degree_as_node_size = FALSE,
                          laplacian = FALSE
) {
    if(length(partition_type) > 1) partition_type <- partition_type[[1]][1]
    partition_type <- match.arg(partition_type)

    #import python modules with reticulate
    numpy <- import("numpy", delay_load = TRUE)
    leidenalg <- import("leidenalg", delay_load = TRUE)
    ig <- import("igraph", delay_load = TRUE)

    #convert matrix input (corrects for sparse matrix input)
    if(is.matrix(object) || is(object, "dgCMatrix")){
        object <- object
    } else{
        object <- as.matrix(object)
    }

    #compute weights if non-binary adjacency matrix given
    is_pure_adj <- all(as.logical(object) == object)
    if (is.null(weights) && !is_pure_adj) {
        if(!is.matrix(object)) object <- as.matrix(object)
        #assign weights to edges (without dependancy on igraph)
        t_mat <- t(object)
        weights <- t_mat[t_mat!=0]
        #remove zeroes from rows of matrix and return vector of length edges
    }

    py_graph <- make_py_graph(object, weights = weights)

    #compute partitions
    partition <- find_partition(py_graph, partition_type = partition_type,
                                initial_membership = initial_membership,
                                weights = weights,
                                node_sizes = node_sizes,
                                resolution_parameter = resolution_parameter,
                                seed = seed,
                                n_iterations = n_iterations,
                                max_comm_size = max_comm_size,
                                degree_as_node_size = degree_as_node_size
    )
    partition
}

##' @export
leiden.data.frame <- leiden.matrix

##' @importFrom igraph graph_from_adjacency_matrix edge_attr set_edge_attr E
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
                          max_comm_size = 0L,
                          degree_as_node_size = FALSE,
                          laplacian = FALSE
) {
    #cast to sparse matrix
    object <- as(object, "dgCMatrix")
    #run as igraph object (passes to reticulate)
    if(is.null(weights)){
        object <- graph_from_adjacency_matrix(adjmatrix = object, weighted = TRUE)
        weights <- edge_attr(object)$weight
    } else {
        object <- graph_from_adjacency_matrix(adjmatrix = object, weighted = TRUE)
        object <- set_edge_attr(object, "weight", index=E(object), weights)
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

##' @importFrom igraph graph_from_adjacency_matrix edge_attr set_edge_attr E is.igraph
##' @importFrom methods as
##' @importClassesFrom Matrix dgCMatrix dgeMatrix
##' @export
leiden.list <- function(object,
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
                          max_comm_size = 0L,
                          degree_as_node_size = FALSE,
                          laplacian = FALSE
) {
    if(length(partition_type) > 1) partition_type <- partition_type[[1]][1]
    partition_type <- match.arg(partition_type)

    if(length(object) == 1 || is.igraph(object)){
        if(!is.igraph(object) && is.list(object)) object <- object[[1]]
        partition <- leiden.igraph(object,
                            partition_type = partition_type,
                            weights = weights,
                            node_sizes = node_sizes,
                            resolution_parameter = resolution_parameter,
                            seed = seed,
                            n_iterations = n_iterations,
                            max_comm_size = max_comm_size,
                            degree_as_node_size = degree_as_node_size,
                            laplacian = laplacian
        )
    } else{

        #import python modules with reticulate
        numpy <- reticulate::import("numpy", delay_load = TRUE)
        leidenalg <- import("leidenalg", delay_load = TRUE)
        ig <- import("igraph", delay_load = TRUE)

        py_list <- r_to_py(lapply(object, function(r_graph){
            make_py_graph(r_graph, weights = weights)
        }))


        if(partition_type == 'ModularityVertexPartition.Bipartite') partition_type <- "ModularityVertexPartition"
        if(partition_type == 'CPMVertexPartition.Bipartite') partition_type <- "CPMVertexPartition"

        #compute partitions with reticulate
        partition <- find_partition_multiplex(py_list, partition_type = partition_type,
                                    initial_membership = initial_membership,
                                    weights = weights,
                                    node_sizes = node_sizes,
                                    resolution_parameter = resolution_parameter,
                                    seed = seed,
                                    n_iterations = n_iterations,
                                    max_comm_size = max_comm_size,
                                    degree_as_node_size = degree_as_node_size
        )
    }
    partition
}

##' @export
leiden.default <- leiden.matrix

##' @importFrom igraph V as_edgelist is.weighted is.named edge_attr as_adjacency_matrix laplacian_matrix vertex_attr is_bipartite bipartite_mapping set_vertex_attr simplify
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
                          max_comm_size = 0L,
                          degree_as_node_size = FALSE,
                          laplacian = FALSE
) {
    if(length(partition_type) > 1) partition_type <- partition_type[[1]][1]
    partition_type <- match.arg(partition_type)

    #import python modules with reticulate
    numpy <- reticulate::import("numpy", delay_load = TRUE)
    leidenalg <- import("leidenalg", delay_load = TRUE)
    ig <- import("igraph", delay_load = TRUE)

    #default partition
    if(length(partition_type) > 1) partition_type <- partition_type[[1]][1]
    partition_type <- match.arg(partition_type)

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

    #derive Laplacian
    if(laplacian == TRUE){
        object <- simplify(object, remove.multiple = TRUE, remove.loops = TRUE)
        laplacian <- laplacian_matrix(object)
        if(!is.weighted(object)){
            edge_attr(object)$weight
            object <- set_edge_attr(object, "weight", value = -as.matrix(laplacian)[as.matrix(laplacian) < 0])
        }
    }

    py_graph <- make_py_graph(object, weights = weights)

    if(length(partition_type) > 1) partition_type <- partition_type[1]
    if(partition_type == "ModularityVertexPartition.Bipartite"){
        if(is.null(vertex_attr(object, "type"))){
            if(bipartite_mapping(object)$res){
                packageStartupMessage("computing bipartite partitions")
                object <- set_vertex_attr(object, "type", value = bipartite_mapping(object)$type)
            } else {
                packageStartupMessage("cannot compute bipartite types, defaulting to partition type ModularityVertexPartition")
                partition_type <- "ModularityVertexPartition"
            }
        }
    }
    if(partition_type == "CPMVertexPartition.Bipartite"){
        if(is.null(vertex_attr(object, "type"))){
            if(bipartite_mapping(object)$res){
                packageStartupMessage("computing bipartite partitions")
                object <- set_vertex_attr(object, "type", value = bipartite_mapping(object)$type)
            } else {
                packageStartupMessage("cannot compute bipartite types, defaulting to partition type CPMVertexPartition")
                partition_type <- "CPMVertexPartition"
            }
        }
    }

    if(!is.null(vertex_attr(object, "type")) || is_bipartite(object)){
        type <- as.integer(unlist(V(object)$type))
        py_graph$vs$set_attribute_values('type', r_to_py(as.integer(type)))
    }

    #compute partitions
    partition <- find_partition(py_graph, partition_type = partition_type,
                                initial_membership = initial_membership ,
                                weights = weights,
                                node_sizes = node_sizes,
                                resolution_parameter = resolution_parameter,
                                seed = seed,
                                n_iterations = n_iterations,
                                max_comm_size = max_comm_size,
                                degree_as_node_size = degree_as_node_size
    )
    partition
}


# global reference to python modules (will be initialized in .onLoad)
leidenalg <- NULL
ig <- NULL
numpy <- NULL

#' @importFrom utils install.packages capture.output

.onAttach <- function(libname, pkgname) {
    if(!reticulate::py_available()){
        tryCatch({
            if(!("r-reticulate" %in% reticulate::conda_list()$name)){
                reticulate::conda_create(envname = "r-reticulate")
                reticulate::conda_install(envname = "r-reticulate", packages = "conda")
            }
            suppressWarnings(suppressMessages(reticulate::use_python(reticulate::conda_python())))
            suppressWarnings(suppressMessages(reticulate::use_condaenv("r-reticulate")))
        }, error = function(e){
            packageStartupMessage("Unable to set up conda environment r-reticulate")
            packageStartupMessage("run in terminal:")
            packageStartupMessage("conda init")
            packageStartupMessage("conda create -n r-reticulate")
        },
        finally = packageStartupMessage("conda environment r-reticulate installed"))
    }
    tryCatch({
        if(reticulate::py_available() || sum("r-reticulate" == reticulate::conda_list()$name) >= 1){
            install_python_modules <- function(method = "auto", conda = "auto") {
                if(!is.null(reticulate::conda_binary())){
                    reticulate::use_python(reticulate::conda_python())
                    if(!("r-reticulate" %in% reticulate::conda_list()$name)){
                        reticulate::conda_create(envname = "r-reticulate", )
                        if(!reticulate::py_module_available("conda")) reticulate::conda_install(envname = "r-reticulate", packages = "conda")
                    }
                    suppressWarnings(suppressMessages(reticulate::use_condaenv("r-reticulate")))
                    if(.Platform$OS.type == "windows"){
                        install.packages("devtools",  quiet = TRUE)
                        devtools::install_github("rstudio/reticulate", ref = "86ebb56",  quiet = TRUE)
                        if(!reticulate::py_module_available("numpy")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = "r-reticulate", packages = "numpy")))
                        if(!reticulate::py_module_available("pandas")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = "r-reticulate", packages = "pandas")))
                        if(!reticulate::py_module_available("igraph")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = "r-reticulate", packages = "python-igraph")))
                        if(!reticulate::py_module_available("mkl")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = "r-reticulate", packages = "mkl", channel = "intel")))
                        if(!reticulate::py_module_available("umap")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = "r-reticulate", packages = "umap-learn", channel = "conda-forge")))
                        if(!reticulate::py_module_available("leidenalg")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = "r-reticulate", packages = "leidenalg", channel = "conda-forge")))
                        install.packages("reticulate",  quiet = TRUE)
                        if(!reticulate::py_module_available("leidenalg")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = "r-reticulate", packages = "leidenalg"))) #, channel = "conda-forge")
                        utils::install.packages("reticulate",  quiet = TRUE)
                    } else {
                        if(!reticulate::py_module_available("numpy")) suppressWarnings(suppressMessages(reticulate::conda_install("r-reticulate", "numpy")))
                        if(!reticulate::py_module_available("pandas")) suppressWarnings(suppressMessages(reticulate::conda_install("r-reticulate", "pandas")))
                        if(!reticulate::py_module_available("igraph")) suppressWarnings(suppressMessages(reticulate::conda_install("r-reticulate", "python-igraph")))
                        if(!reticulate::py_module_available("umap")) suppressWarnings(suppressMessages(reticulate::conda_install("r-reticulate", "umap-learn", forge = TRUE)))
                        if(!reticulate::py_module_available("leidenalg")) suppressWarnings(suppressMessages(reticulate::conda_install("r-reticulate", "leidenalg", forge = TRUE)))
                        #Sys.setenv(PATH = paste0(strsplit(reticulate::py_config()$pythonhome, ":")[[1]][1], "/bin:$PATH"))
                        Sys.setenv(RETICULATE_PYTHON = reticulate::conda_python())
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
                    if(!reticulate::py_module_available("numpy")) suppressWarnings(suppressMessages(reticulate::py_install("numpy")))
                    if(!reticulate::py_module_available("pandas")) suppressWarnings(suppressMessages(reticulate::py_install("pandas")))
                    if(!reticulate::py_module_available("igraph")) suppressWarnings(suppressMessages(reticulate::py_install("python-igraph", method = method, conda = conda)))
                    if(!reticulate::py_module_available("umap")) suppressWarnings(suppressMessages(reticulate::py_install("umap-learn")))
                    if(!reticulate::py_module_available("leidenalg")) suppressWarnings(suppressMessages(reticulate::py_install("leidenalg", method = method, conda = conda, forge = TRUE)))
                    #Sys.setenv(PATH = paste0(strsplit(reticulate::py_config()$pythonhome, ":")[[1]][1], "/bin:$PATH"))
                    Sys.setenv(RETICULATE_PYTHON = reticulate::py_config()$python)
                }
            }
            quiet <- function(expr, all = TRUE) {
                if (Sys.info()['sysname'] == "Windows") {
                    file <- "NUL"
                } else {
                    file <- "/dev/null"
                }

                if (all) {
                    suppressWarnings(suppressMessages(suppressPackageStartupMessages(
                        capture.output(expr, file = file)
                    )))
                } else {
                    capture.output(expr, file = file)
                }

            }
            quiet(install_python_modules())
        }
    }, error = function(e){
        packageStartupMessage("Unable to install python modules igraph and leidenalg")
        packageStartupMessage("run in terminal:")
        packageStartupMessage("conda install -c conda-forge vtraag python-igraph")
    },
    finally = packageStartupMessage("python modules igraph and leidenalg installed"))
    if (suppressWarnings(suppressMessages(requireNamespace("reticulate")))) {
        modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
        if (modules) {
            ## assignment in parent environment!
            numpy <- reticulate::import("numpy", delay_load = TRUE)
            pd <- reticulate::import("pandas", delay_load = TRUE)
            leidenalg <- reticulate::import("leidenalg", delay_load = TRUE)
            ig <- reticulate::import("igraph", delay_load = TRUE)
        }
    }
}
