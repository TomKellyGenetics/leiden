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
##' @param legacy (defaults to FALSE). Force calling python implementation via reticulate. Default behaviour is calling cluster_leiden in igraph with Modularity (for undirected graphs) and CPM cost functions.
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
                   laplacian = FALSE,
                   legacy = FALSE) {
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
                          laplacian = FALSE,
                          legacy = FALSE
) {
    # disable printing numerals in scientific notation
    oo <- options(scipen = 100000000000)
    # restore options when function terminates
    on.exit(options(oo))

    if(length(partition_type) > 1) partition_type <- partition_type[[1]][1]
    partition_type <- match.arg(partition_type)

    #import python modules with reticulate
    numpy <- import("numpy", delay_load = TRUE)
    leidenalg <- import("leidenalg", delay_load = TRUE)
    ig <- import("igraph", delay_load = TRUE)
    pd <- import("pandas", delay_load = TRUE)

    #convert matrix input (corrects for sparse matrix input)
    if(is.matrix(object) || is(object, "dgCMatrix")){
        object <- object
    } else{
        object <- as.matrix(object)
    }

    #compute weights if non-binary adjacency matrix given
    is_pure_adj <- all(as.logical(unlist(object)) == object)
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
                          laplacian = FALSE,
                          legacy = FALSE
) {
    if(length(partition_type) > 1) partition_type <- partition_type[[1]][1]
    partition_type <- match.arg(partition_type)

    # disable printing numerals in scientific notation
    oo <- options(scipen = 100000000000)
    # restore options when function terminates
    on.exit(options(oo))

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
                  laplacian = laplacian,
                  legacy = legacy
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
                        laplacian = FALSE,
                        legacy = FALSE
) {
    # disable printing numerals in scientific notation
    oo <- options(scipen = 100000000000)
    # restore options when function terminates
    on.exit(options(oo))

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
                                   laplacian = laplacian,
                                   legacy = legacy
        )
    } else{

        #import python modules with reticulate
        numpy <- reticulate::import("numpy", delay_load = TRUE)
        leidenalg <- import("leidenalg", delay_load = TRUE)
        ig <- import("igraph", delay_load = TRUE)

        object <- lapply(object, function(graph){
            if(is.matrix(graph) || is(graph, "dgCMatrix")){
                graph_from_adjacency_matrix(graph)
            } else {
                graph
            }
        })

        names(object) <- c()

        py_list <- r_to_py(lapply(object, function(r_graph){
            make_py_graph(r_graph, weights = weights)
        }))


        if(partition_type == 'ModularityVertexPartition.Bipartite') partition_type <- "ModularityVertexPartition"
        if(partition_type == 'CPMVertexPartition.Bipartite') partition_type <- "CPMVertexPartition"
        partition_type <- gsub(".Bipartite", "", partition_type)
        partition_type <- gsub(".Multiplex", "", partition_type)

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

##' @importFrom igraph V as_edgelist is_weighted is_named edge_attr as_adjacency_matrix laplacian_matrix vertex_attr is_bipartite bipartite_mapping set_vertex_attr simplify as.undirected is_directed communities membership cluster_leiden which_mutual which_loop is_named is_weighted
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
                          laplacian = FALSE,
                          legacy = FALSE
) {
    #default partition
    if(length(partition_type) > 1) partition_type <- partition_type[[1]][1]
    partition_type <- match.arg(partition_type)

    # disable printing numerals in scientific notation
    oo <- options(scipen = 100000000000)
    # restore options when function terminates
    on.exit(options(oo))

    #import python modules with reticulate
    numpy <- reticulate::import("numpy", delay_load = TRUE)
    leidenalg <- import("leidenalg", delay_load = TRUE)
    ig <- import("igraph", delay_load = TRUE)

    #pass weights to igraph if not found
    if(!is_weighted(object) && !is.null(weights)){
        #assign weights to edges (without dependancy on igraph)
        if(length(E(object)) == length(weights)){
            set_edge_attr(object, "weight", value = weights)
        } else {
            warning(paste("weights but be same length as number of edges:", length(E(object))))
            weights <- NULL
        }
    }

    #derive Laplacian
    if(laplacian == TRUE){
        object <- simplify(object, remove.multiple = TRUE, remove.loops = TRUE)
        laplacian <- laplacian_matrix(object)
        if(!is_weighted(object)){
            object <- set_edge_attr(object, "weight", value = -as.matrix(laplacian)[as.matrix(laplacian) < 0])
        }
    }

    #check whether compatible with igraph implementations in R
    if(is_directed(object) && !is_bipartite(object)){
        #coerce to undirected graph object if possible
        if(all(which_mutual(object) | which_loop(object)) || partition_type == "CPMVertexPartition"){
            object <- as.undirected(object, mode = "each")
        }
    }
    call_igraph <- !is_directed(object) && !is_bipartite(object) && legacy == FALSE && (partition_type == "CPMVertexPartition" || partition_type == "ModularityVertexPartition")
    #print(call_igraph)

    if(call_igraph == TRUE){
        #call igraph implementation
        if(partition_type == "CPMVertexPartition"){
            objective_function <- "cpm"
        }
        if(partition_type == "ModularityVertexPartition"){
            objective_function <- "modularity"
        }

        #compute partitions with igraph in C
        if(!is.null(seed)) set.seed(seed)
        partition <- membership(cluster_leiden(graph = object,
                                               objective_function = objective_function,
                                               weights = weights,
                                               resolution_parameter = resolution_parameter,
                                               initial_membership = initial_membership,
                                               n_iterations = n_iterations,
                                               vertex_weights = NULL
        ))
        partition <- as.numeric(partition)
    } else {
        #call python reticulate implementation
        #import python modules with reticulate
        numpy <- import("numpy", delay_load = TRUE)
        leidenalg <- import("leidenalg", delay_load = TRUE)
        ig <- import("igraph", delay_load = TRUE)

        py_graph <- make_py_object(object, weights = weights)

        if(is_bipartite(object) && partition_type == "ModularityVertexPartition"){
            partition_type <- "ModularityVertexPartition.Bipartite"
        }
        if(partition_type == "ModularityVertexPartition.Bipartite"){
            if(is.null(vertex_attr(object, "type"))){
                if(bipartite_mapping(object)$res){
                    packageStartupMessage("computing bipartite partitions")
                    object <- set_vertex_attr(object, "type", value = bipartite_mapping(object)$type)
                    partition_type <- "ModularityVertexPartition.Bipartite"
                } else {
                    packageStartupMessage("cannot compute bipartite types, defaulting to partition type ModularityVertexPartition")
                    partition_type <- "ModularityVertexPartition"
                }
            }
        }
        if(is_bipartite(object) && partition_type == "CPMVertexPartition"){
            partition_type <- "CPMVertexPartition.Bipartite"
        }
        if(partition_type == "CPMVertexPartition.Bipartite"){
            if(is.null(vertex_attr(object, "type"))){
                if(bipartite_mapping(object)$res){
                    packageStartupMessage("computing bipartite partitions")
                    object <- set_vertex_attr(object, "type", value = bipartite_mapping(object)$type)
                    partition_type <- "CPMVertexPartition.Bipartite"
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

        #compute partitions with reticulate
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
    }
    partition
}


# global reference to python modules (will be initialized in .onLoad)
leidenalg <- NULL
ig <- NULL
numpy <- NULL
pd <- NULL

#' @importFrom utils install.packages capture.output

.onAttach <- function(libname, pkgname) {
    if(!reticulate::py_available()){
        tryCatch({
            is.reticulate.env <- any(grepl("r-reticulate/bin", reticulate::conda_list()$python))
            # create conda env if no base image found
            if(!(is.reticulate.env)){
                if(interactive()){
                    install.deps <- readline("create conda environment (yes/no)?")
                    packageStartupMessage(install.deps)
                } else {
                    packageStartupMessage("create conda environment (yes/no)?")
                    install.deps <- "no (use interactive mode)"
                    packageStartupMessage("no (use interactive mode)")
                }
                if(install.deps == "yes" || install.deps == "y"){
                    reticulate::miniconda_update()
                    reticulate::conda_create(envname = "r-reticulate")
                    reticulate::conda_install(envname = "r-reticulate", packages = "conda")
                }
            }
            # use "r-reticulate" or "base" image (which ever is used by reticulate if installed already)
            reticulate.env <- reticulate::conda_list()$name[grep("r-reticulate/bin/python", reticulate::conda_list()$python)][1]
            packageStartupMessage(paste(c("using environment:",  reticulate.env), collapse = " "))
            suppressWarnings(suppressMessages(reticulate::use_python(reticulate::conda_python())))
            suppressWarnings(suppressMessages(reticulate::use_condaenv(reticulate.env)))
        }, error = function(e){
            packageStartupMessage("Unable to set up conda environment r-reticulate")
            packageStartupMessage("run in terminal:")
            packageStartupMessage("conda init")
            packageStartupMessage("conda create -n r-reticulate")
        },
        finally = packageStartupMessage("conda environment r-reticulate installed"))
    }
    tryCatch({
        is.reticulate.env <- any(grepl("r-reticulate/bin", reticulate::conda_list()$python))
        if(reticulate::py_available() || is.reticulate.env ){
            if(interactive()){
                install.deps <- readline("install dependencies (yes/no)?")
                packageStartupMessage(install.deps)
            } else {
                packageStartupMessage("install dependencies (yes/no)?")
                install.deps <- "no (use interactive mode)"
                packageStartupMessage("no (use interactive mode)")
            }
            if(install.deps == "yes" || install.deps == "y"){
                reticulate.env <- reticulate::conda_list()$name[grep("r-reticulate/bin/python", reticulate::conda_list()$python)][1]
                packageStartupMessage(paste(c("using environment:",  reticulate.env), collapse = " "))
                install_python_modules <- function(method = "auto", conda = "auto") {
                    if(!is.null(reticulate::conda_binary())){
                        reticulate::use_python(reticulate::conda_python())
                        if(!(is.reticulate.env)){
                            reticulate::conda_create(envname = reticulate.env)
                            if(!reticulate::py_module_available("conda")) reticulate::conda_install(envname = reticulate.env, packages = "conda")
                        }
                        suppressWarnings(suppressMessages(reticulate::use_condaenv(reticulate.env)))
                        if(.Platform$OS.type == "windows"){
                            utils::install.packages("devtools",  quiet = TRUE)
                            devtools::install_github("rstudio/reticulate", ref = "86ebb56",  quiet = TRUE)
                            if(!reticulate::py_module_available("numpy")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = reticulate.env, packages = "numpy")))
                            if(!reticulate::py_module_available("pandas")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = reticulate.env, packages = "pandas")))
                            if(!reticulate::py_module_available("igraph")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = reticulate.env, packages = "python-igraph")))
                            if(!reticulate::py_module_available("mkl")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = reticulate.env, packages = "mkl", channel = "intel")))
                            if(!reticulate::py_module_available("umap")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = reticulate.env, packages = "umap-learn", channel = "conda-forge")))
                            if(!reticulate::py_module_available("leidenalg")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = reticulate.env, packages = "leidenalg", channel = "conda-forge")))
                            install.packages("reticulate",  quiet = TRUE)
                            if(!reticulate::py_module_available("leidenalg")) suppressWarnings(suppressMessages(reticulate::conda_install(envname = reticulate.env, packages = "leidenalg"))) #, channel = "conda-forge")
                            utils::install.packages("reticulate",  quiet = TRUE)
                        } else {
                            if(!reticulate::py_module_available("numpy")) suppressWarnings(suppressMessages(reticulate::conda_install(reticulate.env, "numpy")))
                            if(!reticulate::py_module_available("pandas")) suppressWarnings(suppressMessages(reticulate::conda_install(reticulate.env, "pandas")))
                            if(!reticulate::py_module_available("igraph")) suppressWarnings(suppressMessages(reticulate::conda_install(reticulate.env, "python-igraph")))
                            if(!reticulate::py_module_available("umap")) suppressWarnings(suppressMessages(reticulate::conda_install(reticulate.env, "umap-learn", forge = TRUE)))
                            if(!reticulate::py_module_available("leidenalg")) suppressWarnings(suppressMessages(reticulate::conda_install(reticulate.env, "leidenalg", forge = TRUE)))
                            #Sys.setenv(PATH = paste0(strsplit(reticulate::py_config()$pythonhome, ":")[[1]][1], "/bin:$PATH"))
                            Sys.setenv(RETICULATE_PYTHON = reticulate::conda_python())
                        }
                    } else {
                        if(!reticulate::py_module_available("numpy")) suppressWarnings(suppressMessages(reticulate::py_install("numpy")))
                        if(!reticulate::py_module_available("pandas")) suppressWarnings(suppressMessages(reticulate::py_install("pandas")))
                        if(!reticulate::py_module_available("igraph")) suppressWarnings(suppressMessages(reticulate::py_install("python-igraph", method = method, conda = conda)))
                        if(!reticulate::py_module_available("umap")) suppressWarnings(suppressMessages(reticulate::py_install("umap-learn")))
                        if(!reticulate::py_module_available("leidenalg")) suppressWarnings(suppressMessages(reticulate::py_install("leidenalg", method = method, conda = conda, forge = TRUE)))
                        Sys.setenv(RETICULATE_PYTHON = reticulate::py_config()$python)
                    }
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
        packageStartupMessage("conda install -n r-reticulate -c conda-forge leidenalg python-igraph pandas umap-learn")
    },
    finally = packageStartupMessage("python modules igraph and leidenalg installed"))
    if (suppressWarnings(suppressMessages(requireNamespace("reticulate")))) {
        modules <- reticulate::py_module_available("pandas")
        if (modules) {
            ## assignment in parent environment!
            pd <- reticulate::import("pandas", delay_load = TRUE)
        }
        modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")
        if (modules) {
            ## assignment in parent environment!
            numpy <- reticulate::import("numpy", delay_load = TRUE)
            leidenalg <- reticulate::import("leidenalg", delay_load = TRUE)
            ig <- reticulate::import("igraph", delay_load = TRUE)
        }
    }
}

