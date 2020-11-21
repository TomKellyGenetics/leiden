##' convert to python object
##' @param object an igraph object or matrix
##' @param weights Parameters to pass to the Python leidenalg function (defaults initial_membership=None, weights=None). Weights are derived from weighted igraph objects and non-zero integer values of adjacency matrices.
##' @noRd
##' @description internal function to compute partitions by calling Python with reticulate
##' @keywords internal
make_py_object <- function(object, weights = NULL) {
  UseMethod("make_py_object", object)
}

make_py_object.matrix <- function(object, weights = NULL){
  #import python modules with reticulate
  numpy <- reticulate::import("numpy", delay_load = TRUE)
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
  adj_mat_py <- r_to_py(adj_mat, convert = TRUE)
  if(is(object, "dgCMatrix")){
    adj_mat_py <- adj_mat_py$toarray()
  }
  adj_mat_py <- adj_mat_py$tolist()

  adj_mat_py
}

make_py_object.data.frame <- function(object, weights = NULL){
  pd <- reticulate::import("pandas", delay_load = TRUE)
  adj_df_py <- pd$DataFrame(data = r_to_py(object, convert = TRUE))

  adj_df_py
}

make_py_object.igraph <- function(object, weights = NULL){
  #import python modules with reticulate
  numpy <- import("numpy", delay_load = TRUE)
  leidenalg <- import("leidenalg", delay_load = TRUE)
  ig <- import("igraph", delay_load = TRUE)

  ##convert to python numpy.ndarray, then a list
  if(!is_named(object)){
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

  py_graph <- ig$Graph()
  py_graph$add_vertices(r_to_py(vertices))
  py_graph$add_edges(r_to_py(edgelist))

  #compute weights if weighted graph given
  if(is_weighted(object)){
    #assign weights to edges (without dependancy on igraph)
    weights <- r_to_py(edge_attr(object)$weight)
    py_graph$es$set_attribute_values('weight', weights)
  }
  py_graph
}

make_py_graph <- function(object, weights = NULL) {
  UseMethod("make_py_graph", object)
}

make_py_graph.data.frame <- function(object, weights = NULL){
  py_graph <- make_py_graph(as.matrix(object), weights = weights)
}

make_py_graph.matrix <- function(object, weights = NULL){
  #compute weights if non-binary adjacency matrix given
  is_pure_adj <- all(as.logical(adj_mat) == adj_mat)
  if (is.null(weights) && !is_pure_adj) {
    if(!is.matrix(object)) adj_mat <- as.matrix(adj_mat)
    #assign weights to edges (without dependancy on igraph)
    t_mat <- t(adj_mat)
    weights <- t_mat[t_mat!=0]
    #remove zeroes from rows of matrix and return vector of length edges
  }

  #import python modules with reticulate
  numpy <- reticulate::import("numpy", delay_load = TRUE)
  leidenalg <- import("leidenalg", delay_load = TRUE)
  ig <- import("igraph", delay_load = TRUE)

  adj_mat_py <- make_py_object(object, weights = weights)

  #convert graph structure to a Python compatible object
  GraphClass <- if (!is.null(weights) && !is_pure_adj){
    ig$Graph$Weighted_Adjacency
  } else {
    ig$Graph$Adjacency
  }
  py_graph <- GraphClass(adj_mat_py)
}

make_py_graph.igraph <- make_py_object.igraph

