# module PartitionedGraphs
#
# using SparseArrays
# using LinearAlgebra

##' @importClassesFrom Matrix dgCMatrix

Partition <- function() as.integer(c(c()))
WeightMatrix <- function() as(matrix(), "dgCMatrix")

# setClass("PartitionedGraph",
#          slots = c(
#            edge_weight = "matrix",
#            cardinality = "integer",
#            partition = "integer",
#            size = "integer",
#            membership = "integer"
#          )
# )

##' @import R6

PartitionedGraph <- R6Class("PartitionedGraph",
                            list(
                              edge_weight = c(),
                              cardinality = NULL,
                              partition = NULL,
                              size = NULL,
                              membership = NULL,
                              initialize = function(adjmat,
                                                    cardinality = create_unit_cardinality(nv(adjmat)),
                                                    partition = create_singleton_partition(nv(adjmat))
                                                    ){
                                adjmat <- as(adjmat, "dgCMatrix")
                                stopifnot((is.matrix(adjmat) || is(adjmat, "dgCMatrix")))
                                cardinality <- as.integer(cardinality)
                                stopifnot(is.integer(cardinality), length(cardinality) == nv(adjmat))
                                partition <- lapply(partition, as.integer)
                                #print(partition)
                                stopifnot(is.integer(partition[[1]]), length(unlist(partition)) == nv(adjmat))

                                check_adjacent_matrix(adjmat)
                                adjmat <- as(adjmat, "dgCMatrix")
                                self$edge_weight <- adjmat#@x #sparse (edgelist)

                                n <- nv(adjmat)
                                check_cardinality(n, cardinality)
                                self$cardinality <- cardinality

                                check_partition(n, partition)
                                self$partition <- partition

                                self$size <- rep(0L, length(partition))
                                self$membership <- rep(0L, n)
                                for(i in 1:length(partition)){
                                  community <- partition[[i]]
                                  self$size[i] = weighted_sum(community, cardinality)
                                  self$membership[community] <- i
                                }
                                return(self)
                              }
                            )
)

#graph <- new("PartitionedGraph", edge_weight = matrix(c(1L, 2L, 3L, 4L), 2, 2), cardinality = c(1L, 3L))
# graph <- PartitionedGraph$new(adjmat)

# number of vertices
# A::AbstractMatrix
nv <- function(A) UseMethod("nv", A)
nv.matrix <-  function(A) nrow(A)
nv.double <-  function(A) nrow(A)
nv.numeric <-  function(A) length(A)

nv.Matrix <- function(A){
  A <- as(A, "dgCMatrix")
  nrow(A)
}
#graph::PartitionedGraph
nv.PartitionedGraph <- function(graph) nrow(graph$edge_weight)
nv.igraph <- function(G) length(V(G))

# number of communities
nc <- function(graph) UseMethod("nc", graph)
nc.PartitionedGraph <- function(graph) length(graph$partition)
nc.igraph <- function(graph) length(unique(get.vertex.attribute(ig, "partition")))

# neighbors of node `u`
#graph::PartitionedGraph, u::Int
neighbors <- function(graph, u){
  #print(u)
  u <- as.integer(u)
  A <- graph$edge_weight
  A <- as(A, "dgCMatrix")
  na.omit(as.integer(setdiff(as.character(which(apply(A[,max(1, u-1):min(u+1, ncol(A))], 1, sum) > 0)), as.character(u))))
}

# communities connected with `u`
#graph::PartitionedGraph, u::Int
connected_communities <- function (graph, u){
  u <- as.integer(u)
  A <- graph$total_weight
  A <- as(A, "dgCMatrix")
  as.integer(setdiff(as.character(which(apply(A[,max(1, u-1):min(u+1, ncol(A))], 1, sum) > 0)), as.character(u)))
}

# generators of the default parameter
create_unit_cardinality <- function(n) rep(1, as.integer(n))
create_singleton_partition <- function(n) as.list(as.integer(c(1:n)))

#adjmat::AbstractMatrix
check_adjacent_matrix <- function(adjmat){
  n <- nrow(adjmat)
  m <- ncol(adjmat)
  if(m != n) warning("invalid adjacent matrix: not a square matrix")
  if(n == 0) warning("invalid adjacent matrix: empty matrix")
  if(!isSymmetric(adjmat)) warning("invalid adjacent matrix: not an symmetric matrix")
  if(is(adjmat, "Matrix")){
    if(any(adjmat@x < 0)) warning("invalid adjacent matrix: found negative weight(s)")
  } else {
    if(any(adjmat < 0)) warning("invalid adjacent matrix: found negative weight(s)")
  }
  return(NULL)
}

#n::Integer, cardinality::Vector{Int}
check_cardinality <- function(n, cardinality){
  n <- as.integer(n)
  if(n != length(cardinality)) warning("invalid cardinality: mismatching length")
  if(any(cardinality <= 0)) warning("invalid cardinality: found non-positive value")
  return(NULL)
}

#n::Integer, partition::Partition
check_partition <- function(n, partition){
  n <- as.integer(n)
  found <- c()
  for(community in partition){
    for(u in community){
      if(u %in% found) warning("invalid partition: found duplicated node")
      if(!(1 <= u && u <= n)) warning("invalid partition: found out-of-bounds node")
      found <- c(found, u)
    }
  }
  if(length(found) != n) warning("invalid partition: found missing node")
  return(NULL)
}

# Move node `u` to `dst`.
#graph::PartitionedGraph, (u, dst)::Pair{Int,Int}
move_node <- function(graph, u, dst){
  if(!(1 <= dst && dst <= nc(graph) + 1)) warning("(1 <= dst && dst <= nc(graph) + 1)")
  src <- graph$membership[u]
  if(src == dst){
    # no movement
    return(graph)
  }
    cardinality <- graph$cardinality[u]
    community_src <- graph$partition[[src]]
    pos <- grep(u, community_src)
    if(dst > nc(graph)){
      community_dst <- as.integer(c())
      graph$partition <- list(graph$partition, community_dst)
      graph$size <- c(graph$size, 0)
    } else {
      community_dst <- graph$partition[[dst]]
    }
    community_src <- setdiff(community_src, u) #community_src[1:length(community_src) != pos]
    community_dst <- unlist(c(community_dst, u))
    graph$partition[[src]] <-  community_src
    graph$partition[[dst]] <- community_dst
    graph$size[src] <- graph$size[src] - cardinality
    graph$size[dst] <- graph$size[dst] + cardinality
    graph$membership[u] <- dst
    return(graph)
  }


# Drop empty communities from the graph.
#graph::PartitionedGraph
drop_empty_communities <- function(graph){
  empty = as.integer(c())
  for(i in 1:length(graph$partition)){
    community <- graph$partition[[i]]
    if(!((is.null(community) | length(community) == 0) == (graph$size[i] == 0))) warning("empty community must have graph size of zero")
    if(is.null(community) | length(community) == 0){
      empty <- c(empty, i)
    }
  }
  graph$partition<- graph$partition[!(1:length(graph$partition) %in% empty)]
  graph$size <- graph$size[!(1:length(graph$size) %in% empty)]

  for(i in 1:length(graph$partition)){
    community <- graph$partition[[i]]
    graph$membership[community] <- i
  }
  return(graph)
}

#graph::PartitionedGraph, partition::Partition
reset_partition <- function(graph, partition){
  check_partition(nv(graph), partition)
  graph$partition <- as.list(rep(NA, length(partition)))
  graph$size <- rep(NA, length(partition))
  for(i in 1:length(partition)){
    community <- partition[[i]]
    graph$partition[[i]] <- community
    graph$size[i] <- weighted_sum(community, graph$cardinality)
    graph$membership[community] <- i
  }
  return(graph)
}

#xs::Vector{Int}, weights::Vector{Int}
weighted_sum <- function(xs, weights){
  if(is.null(xs) | length(xs) == 0){
    return(0)
  } else {
    return(sum(sapply(as.list(xs), function(x) weights[x])))
  }
}

