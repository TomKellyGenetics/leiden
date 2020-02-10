#module Leiden

# using Random:
#   shuffle,
# shuffle!

##' @importFrom permute shuffle

#   using SparseArrays:
#   sparse,


##' @importFrom Matrix sparseMatrix
##' @importClassesFrom Matrix sparseMatrix
# spzeros
# sparseMatrix(i, j)

# using StatsFuns:
#   logsumexp
#log(sum(exp(x)))

#source("src/PartitionedGraphs.R")
# using .PartitionedGraphs:
#   PartitionedGraph,
# Partition,
# WeightMatrix,
# nv, nc,
# neighbors,
# move_node!,
# drop_empty_communities!,
# reset_partition!,
# create_singleton_partition


# The Leiden algorithm
# --------------------

##' @param resolution numeric 1.0,
##' @param randomness numeric 0.01,
##' @param partition Partition Object

library("Matrix")
A <- as(matrix(0, 6, 6), "dgCMatrix")
edges <- rbind(
  c(1, 2),
  c(2, 3),
  c(3, 1),
  c(4, 5),
  c(5, 6),
  c(6, 4),
  c(1, 4)
)
apply(edges, 1, function(edge){
  u <- edge[1]
  v <- edge[2]
  A[u, v] <<- A[v, u] <<- 1
})
set.seed(1234)
# leiden_naive(A, resolution = 0.25)
set.seed(1234)
# louvain_naive(A, resolution = 0.25)

# ig <- graph_from_adjacency_matrix(A)
# ig <- as.undirected(ig)
# cluster_louvain(ig)

leiden_naive <- function(adjmat,
                         resolution = 1.0,
                         randomness = 0.01,
                         partition = create_singleton_partition(nrow(adjmat))){
  UseMethod("leiden_naive", adjmat)
}

leiden_naive.matrix <- function(adjmat,
                                resolution = 1.0,
                                randomness = 0.01,
                                partition = create_singleton_partition(nrow(adjmat))){
  adjmat <- as(adjmat, "dgCMatrix")
  graph <- PartitionedGraph$new(adjmat, partition = partition)
  return(leiden_naive.default(graph, as.numeric(resolution), as.numeric(randomness)))
}

leiden_naive.Matrix <- function(adjmat,
                                resolution = 1.0,
                                randomness = 0.01,
                                partition = create_singleton_partition(nrow(adjmat))){
  graph <- PartitionedGraph$new(adjmat, partition = partition)
  return(leiden_naive.default(graph, as.numeric(resolution), as.numeric(randomness)))
}


##' @param graph PartitionedGraph Object
##' @param  γ numeric
##' @param θ numeric
leiden_naive.default <- function(graph, γ, θ){
  stack <- Partition()
  loop <- 0
  while(loop < 1){
    # loop <- FALSE
    graph <- move_nodes_fast(graph, γ)
    print(paste0("nv = ", nv(graph),"; nc = ", nc(graph), "; H = ", H(graph, γ)))
    if(nc(graph) != nv(graph)){
      refined <- refine_partition(graph, γ, θ)
      if(nc(refined) == nv(refined)){
        warning(paste("not refined",  nc(refined)))
      } else {
      stack <- list(stack, refined$partition)
      graph_prime <- aggregate_graph(refined)
      partition <- as.list(as.integer(rep(NA, nc(graph))))
      for(i in 1:length(refined$partition)){
        community <- refined$partition[[i]]
        u <- community[1]
        j <- graph$membership[u]
        partition[[j]]<-i
      }
      graph_prime <- reset_partition(graph_prime, partition)
      graph <- graph_prime
      }
      loop <- loop + 1
    }
  }
  stack <- list(stack, graph$partition)
  return(list(quality = H(graph, γ), partition = flatten(stack), community = graph$membership))
}


# The Louvain algorithm
# ---------------------

##' @param resolution numeric 1.0,
##' @param randomness numeric 0.01,
##' @param partition Partition Object
##'


louvain_naive <- function(adjmat,
                          resolution = 1.0,
                          randomness = 0.01,
                          partition = create_singleton_partition(nrow(adjmat))){
  UseMethod("louvain_naive", adjmat)
}

louvain_naive.matrix <- function(adjmat,
                                 resolution = 1.0,
                                 randomness = 0.01,
                                 partition = create_singleton_partition(nrow(adjmat))){
  adjmat <- as(adjmat, "dgCMatrix")
  graph <- PartitionedGraph$new(adjmat, partition = partition)
  return(louvain_naive.default(graph, as.numeric(resolution)))
}

louvain_naive.Matrix <- function(adjmat,
                                 resolution = 1.0,
                                 randomness = 0.01,
                                 partition = create_singleton_partition(nrow(adjmat))){
  graph <- PartitionedGraph$new(adjmat, partition = partition)
  return(louvain_naive.default(graph, as.numeric(resolution)))
}

##' @param graph PartitionedGraph Object
##' @param  γ numeric
louvain_naive.default <- function(graph, γ){
  stack <- Partition()
  loop <- 0
  while(loop < 5){
    graph <- move_nodes(graph, γ)
    print(paste0("nv = ", nv(graph),"; nc = ", nc(graph), "; H = ", H(graph, γ)))
    stack <- list(stack, graph$partition)
    if(nc(graph) != nv(graph)){
      graph <- aggregate_graph(graph)
      loop <- loop + 1
    }
  }
  return(list(quality = H(graph, γ), partition = flatten(stack), community = graph$membership))
}

##' @param graph PartitionedGraph Object
##' @param  γ numeric

move_nodes <- function(graph, γ){
  H_old <- H(graph, γ)
  nodes <- 1:nv(graph)
  connected <- as.integer(NA)
  total_weights <- rep(0.0, nc(graph))
  loop <- TRUE
  while(loop){
    loop <- FALSE
    nodes <- shuffle(as.integer(length(nodes)))
    for(u in nodes){
      # compute total edge weights for each community
      connected <- c()
      for(v in neighbors(graph, u)){
        i <- graph$membership[v]
        if(total_weights[i] == 0){
          connected <- c(connected, i)
        }
        total_weights[i] <- total_weights[i] + graph$edge_weight[u,v]


        # find the best community to which `u` belongs
        c_u <- graph$cardinality[u]
        weight_u <- graph$edge_weight[u,u]
        src <- dst <- graph$membership[[u]]
        weight_src <- total_weights[src]
        size_src <- graph$size[src]
        maxgain <- 0.0
        for(i in connected){
          if(i != src){
            gain <- total_weights[i] + weight_u - weight_src - γ * (graph$size[i] - size_src + c_u) * c_u
            if(gain > maxgain){
              dst <- i
              maxgain = gain
            }
            total_weights[i] = 0
          }
        }
        total_weights[src] = 0

        #@assert all(x == 0 for x in total_weights)
        # if((all(total_weights == 0))){
        #   warning("all total_weights are 0")
        # } else {
        #   print("total weights:\n", total_weights)
        # }

        if(src != dst){
          # move `u` to the best community and add its neighbors to the queue if needed
          graph <- move_node(graph, u, dst)
        }
      }

      H_new <- H(graph, γ)
      if(H_new > H_old){
        H_old <- H_new
        loop <- TRUE
      }
    }
  }
  return(drop_empty_communities(graph))
}

##' @param graph PartitionedGraph Object
##' @param  γ numeric
##' @importFrom permute shuffle

move_nodes_fast <- function(graph, γ){
  n <- nv(graph)
  queue <- shuffle(n)
  queued <- 1:n #BitSet(1:n)
  connected <- as.integer(c())
  total_weights <- rep(0.0, nc(graph))
  while(length(queue) > 0){
    #print(queue)
    u <- queue[1]
    queue <- setdiff(queue, u)
    queued <- setdiff(queued, u)

    # compute total edge weights for each community
    connected <- as.integer(c())
    for(v in neighbors(graph, u)){
      i <- graph$membership[v]
      if(total_weights[i] == 0){
        connected <- c(connected, i)
      }
      total_weights[i] <- total_weights[i] + graph$edge_weight[u,v]
    }

    # find the best community to which `u` belongs
    c_u <- graph$cardinality[u]
    weight_u <- graph$edge_weight[u,u]
    src <- dst <- graph$membership[u]
    weight_src <- total_weights[src]
    size_src <- graph$size[src]
    maxgain <- 0.0
    for(i in connected){
      if(i != src){
        gain <- total_weights[i] + weight_u - weight_src - γ * (graph$size[i] - size_src + c_u) * c_u
        if(gain > maxgain){
          dst <- i
          maxgain <- gain
        }
        total_weights[i] = 0
      }
    }
    total_weights[src] = 0
    #@assert all(x == 0 for x in total_weights)

    if(src != dst){
      # move `u` to the best community and add its neighbors to the queue if needed
      graph <- move_node(graph, u, dst)
      for(v in neighbors(graph, u)){
        if((graph$membership[v] != graph$membership[u]) && !(v %in% queued)){
          queue <- c(queue, v)
          queued <- c(queued, v)
        }
      }
    }
    #print(queued)
  }
  if(length(queued) == 0){# | is.na(queued[1])){
    print("passes")
  } else{
    warning("error")
  }
  return(drop_empty_communities(graph))
}


##' @param graph PartitionedGraph Object
##' @param  γ numeric
##' @param θ numeric

refine_partition <- function(graph, γ, θ){
  if(!(γ > 0)) warning("assert γ > 0")
  if(!(θ > 0)) warning("assert θ > 0")
  refined <- PartitionedGraph$new(graph$edge_weight, cardinality = graph$cardinality)
  is_well_connected_node <- function (u){ #integer
    u <- as.integer(u)
    i <- graph$membership[u]
    c <- graph$cardinality[u]
    threshold <- γ * c * (graph$size[i] - c)
    x <- 0.0
    for(v in graph$partition[[i]]){
      if(v != u){
        x <- x + graph$edge_weight[u,v]
        if(x >= threshold){ # return as early as possible
          return(TRUE)
        }
      } else {
        print("computing well-connected graph")
      }
    }
    return(FALSE)
  }
  is_well_connected_community <- function(u, i, between_weights){
    #u::Int, i::Int, between_weights::Vector{Float64}
    sz <- refined$size[i]
    return(between_weights[i] >= γ * sz * (graph$size[graph$membership[u]] - sz))
  }
  is_singleton <- function(u){
    #u::Int
    u <- as.integer(u)
    return(length(refined$partition[refined$membership[u]]) == 1)
  }
  total_weights <- rep(0, nc(refined))
  between_weights <- rep(0, nc(refined))
  for(subset in graph$partition){
    communities <- as.integer(c())
    logprobs <- as.numeric(c())
    indexes <- as.integer(c())
    for(u in subset){
      weight <- 0.0
      for(v in subset){
        if(v != u){
          weight <- weight + refined$edge_weight[u,v]
        }
      }
      between_weights[refined$membership[u]] <- weight
    }
    for(u in shuffle(subset)){
      if((!is_well_connected_node(u) || !is_singleton(u))){

        communities <- as.integer(c())
        for(v in neighbors(refined, u)){
          if(v %in% subset){
            i <- refined$membership[v]
            if(total_weights[[i]] == 0){
              communities <- c(communities, i)
            }
            total_weights[[i]] <- total_weights[[i]] + refined$edge_weight[u,v]
          }
        }
      }

      c_u <- refined$cardinality[u]
      weight_u <- refined$edge_weight[u,u]
      src <- refined$membership[u]
      weight_src <- total_weights[src]
      size_src <- refined$size[src]
      logprobs <- as.numeric(c())
      indexes <- as.integer(c())
      for(i in communities){
        if(!(i == src || !is_well_connected_community(u, i, between_weights))){
          gain <- total_weights[[i]] + weight_u - weight_src - γ * (refined$size[i] - size_src + c_u) * c_u
          if(gain >= 0){
            logprobs <- c(logprobs, 1/θ * gain)
            indexes <- c(indexes, i)
          }
        }
      }
      total_weights[communities] <- 0
      if(!(length(indexes) == 0)){# | is.na(indexes))){

        probs <- exp(logprobs - logsumexp(logprobs))
        dst <- indexes[sample(prob = probs)]
        refined <- move_node(refined, u, dst)

        for(v in neighbors(refined, u)){
          if(!(v == u || !(v %in% subset))){
            i <- refined$membership[v]
            weight <- refined$edge_weight[u,v]
            if(i == dst){
              between_weights[src] <- between_weights[src] - weight
              between_weights[dst] <- between_weights[dst] - weight
            } else if(i == src){
              between_weights[src] <- between_weights[src] + weight
              between_weights[dst] <- between_weights[dst] + weight
            } else {
              between_weights[src] <- between_weights[src] -  weight
              between_weights[dst] <- between_weights[dst] + weight
            }
          }
        }
      }
      for(u in subset){
        between_weights[refined$membership[u]] <- 0
      }
      #@assert all(x == 0 for x in between_weights)
    }
    return(drop_empty_communities(refined))
  }
}

# function sample(probs::Vector{Float64})
# r = rand()
# p = 0.0
# i = 1
# while i < lastindex(probs)
# p += probs[i]
# if p > r
# return i
# end
# i += 1
# end
# return lastindex(probs)
# end

##' @param graph::PartitionedGraph
aggregate_graph <- function(graph){
  n <- nv(graph)
  I <- as.integer(c())
  J <- as.integer(c())
  V <- as.numeric(c())
  cardinality <- rep(1L, n)
  connected <- as.integer(c())
  total_weights <- rep(0.0, n)
  for(i in 1:length(graph$partition)){
    community <- graph$partition[[i]]
    cardinality[i] = graph$size[i]
    connected <- as.integer(c())
    for(u in community){
      for(v in neighbors(graph, u)){
        j <- graph$membership[v]
        if(!(i == j && u > v)){

          if(total_weights[j] == 0){
            connected <- c(connected, j)
          }
          total_weights[j] <- total_weights[j] + graph$edge_weight[u,v]
        }
      }
    }
    for(j in connected){
      v <- total_weights[j]
      I <- c(I, i)
      J <- c(J, j)
      V <- c(V, v)
      total_weights[j] <- 0
    }
    #@assert all(x == 0 for x in total_weights)
  }
  #print(I)
  #print(J)
  adj <- sparseMatrix(I, J, x = V, dims = c(n, n))
  graph <- PartitionedGraph$new(adj, cardinality = cardinality, partition = graph$partition)
  #graph$size <- sapply(graph$partition, length)
  return(graph)
}

#graph::PartitionedGraph, γ::Float64
H <- function (graph, γ){
  quality <- 0.0
  if(length(graph$partition) > 0){
    for(i in 1:length(graph$partition)){
      community <- graph$partition[[i]]
      #print(community)
      for(u in community){
        #print(u)
        for(v in community){
          #print(v)
          if(u <= v){
            quality <- quality + graph$edge_weight[u,v]
          }
        }
      }
      n <- graph$size[i]
      quality <- quality - γ * as.integer(round(n * (n - 1) / 2, 0))
    }
  }
  return(quality)
}


#stack::Vector{Partition}
flatten <- function(stack){
  k <- length(stack)
  result <- stack[[k]]
  result <- result[order(sapply(result, length))]
  result <- lapply(result, sort)
  names(result) <- 1:length(result)
  return(result)
}

# end # module

# ##' @importFrom igraph V E
# nv <- function(graph){
#   length(V(graph))
# }
# nc <- function(graph){
#   length(E(graph))
# }
