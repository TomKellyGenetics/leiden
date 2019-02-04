## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library("reticulate")
module <- py_available()

## ---- eval=FALSE---------------------------------------------------------
#  if (!requireNamespace("devtools"))
#      install.packages("devtools")
#  devtools::install_github("TomKellyGenetics/leiden")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("leiden")

## ------------------------------------------------------------------------
library("leiden")

## ------------------------------------------------------------------------
adjacency_matrix <- rbind(cbind(matrix(round(rbinom(4000, 1, 0.8)), 20, 20),
                                matrix(round(rbinom(4000, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.1)), 20, 20)),
                          cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(400, 1, 0.8)), 20, 20), 
                                matrix(round(rbinom(4000, 1, 0.2)), 20, 20)),
                          cbind(matrix(round(rbinom(400, 1, 0.3)), 20, 20), 
                                matrix(round(rbinom(4000, 1, 0.1)), 20, 20), 
                                matrix(round(rbinom(4000, 1, 0.9)), 20, 20)))
str(adjacency_matrix)
dim(adjacency_matrix )

## ------------------------------------------------------------------------
library("igraph")
rownames(adjacency_matrix) <- 1:60
colnames(adjacency_matrix) <- 1:60
graph_object <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
graph_object

## ---- warning=FALSE, message=FALSE---------------------------------------
plot(graph_object, vertex.color = "grey75")

## ---- eval=FALSE---------------------------------------------------------
#  library("igraph")
#  adjacency_matrix <- igraph::as_adjacency_matrix(graph_object)

## ------------------------------------------------------------------------
partition <- leiden(adjacency_matrix)
table(partition)

## ---- warning=FALSE, message=FALSE---------------------------------------
library("RColorBrewer")
node.cols <- brewer.pal(max(partition),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols)

## ---- eval=module--------------------------------------------------------
#  #run with defaults
#    partition <- leiden(adjacency_matrix)
#  
#  
#  #run with ModularityVertexPartition"
#    partition <- leiden(adjacency_matrix, partition_type = "ModularityVertexPartition")
#  
#  
#  #run with resolution parameter
#    partition <- leiden(adjacency_matrix, resolution_parameter = 0.95)

## ---- warning=FALSE, message=FALSE---------------------------------------
partition <- leiden(adjacency_matrix, resolution_parameter = 0.5)
node.cols <- brewer.pal(max(partition),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols)

## ---- warning=FALSE, message=FALSE---------------------------------------
partition <- leiden(adjacency_matrix, resolution_parameter = 1.8)
node.cols <- brewer.pal(max(partition),"Pastel1")[partition]
plot(graph_object, vertex.color = node.cols)

## ---- eval=FALSE---------------------------------------------------------
#  adjacency_matrix <- as.matrix(object@snn)
#  membership <- leiden(adjacency_matrix)
#  object@ident <- as.factor(membership)
#  names(test@ident) <- rownames(test@meta.data)
#  object@meta.data$ident <- as.factor(membership)

