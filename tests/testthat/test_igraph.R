library("leiden")
library("reticulate")
library("igraph")
context("running Leiden on an igraph object")

adj_mat <- matrix(round(runif(10000, 0, 1)), 100, 100)
snn_graph <- graph_from_adjacency_matrix(adj_mat)

modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")

skip_if_no_python <- function() {
  if (!modules)
    testthat::skip("scipy not available for testing")
}


test_that("run with defaults", {
  skip_if_no_python()
  partition <- leiden(snn_graph)
  expect_length(partition, 100)
})

test_that("run with ModularityVertexPartition", {
  skip_if_no_python()
  partition <- leiden(snn_graph, partition_type = "ModularityVertexPartition")
  expect_length(partition, 100)
})

test_that("run with resolution parameter", {
  skip_if_no_python()
  partition <- leiden(snn_graph, resolution_parameter = 0.95)
  expect_length(partition, 100)
})

weights <- sample(1:10, sum(adj_mat!=0), replace=TRUE)

test_that("run with non-wieghted adjacency matrix and weights vector", {
  skip_if_no_python()
  partition <- leiden(snn_graph, weights = weights)
  expect_length(partition, 100)
})

adj_mat <- ifelse(adj_mat == 1, weights, 0)

test_that("run with wieghted adjacency matrix", {
  skip_if_no_python()
  partition <- leiden(snn_graph)
  expect_length(partition, 100)
})

rownames(adj_mat) <- 1:nrow(adj_mat)
colnames(adj_mat) <- 1:ncol(adj_mat)

test_that("run with named adjacency matrix", {
  skip_if_no_python()
  partition <- leiden(snn_graph)
  expect_length(partition, 100)
})

