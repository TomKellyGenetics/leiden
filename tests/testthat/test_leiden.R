library("leiden")
library("reticulate")
context("running Leiden on an adjacency matrix")

adj_mat <- matrix(round(runif(10000, 0, 1)), 100, 100)

modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")

skip_if_no_python <- function() {
  if (!modules)
    testthat::skip("scipy not available for testing")
}


test_that("run with defaults", {
  skip_if_no_python()
  partition <- leiden(adj_mat)
  expect_length(partition, 100)
})

test_that("run with ModularityVertexPartition", {
  skip_if_no_python()
  partition <- leiden(adj_mat, partition_type = "ModularityVertexPartition")
  expect_length(partition, 100)
})

test_that("run with resolution parameter", {
  skip_if_no_python()
  partition <- leiden(adj_mat, resolution_parameter = 0.95)
  expect_length(partition, 100)
})

weights <- sample(1:10, sum(adj_mat!=0), replace=TRUE)

test_that("run with non-wieghted adjacency matrix and weights vector", {
  skip_if_no_python()
  partition <- leiden(adj_mat, weights = weights)
  expect_length(partition, 100)
})

adj_mat[adj_mat == 1] <- weights

test_that("run with wieghted adjacency matrix", {
  skip_if_no_python()
  partition <- leiden(adj_mat)
  expect_length(partition, 100)
})
