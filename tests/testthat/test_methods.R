library("leiden")
library("reticulate")
library("igraph")
context("running Leiden on an dense and sparse matrices")

adj_mat <- matrix(round(runif(10000, 0, 1)), 100, 100)
snn_graph <- graph_from_adjacency_matrix(adj_mat)

modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")

skip_if_no_python <- function() {
  if (!modules)
    testthat::skip("leidenalg not available for testing")
}

test_that("run dense matrix with defaults", {
  skip_if_no_python()
  partition <- leiden(adj_mat)
  expect_length(partition, 100)
})

test_that("run graph object with defaults", {
  skip_if_no_python()
  partition <- leiden(snn_graph)
  expect_length(partition, 100)
})

library("Matrix")
adj_mat_sparse <- as(adj_mat, "Matrix")

test_that("run dsCMatrix with defaults", {
  skip_if_no_python()
  partition <- leiden(adj_mat_sparse)
  expect_length(partition, 100)
})


adj_mat_sparse <- as(adj_mat, "dgCMatrix")

test_that("run sparse dgCMatrix matrix with defaults", {
  skip_if_no_python()
  partition <- leiden(adj_mat_sparse)
  expect_length(partition, 100)
})

adj_mat_sparse <- as(adj_mat, "dgeMatrix")

test_that("run sparse dgeMatrix matrix with defaults", {
  skip_if_no_python()
  partition <- leiden(adj_mat_sparse)
  expect_length(partition, 100)
})

adj_mat_df<- as.data.frame(adj_mat)

test_that("run data.frame matrix with defaults", {
  skip_if_no_python()
  partition <- leiden(adj_mat_df)
  expect_length(partition, 100)
})

library("data.table")
adj_mat_dt <- as.data.table(adj_mat)

test_that("run data.table matrix with defaults", {
  skip_if_no_python()
  partition <- leiden(adj_mat_dt)
  expect_length(partition, 100)
})

library("tibble")
adj_mat_tb <- as_tibble(adj_mat, .name_repair = "unique")

test_that("run data.table matrix with defaults", {
  skip_if_no_python()
  partition <- leiden(adj_mat_tb)
  expect_length(partition, 100)
})


