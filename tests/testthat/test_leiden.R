library("leiden")
context("running Leiden on an adjacency matrix")

adj_mat <- matrix(round(runif(10000, 0, 1)), 100, 100)

test_that("run with defaults", {
  partition <- leiden(adj_mat)
  expect_length(partition, 100)
})

test_that("run with ModularityVertexPartition", {
  partition <- leiden(adj_mat, partition_type = "ModularityVertexPartition")
  expect_length(partition, 100)
})

test_that("run with resolution parameter", {
  partition <- leiden(adj_mat, resolution_parameter = 0.95)
  expect_length(partition, 100)
})
