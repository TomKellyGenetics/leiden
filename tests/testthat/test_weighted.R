library("leiden")
library("reticulate")
library("igraph")
context("running Leiden on unweighted objects")

modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")

skip_if_no_python <- function() {
  if (!modules)
    testthat::skip("leidenalg not available for testing")
}

library("leiden")
library("igraph")
set.seed(9000)
mat <- round(matrix(runif(100, max = 1), 10, 10, ), 0)
sp.mat <- as(mat, Class = "dgCMatrix")
graph <- igraph::graph_from_adjacency_matrix(mat, weighted = NULL)


test_that("run with unweighted dense matrix", {
  skip_if_no_python()
  set.seed(9000)
  mat <- round(matrix(runif(100, max = 1), 10, 10, ), 0)
  mat
  part_mat0 <- leiden::leiden(mat, seed = 9000)
  expect_equal(part_mat0, c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1))
})

test_that("run with unweighted sparse matrix", {
  skip_if_no_python()
  sp.mat <- as(mat, Class = "dgCMatrix")
  part_sp.mat0 <- leiden::leiden(sp.mat, seed = 9000)
  expect_equal(part_sp.mat0, c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1))
})


test_that("run with unweighted graph object", {
  skip_if_no_python()
  graph <- igraph::graph_from_adjacency_matrix(mat, weighted = NULL)
  graph
  part_graph0 <- leiden::leiden(graph, seed = 9000)
  part_graph0
  expect_equal(part_graph0,  c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1))
})

test_that("same output with different input class", {
  skip_if_no_python()
  part_mat0 <- leiden::leiden(mat, seed = 9000)
  part_graph0 <- leiden::leiden(graph, seed = 9000)
  expect_equivalent(part_mat0, c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1))
  expect_equivalent(part_graph0, c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1))
  expect_equivalent(part_mat0, part_graph0)
  expect_equal(sum(diag(table(part_mat0, part_graph0))), length(part_mat0))
  expect_equal(sum(diag(table(part_mat0, part_graph0))), nrow(mat))
  expect_equal(sum(diag(table(part_mat0, part_graph0))), length(V(graph)))
})

# # Test Seurat (development version)
# devtools::install_github("TomKellyGenetics/seurat", ref = "pr", quiet = TRUE)
# devtools::install_github("satijalab/seurat", ref = "develop", quiet = TRUE)
# library("Seurat")
# rownames(mat) <- colnames(mat) <- 1:10
# obj <- Seurat::CreateSeuratObject(mat)
# obj <- Seurat::FindVariableFeatures(obj)
# obj <- Seurat::ScaleData(obj)
# obj <- Seurat::RunPCA(obj, dim = 1:3)
# obj <- Seurat::FindNeighbors(obj, dim = 1:3,)
# obj@graphs$RNA_nn <- Seurat:::as.Graph.matrix(mat)
# obj@graphs$RNA_snn <- Seurat:::as.Graph.matrix(mat)
# obj <- Seurat::FindClusters(obj, algorithm = "Leiden", method = "matrix", random.seed = 9000, weights =  mat[mat !=0], resolution = 1)
# part_mat <- obj@active.ident
# obj <- Seurat::FindClusters(obj, algorithm = "Leiden", method = "igraph", random.seed = 9000, weights =  mat[mat !=0], resolution = 1)
# part_graph <- obj@active.ident
# table(part_mat, part_graph)
# all(part_mat == c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1))
# all(part_graph == c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1))
# sum(diag(table(part_mat, part_graph))) == length(part_mat)
# sum(diag(table(part_mat, part_graph))) == nrow(mat)
# sum(diag(table(part_mat, part_graph))) == length(V(graph))


context("running Leiden on weighted objects")

set.seed(9000)
mat <- round(matrix(runif(100, max = 3), 10, 10), 0)
sp.mat <- as(mat, Class = "dgCMatrix")
graph <- igraph::graph_from_adjacency_matrix(mat, weighted = NULL)

test_that("run with weighted dense matrix", {
  skip_if_no_python()
  set.seed(9000)
  mat <- round(matrix(runif(100, max = 3), 10, 10, ), 0)
  mat
  #passing weights
  part_mat1 <- leiden::leiden(mat, seed = 9000, weights = t(mat)[t(mat) != 0])
  expect_equal(part_mat1, c(2, 1, 1, 1, 2, 1, 1, 2, 2, 2))
  #detecting weights from matrix
  part_mat1 <- leiden::leiden(mat, seed = 9000)
  expect_equal(part_mat1, c(2, 1, 1, 1, 2, 1, 1, 2, 2, 2))
})

test_that("run with weighted sparse matrix", {
  skip_if_no_python()
  sp.mat <- as(mat, Class = "dgCMatrix")
  #passing weights
  part_sp.mat1 <- leiden::leiden(sp.mat, seed = 9000, weights = mat[mat != 0])
  expect_equal(part_sp.mat1, c(2, 1, 1, 1, 2, 1, 1, 2, 2, 2))
  #detecting weights from matrix
  part_sp.mat1 <- leiden::leiden(sp.mat, seed = 9000)
  expect_equal(part_sp.mat1, c(2, 1, 1, 1, 2, 1, 1, 2, 2, 2))
})


test_that("run with weighted graph object", {
  skip_if_no_python()
  graph <- igraph::graph_from_adjacency_matrix(mat, weighted = NULL)
  graph
  #passing weights
  part_graph1 <- leiden::leiden(graph, seed = 9000, weights = edge.attributes(graph)$weight)
  expect_equal(part_graph1, c(2, 1, 1, 1, 2, 1, 1, 2, 2, 2))
  #detecting weights from matrix
  part_graph1 <- leiden::leiden(graph, seed = 9000)
  expect_equal(part_graph1, c(2, 1, 1, 1, 2, 1, 1, 2, 2, 2))
})

test_that("same output with different input class", {
  skip_if_no_python()
  part_mat1 <- leiden::leiden(mat, seed = 9000)
  part_graph1 <- leiden::leiden(graph, seed = 9000)
  expect_equivalent(part_mat1, c(2, 1, 1, 1, 2, 1, 1, 2, 2, 2))
  expect_equivalent(part_graph1, c(2, 1, 1, 1, 2, 1, 1, 2, 2, 2))
  expect_equivalent(part_mat1, part_graph1)
  expect_equal(sum(diag(table(part_mat1, part_graph1))), length(part_mat1))
  expect_equal(sum(diag(table(part_mat1, part_graph1))), nrow(mat))
  expect_equal(sum(diag(table(part_mat1, part_graph1))), length(V(graph)))
})

# # Test Seurat (development version)
# devtools::install_github("TomKellyGenetics/seurat", ref = "pr", quiet = TRUE)
# devtools::install_github("satijalab/seurat", ref = "develop", quiet = TRUE)
# library("Seurat")
# rownames(mat) <- colnames(mat) <-1:10
# obj <- Seurat::CreateSeuratObject(mat)
# obj <- Seurat::FindVariableFeatures(obj)
# obj <- Seurat::ScaleData(obj)
# obj <- Seurat::RunPCA(obj, dim = 1:3)
# obj <- Seurat::FindNeighbors(obj, dim = 1:3,)
# obj@graphs$RNA_snn <- Seurat:::as.Graph.matrix(mat)
# weights <- mat[mat !=0]
# obj <- Seurat::FindClusters(obj, algorithm = "Leiden", method = "matrix", random.seed = 9000, weights = weights, resolution = 1)
# part_mat2 <- obj@active.ident
# obj <- Seurat::FindClusters(obj, algorithm = "Leiden", method = "igraph",  random.seed = 9000, weights = weights, resolution = 1)
# part_graph2 <- obj@active.ident
# table(part_mat, part_graph)
# all(part_mat2 == c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1))
# all(part_graph2 == c(1, 2, 2, 2, 1, 2, 2, 1, 1, 1))
# sum(diag(table(part_mat2, part_graph2))) == length(part_mat)
# sum(diag(table(part_mat2, part_graph2))) == nrow(mat)
# sum(diag(table(part_mat2, part_graph2))) == length(V(graph))

