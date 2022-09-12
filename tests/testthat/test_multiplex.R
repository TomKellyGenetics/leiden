library("reticulate")
library("igraph")
library("leiden")
library("multiplex")
set.seed(9000)

context("running Leiden on multiplex igraph objects")

suppressWarnings(suppressMessages({
  #imported from Achab94/mplex::aarhus_mplex
  multiplex_mplex <-  list(nodes = structure(list(nodeID = 1:61, nodeLabel = structure(c(3L,
                                                                       18L, 33L, 4L, 5L, 10L, 11L, 1L, 26L, 27L, 29L, 30L, 32L, 36L,
                                                                       37L, 40L, 43L, 60L, 12L, 23L, 51L, 56L, 58L, 6L, 13L, 15L, 16L,
                                                                       24L, 31L, 38L, 42L, 45L, 54L, 55L, 57L, 61L, 2L, 14L, 22L, 19L,
                                                                       25L, 28L, 34L, 35L, 53L, 7L, 9L, 17L, 41L, 47L, 48L, 52L, 8L,
                                                                       39L, 49L, 50L, 46L, 44L, 59L, 20L, 21L),
                                                                     .Label = c("U1", "U10",
                                                                                "U102", "U106", "U107", "U109", "U110", "U112", "U113", "U118",
                                                                                "U123", "U124", "U126", "U13", "U130", "U134", "U138", "U139",
                                                                                "U14", "U140", "U141", "U142", "U17", "U18", "U19", "U21", "U22",
                                                                                "U23", "U26", "U29", "U3", "U32", "U33", "U37", "U4", "U41",
                                                                                "U42", "U47", "U48", "U49", "U53", "U54", "U59", "U6", "U62",
                                                                                "U63", "U65", "U67", "U68", "U69", "U71", "U72", "U73", "U76",
                                                                                "U79", "U86", "U90", "U91", "U92", "U97", "U99"), class = "factor")),
                           class = "data.frame",
                           row.names = c(NA, -61L)), layerNames = c("lunch", "facebook", "coauthor", "leisure",  "work"),
         L1 = structure(list(ID_Start = c(1L, 1L, 10L, 10L, 10L,
                                          10L, 11L, 11L, 11L, 12L, 12L, 12L, 13L, 13L, 13L, 13L, 14L, 15L,
                                          15L, 17L, 17L, 17L, 17L, 17L, 17L, 17L, 18L, 18L, 18L, 18L, 18L,
                                          18L, 19L, 19L, 19L, 19L, 20L, 20L, 20L, 21L, 21L, 21L, 22L, 23L,
                                          23L, 23L, 23L, 23L, 23L, 24L, 24L, 24L, 24L, 24L, 24L, 24L, 24L,
                                          24L, 24L, 24L, 25L, 25L, 25L, 25L, 25L, 25L, 25L, 25L, 25L, 26L,
                                          26L, 26L, 26L, 27L, 27L, 27L, 27L, 27L, 27L, 27L, 27L, 27L, 27L,
                                          28L, 28L, 28L, 28L, 28L, 28L, 28L, 29L, 29L, 29L, 29L, 29L, 3L,
                                          3L, 3L, 3L, 3L, 30L, 30L, 31L, 31L, 31L, 31L, 32L, 32L, 32L,
                                          33L, 33L, 33L, 34L, 34L, 34L, 35L, 37L, 37L, 38L, 38L, 39L, 4L,
                                          4L, 4L, 4L, 4L, 4L, 4L, 40L, 40L, 40L, 41L, 41L, 42L, 43L, 43L,
                                          43L, 44L, 44L, 44L, 44L, 44L, 44L, 44L, 46L, 46L, 46L, 47L, 47L,
                                          47L, 48L, 49L, 49L, 5L, 5L, 5L, 50L, 50L, 51L, 51L, 52L, 53L,
                                          53L, 54L, 54L, 54L, 55L, 55L, 56L, 59L, 6L, 6L, 6L, 6L, 6L, 7L,
                                          7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 9L,
                                          9L), ID_Arrive = c(2L, 3L, 11L, 14L, 15L, 16L, 14L, 15L, 16L,
                                                             13L, 20L, 40L, 20L, 26L, 40L, 45L, 15L, 16L, 30L, 23L, 25L, 46L,
                                                             49L, 50L, 51L, 52L, 21L, 25L, 43L, 44L, 56L, 58L, 24L, 25L, 30L,
                                                             36L, 40L, 42L, 45L, 44L, 51L, 57L, 26L, 25L, 46L, 49L, 50L, 51L,
                                                             52L, 25L, 27L, 28L, 29L, 30L, 31L, 32L, 33L, 34L, 35L, 36L, 31L,
                                                             35L, 43L, 46L, 47L, 48L, 52L, 56L, 58L, 27L, 34L, 44L, 46L, 31L,
                                                             32L, 33L, 34L, 35L, 36L, 38L, 44L, 54L, 59L, 29L, 31L, 32L, 33L,
                                                             34L, 35L, 36L, 31L, 32L, 33L, 34L, 35L, 21L, 44L, 51L, 57L, 7L,
                                                             34L, 36L, 32L, 33L, 34L, 35L, 33L, 34L, 35L, 34L, 35L, 36L, 35L,
                                                             36L, 57L, 58L, 40L, 41L, 39L, 44L, 44L, 10L, 11L, 14L, 15L, 16L,
                                                             18L, 6L, 41L, 42L, 45L, 42L, 45L, 45L, 51L, 56L, 58L, 51L, 53L,
                                                             54L, 55L, 57L, 59L, 61L, 47L, 49L, 52L, 48L, 51L, 52L, 52L, 50L,
                                                             52L, 12L, 13L, 20L, 51L, 52L, 52L, 57L, 56L, 54L, 55L, 55L, 59L,
                                                             61L, 59L, 61L, 58L, 61L, 10L, 11L, 14L, 15L, 16L, 21L, 34L, 44L,
                                                             51L, 57L, 12L, 13L, 37L, 40L, 41L, 42L, 45L, 18L, 25L, 43L, 51L,
                                                             56L, 58L), Weight = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                                                                   1L, 1L, 1L, 1L, 1L, 1L, 1L)), row.names = c(NA, 193L), class = "data.frame"),
         L2 = structure(list(ID_Start = c(12L, 12L, 13L, 13L, 13L,
                                          13L, 13L, 13L, 15L, 15L, 15L, 17L, 17L, 17L, 17L, 19L, 19L,
                                          19L, 19L, 19L, 19L, 19L, 19L, 19L, 19L, 19L, 19L, 21L, 21L,
                                          21L, 21L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 24L,
                                          24L, 24L, 24L, 24L, 26L, 26L, 26L, 26L, 26L, 26L, 26L, 26L,
                                          26L, 28L, 28L, 28L, 28L, 28L, 29L, 29L, 29L, 29L, 30L, 30L,
                                          30L, 30L, 30L, 31L, 31L, 31L, 31L, 31L, 33L, 33L, 33L, 34L,
                                          34L, 34L, 37L, 37L, 39L, 39L, 39L, 4L, 4L, 4L, 4L, 4L, 4L,
                                          4L, 44L, 44L, 46L, 46L, 47L, 47L, 47L, 5L, 5L, 5L, 5L, 50L,
                                          51L, 51L, 56L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 8L,
                                          8L, 8L, 8L, 8L, 9L, 9L, 9L),
                             ID_Arrive = c(13L, 23L, 21L,
                                           23L, 26L, 44L, 46L, 51L, 16L, 34L, 39L, 23L, 46L, 47L, 51L,
                                           23L, 24L, 26L, 27L, 28L, 29L, 30L, 31L, 33L, 34L, 56L, 58L,
                                           34L, 44L, 46L, 51L, 34L, 37L, 39L, 44L, 46L, 47L, 50L, 51L,
                                           56L, 28L, 30L, 31L, 33L, 34L, 27L, 28L, 29L, 30L, 33L, 34L,
                                           39L, 44L, 51L, 29L, 30L, 31L, 33L, 34L, 30L, 31L, 33L, 34L,
                                           31L, 33L, 34L, 39L, 44L, 33L, 34L, 37L, 39L, 44L, 34L, 44L,
                                           53L, 39L, 46L, 50L, 39L, 44L, 44L, 46L, 51L, 12L, 13L, 17L,
                                           5L, 7L, 8L, 9L, 46L, 51L, 47L, 51L, 50L, 51L, 56L, 12L, 13L,
                                           19L, 21L, 53L, 56L, 58L, 58L, 13L, 15L, 16L, 21L, 23L, 26L,
                                           30L, 39L, 44L, 51L, 12L, 13L, 21L, 34L, 37L, 51L, 56L, 58L
                             ),
                             Weight = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L)), row.names = 194:317, class = "data.frame"),
         L3 = structure(list(ID_Start = c(10L, 12L, 18L, 23L, 23L,
                                          23L, 26L, 26L, 26L, 26L, 26L, 28L, 30L, 38L, 39L, 4L, 46L,
                                          46L, 49L, 6L, 8L),
                             ID_Arrive = c(11L, 13L, 46L, 46L, 49L,
                                           52L, 27L, 28L, 30L, 33L, 36L, 33L, 36L, 54L, 55L, 6L, 48L,
                                           49L, 52L, 14L, 37L),
                             Weight = c(1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L)), row.names = 318:338, class = "data.frame"),
         L4 = structure(list(ID_Start = c(10L, 12L, 12L, 15L, 15L,
                                          15L, 17L, 17L, 17L, 17L, 19L, 19L, 19L, 19L, 20L, 20L, 20L,
                                          20L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 23L, 24L,
                                          24L, 24L, 24L, 25L, 25L, 25L, 25L, 28L, 28L, 28L, 29L, 29L,
                                          30L, 31L, 31L, 31L, 31L, 32L, 33L, 33L, 34L, 34L, 34L, 34L,
                                          35L, 37L, 37L, 37L, 37L, 39L, 39L, 4L, 4L, 40L, 40L, 40L,
                                          41L, 41L, 42L, 44L, 46L, 46L, 46L, 46L, 47L, 47L, 48L, 5L,
                                          5L, 5L, 50L, 55L, 6L, 8L, 8L, 8L, 8L, 8L, 9L),
                             ID_Arrive = c(15L,
                                           13L, 25L, 16L, 34L, 39L, 23L, 46L, 50L, 52L, 23L, 25L, 28L,
                                           36L, 23L, 40L, 42L, 45L, 25L, 31L, 35L, 46L, 47L, 48L, 49L,
                                           50L, 52L, 56L, 25L, 31L, 33L, 35L, 31L, 35L, 46L, 56L, 32L,
                                           33L, 36L, 31L, 35L, 33L, 33L, 34L, 35L, 37L, 33L, 34L, 35L,
                                           35L, 36L, 45L, 50L, 58L, 38L, 39L, 43L, 45L, 43L, 55L, 14L,
                                           6L, 41L, 42L, 45L, 42L, 45L, 45L, 51L, 47L, 48L, 50L, 51L,
                                           48L, 50L, 52L, 13L, 20L, 23L, 52L, 61L, 14L, 11L, 37L, 40L,
                                           42L, 45L, 25L),
                             Weight = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L)), row.names = 339:426, class = "data.frame"),
         L5 = structure(list(ID_Start = c(10L, 10L, 10L, 10L, 11L,
                                          11L, 11L, 11L, 11L, 11L, 11L, 11L, 11L, 12L, 12L, 13L, 13L,
                                          13L, 13L, 13L, 13L, 13L, 15L, 17L, 17L, 17L, 18L, 18L, 18L,
                                          19L, 19L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 20L,
                                          20L, 20L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 21L, 23L, 23L,
                                          23L, 23L, 23L, 23L, 23L, 24L, 24L, 24L, 25L, 25L, 25L, 26L,
                                          26L, 26L, 26L, 26L, 26L, 26L, 26L, 26L, 26L, 27L, 27L, 28L,
                                          29L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 30L, 30L, 31L, 31L,
                                          31L, 31L, 31L, 32L, 33L, 34L, 34L, 34L, 35L, 36L, 37L, 37L,
                                          37L, 37L, 38L, 38L, 38L, 39L, 39L, 4L, 4L, 4L, 4L, 40L, 40L,
                                          40L, 41L, 41L, 42L, 43L, 44L, 44L, 44L, 44L, 44L, 44L, 44L,
                                          46L, 46L, 46L, 46L, 46L, 47L, 47L, 48L, 5L, 5L, 5L, 5L, 5L,
                                          50L, 50L, 51L, 51L, 51L, 51L, 53L, 54L, 55L, 6L, 6L, 6L,
                                          6L, 6L, 6L, 6L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L,
                                          7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L,
                                          8L, 8L, 8L, 8L, 8L, 8L, 8L, 8L, 9L),
                             ID_Arrive = c(11L, 14L,
                                           15L, 16L, 13L, 15L, 16L, 18L, 21L, 26L, 34L, 46L, 60L, 13L,
                                           20L, 18L, 20L, 22L, 26L, 40L, 45L, 46L, 34L, 23L, 46L, 51L,
                                           21L, 44L, 51L, 26L, 44L, 11L, 12L, 18L, 21L, 37L, 40L, 41L,
                                           5L, 51L, 7L, 8L, 40L, 42L, 45L, 40L, 41L, 44L, 46L, 51L,
                                           57L, 58L, 60L, 46L, 47L, 48L, 49L, 50L, 51L, 52L, 26L, 31L,
                                           32L, 46L, 51L, 56L, 27L, 28L, 30L, 31L, 32L, 33L, 34L, 36L,
                                           37L, 44L, 36L, 44L, 32L, 35L, 11L, 21L, 44L, 46L, 51L, 57L,
                                           6L, 7L, 9L, 32L, 36L, 34L, 35L, 38L, 41L, 44L, 33L, 44L,
                                           36L, 44L, 46L, 44L, 44L, 40L, 41L, 42L, 45L, 44L, 54L, 59L,
                                           44L, 55L, 11L, 14L, 6L, 7L, 41L, 42L, 45L, 42L, 45L, 45L,
                                           51L, 51L, 53L, 54L, 55L, 57L, 59L, 61L, 47L, 48L, 49L, 51L,
                                           52L, 48L, 51L, 51L, 12L, 13L, 20L, 21L, 22L, 51L, 52L, 52L,
                                           56L, 57L, 58L, 55L, 56L, 61L, 11L, 14L, 18L, 21L, 35L, 51L,
                                           7L, 10L, 11L, 15L, 16L, 18L, 19L, 20L, 21L, 26L, 28L, 29L,
                                           30L, 31L, 32L, 33L, 34L, 35L, 36L, 39L, 44L, 46L, 51L, 57L,
                                           11L, 13L, 19L, 21L, 26L, 34L, 37L, 40L, 41L, 42L, 45L, 51L
                             ),
                             Weight = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                        1L, 1L, 1L)), row.names = 427:620, class = "data.frame"))

  lunch <- graph_from_edgelist(as.matrix(multiplex_mplex[[3]][,1:2]))
  V(lunch)$name <- as.character(multiplex_mplex[[1]][,2])
  facebook <- graph_from_edgelist(as.matrix(multiplex_mplex[[4]][,1:2]))
  facebook <- graph_from_edgelist(
    rbind(as.matrix(multiplex_mplex[[4]][,1:2]),
          cbind(setdiff(as.integer(V(lunch)), as.integer(V(facebook))),
                           setdiff(as.integer(V(lunch)), as.integer(V(facebook))))
    )
  )
  V(facebook)$name <- as.character(multiplex_mplex[[1]][,2])
  coauthor <- graph_from_edgelist(as.matrix(multiplex_mplex[[5]][,1:2]))
  coauthor <- graph_from_edgelist(
    rbind(as_edgelist(coauthor),
          cbind(setdiff(as.integer(V(lunch)), as.integer(V(coauthor))),
                setdiff(as.integer(V(lunch)), as.integer(V(coauthor))))
    )
  )
  leisure <- graph_from_edgelist(as.matrix(multiplex_mplex[[6]][,1:2]))
  V(leisure)$name <- as.character(multiplex_mplex[[1]][,2])
  work <- graph_from_edgelist(as.matrix(multiplex_mplex[[7]][,1:2]))
  V(work)$name <- as.character(multiplex_mplex[[1]][,2])

  multiplex_graph <- list(lunch, facebook, coauthor, leisure, work)
  multiplex_graph <- lapply(multiplex_graph, upgrade_graph)
}))

multiplex_adj_mat <- lapply(multiplex_graph, igraph::as_adjacency_matrix)

modules <- reticulate::py_module_available("leidenalg") && reticulate::py_module_available("igraph")

skip_if_no_python <- function() {
  if (!modules)
    testthat::skip("leidenalg not available for testing")
}

<<<<<<< HEAD

=======
set.seed(9000)
>>>>>>> 7d904bd98ec75ea8cc9bda7b50896458b3c26522
test_that("run with CPMVertexPartition multiplexed", {
  skip_if_no_python()
  partition <- leiden(multiplex_graph,
                      partition_type = "CPMVertexPartition",
                      resolution_parameter = 0.1,
                      seed = 9001)
  expect_length(partition, length(V(multiplex_graph[[1]])))
  expect_equal(sort(unique(partition)), c(1:10))
  expect_equal(partition,
               c(8, 2, 7, 5, 2, 5, 1, 2, 6, 5, 5, 2, 2, 5, 5, 5, 3, 6, 1, 2,
                 7, 10, 3, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 4, 4, 2,
<<<<<<< HEAD
                 2, 2, 6, 4, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 6, 7, 6, 4, 9, 4
=======
                 2, 2, 6, 4, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 3, 7, 6, 4, 9, 4
>>>>>>> 7d904bd98ec75ea8cc9bda7b50896458b3c26522
               ))
  multiplex_graph
})

<<<<<<< HEAD
=======
set.seed(9000)
>>>>>>> 7d904bd98ec75ea8cc9bda7b50896458b3c26522
test_that("run with ModularityVertexPartition multiplexed", {
  skip_if_no_python()
  partition <- leiden(multiplex_graph,
                      partition_type = "ModularityVertexPartition",
                      resolution_parameter = 0.1,
                      degree_as_node_size = TRUE,
                      seed = 9001)
  expect_length(partition, length(V(multiplex_graph[[1]])))
  expect_equal(sort(unique(partition)), 1:6)
<<<<<<< HEAD
  expect_equal(partition,
               c(3, 3, 1, 4, 3, 4, 4, 3, 1, 4, 4, 3, 3, 4, 4, 4, 1, 1, 2, 3,
                 1, 3, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 5, 5, 3, 3,
                 3, 1, 5, 3, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 1, 1, 1, 5, 6, 5))
})

=======
  expect_equal(table(partition),
               structure(c(`1` = 15L, `2` = 13L, `3` = 13L, `4` = 11L, `5` = 8L,
                           `6` = 1L), .Dim = 6L, .Dimnames = list(partition = c("1", "2",
                           "3", "4", "5", "6")), class = "table"))
  expect_equal(partition,
               c(4, 2, 4, 4, 2, 4, 4, 2, 1, 4, 4, 2, 2, 4, 4, 4, 1, 1, 3, 2,
                 2, 2, 1, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 5, 5, 2, 2,
                 2, 1, 5, 2, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 1, 4, 1, 5, 6, 5))
})

set.seed(9000)
>>>>>>> 7d904bd98ec75ea8cc9bda7b50896458b3c26522
test_that("run with ModularityVertexPartition multiplexed and max_comm_size", {
  skip_if_no_python()
  partition <- leiden(multiplex_graph,
                      partition_type = "ModularityVertexPartition",
                      resolution_parameter = 0.1,
                      max_comm_size = 8,
                      degree_as_node_size = TRUE,
                      seed = 9001)
  expect_length(partition, length(V(multiplex_graph[[1]])))
  expect_equal(sort(unique(partition)), 1:10)
  expect_equal(max(table(partition)), 8)
<<<<<<< HEAD
  expect_equal(partition,
               c(4, 4, 2, 3, 8, 3, 3, 4, 2, 3, 3, 8, 8, 3, 3, 3, 6, 6, 5, 8,
                 2, 8, 6, 7, 2, 5, 5, 5, 7, 5, 7, 7, 5, 5, 7, 5, 4, 1, 1, 4, 4,
                 4, 2, 1, 4, 6, 6, 6, 6, 9, 2, 6, 1, 1, 1, 2, 2, 7, 1, 10, 1))
=======
  expect_equal(table(partition),
               structure(c(`1` = 8L, `2` = 8L, `3` = 8L, `4` = 8L, `5` = 8L,
                           `6` = 8L, `7` = 6L, `8` = 5L, `9` = 1L, `10` = 1L), .Dim = 10L, .Dimnames = list(
                             partition = c("1", "2", "3", "4", "5", "6", "7", "8", "9",
                                           "10")), class = "table"))
>>>>>>> 7d904bd98ec75ea8cc9bda7b50896458b3c26522
})
