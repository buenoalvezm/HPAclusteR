test_that("cluster_genes returns one row per gene, numbered from zero", {
  graph <- snn_to_igraph(adata_snn()[["uns"]][["neighbors"]][["snn"]])
  genes <- adata_snn()[["obs_names"]]

  partition <- cluster_genes(genes, graph, method = "louvain", resolution = 8, seed = 1)

  expect_equal(nrow(partition), length(genes))
  expect_setequal(names(partition), c("gene", "cluster"))
  expect_equal(partition[["gene"]], genes)
  expect_equal(min(partition[["cluster"]]), 0L)
  expect_type(partition[["cluster"]], "integer")
})

test_that("cluster_genes is reproducible for a given seed and varies across seeds", {
  graph <- snn_to_igraph(adata_snn()[["uns"]][["neighbors"]][["snn"]])
  genes <- adata_snn()[["obs_names"]]

  first <- cluster_genes(genes, graph, resolution = 8, seed = 3)
  again <- cluster_genes(genes, graph, resolution = 8, seed = 3)
  other <- cluster_genes(genes, graph, resolution = 8, seed = 4)

  expect_identical(first, again)
  # Different seeds must be able to disagree, otherwise the consensus step
  # would be aggregating 100 copies of the same partition.
  expect_false(identical(first[["cluster"]], other[["cluster"]]))
})

test_that("cluster_genes supports both algorithms", {
  graph <- snn_to_igraph(adata_snn()[["uns"]][["neighbors"]][["snn"]])
  genes <- adata_snn()[["obs_names"]]

  leiden <- cluster_genes(genes, graph, method = "leiden", resolution = 8, seed = 1)
  expect_equal(nrow(leiden), length(genes))
  expect_gt(length(unique(leiden[["cluster"]])), 1L)

  expect_error(
    cluster_genes(genes, graph, method = "walktrap"),
    "'arg' should be one of"
  )
})

test_that("higher resolution yields more communities", {
  graph <- snn_to_igraph(adata_snn()[["uns"]][["neighbors"]][["snn"]])
  genes <- adata_snn()[["obs_names"]]

  low <- length(unique(cluster_genes(genes, graph, resolution = 4, seed = 1)[["cluster"]]))
  high <- length(unique(cluster_genes(genes, graph, resolution = 12, seed = 1)[["cluster"]]))

  expect_gt(high, low)
})

test_that("hc_cluster_consensus requires an SNN graph", {
  expect_error(
    hc_cluster_consensus(adata_pca()),
    "Call `hc_snn\\(\\)`"
  )
})

test_that("hc_cluster_consensus populates uns and obs", {
  result <- adata_clustered()

  expect_true(all(
    c("cluster_data", "consensus_clustering", "membership_matrix", "mapping_table") %in%
      names(result[["uns"]])
  ))

  consensus <- result[["uns"]][["consensus_clustering"]]
  expect_setequal(names(consensus), c("gene", "cluster"))
  expect_equal(nrow(consensus), result[["n_obs"]])
  expect_false(anyNA(consensus[["cluster"]]))

  expect_true(all(c("cluster", "cluster_colors") %in% names(result[["obs"]])))
  expect_false(anyNA(result[["obs"]][["cluster_colors"]]))
  expect_equal(nrow(result[["obs"]]), example_adata[["n_obs"]])
})

test_that("every cluster is assigned a colour, recycling the palette", {
  obs <- adata_clustered()[["obs"]]
  mapping <- unique(obs[, c("cluster", "cluster_colors")])
  expect_equal(nrow(mapping), dplyr::n_distinct(obs[["cluster"]]))
  expect_true(all(grepl("^#", mapping[["cluster_colors"]])))
})

test_that("cluster_data holds one column per seed", {
  cluster_data <- adata_clustered()[["uns"]][["cluster_data"]]
  seed_columns <- grep("^seed_", names(cluster_data), value = TRUE)
  expect_equal(length(seed_columns), 20L)
  expect_equal(nrow(cluster_data), example_adata[["n_obs"]])
})

test_that("hc_cluster_consensus rejects resolutions outside the k window", {
  expect_error(
    suppressMessages(hc_cluster_consensus(
      adata_snn(),
      resolution = 1,
      n_seeds = 5,
      verbose = FALSE
    )),
    "outside the accepted range"
  )

  # The window is configurable rather than hard-coded.
  expect_no_error(suppressMessages(hc_cluster_consensus(
    adata_snn(),
    resolution = 1,
    n_seeds = 5,
    min_k = 1,
    max_k = 500,
    verbose = FALSE
  )))
})

test_that("hc_cluster_consensus can be re-run on an already clustered object", {
  # Regression test: the join used to collide with the existing `cluster`
  # column and fail with "Join columns in `x` must be present in the data".
  reclustered <- suppressMessages(hc_cluster_consensus(
    adata_clustered(),
    resolution = 10,
    n_seeds = 10,
    verbose = FALSE
  ))

  expect_true("cluster" %in% names(reclustered[["obs"]]))
  expect_false(any(grepl("\\.(x|y)$", names(reclustered[["obs"]]))))
  expect_equal(sum(names(reclustered[["obs"]]) == "cluster"), 1L)
  expect_false(anyNA(reclustered[["obs"]][["cluster"]]))
})

test_that("hc_cluster_consensus is reproducible", {
  first <- suppressMessages(hc_cluster_consensus(
    adata_snn(), resolution = 8, n_seeds = 10, verbose = FALSE
  ))
  again <- suppressMessages(hc_cluster_consensus(
    adata_snn(), resolution = 8, n_seeds = 10, verbose = FALSE
  ))
  expect_identical(
    first[["uns"]][["consensus_clustering"]],
    again[["uns"]][["consensus_clustering"]]
  )
})

test_that("calculate_ari matches known reference values", {
  expect_equal(calculate_ari(c(1, 1, 2, 2), c(1, 1, 2, 2)), 1)
  expect_equal(calculate_ari(c(1, 1, 2, 2), c(2, 2, 1, 1)), 1)

  # Independent labellings sit near zero.
  set.seed(21)
  a <- sample(1:4, 200, replace = TRUE)
  b <- sample(1:4, 200, replace = TRUE)
  expect_lt(abs(calculate_ari(a, b)), 0.1)

  # Worked example: ARI of a partial agreement.
  expect_equal(
    calculate_ari(c(1, 1, 1, 2, 2, 2), c(1, 1, 2, 2, 3, 3)),
    0.2424242,
    tolerance = 1e-6
  )
})

test_that("hc_cluster_stability needs cluster_data and scores every seed pair", {
  expect_error(hc_cluster_stability(adata_snn()), "Call `hc_cluster_consensus\\(\\)`")

  stability <- hc_cluster_stability(adata_clustered())

  expect_setequal(names(stability), c("cluster_stability_hist", "ari_scores"))
  expect_equal(nrow(stability[["ari_scores"]]), choose(20, 2))
  expect_setequal(names(stability[["ari_scores"]]), c("seed1", "seed2", "ARI"))
  expect_true(all(stability[["ari_scores"]][["seed1"]] < stability[["ari_scores"]][["seed2"]]))
  expect_true(all(stability[["ari_scores"]][["ARI"]] <= 1))
  expect_s3_class(stability[["cluster_stability_hist"]], "ggplot")
})
