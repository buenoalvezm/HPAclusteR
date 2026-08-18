clustered_b <- function() {
  cached("clustered_b", function() {
    suppressMessages(hc_cluster_consensus(
      adata_snn(),
      resolution = 10,
      n_seeds = 20,
      verbose = FALSE
    ))
  })
}

test_that("hc_cluster_compare requires consensus clustering on both objects", {
  expect_error(
    hc_cluster_compare(adata_snn(), adata_clustered()),
    "first AnnDatR object"
  )
  expect_error(
    hc_cluster_compare(adata_clustered(), adata_snn()),
    "second AnnDatR object"
  )
})

test_that("hc_cluster_compare reports overlaps with valid statistics", {
  result <- suppressMessages(
    hc_cluster_compare(adata_clustered(), clustered_b())
  )

  expect_setequal(names(result), c("matches", "heatmap", "network"))

  matches <- result[["matches"]]
  expect_true(all(
    c("cluster_A", "cluster_B", "n_overlap", "percentage_overlap", "p_val", "adj_p_val") %in%
      names(matches)
  ))
  expect_true(all(matches[["adj_p_val"]] < 0.05))
  expect_true(all(matches[["p_val"]] >= 0 & matches[["p_val"]] <= 1))
  expect_true(all(matches[["n_overlap"]] > 0))
  expect_true(all(
    matches[["percentage_overlap"]] >= 0 & matches[["percentage_overlap"]] <= 100
  ))

  expect_s3_class(result[["heatmap"]], "ggplot")
  expect_s3_class(result[["network"]], "ggplot")
})

test_that("hc_cluster_compare optimized mode adds a match score", {
  result <- suppressMessages(
    hc_cluster_compare(adata_clustered(), clustered_b(), graph_type = "optimized")
  )

  expect_true("match_score" %in% names(result))
  expect_gte(result[["match_score"]], 0)
  expect_lte(result[["match_score"]], 100)
  expect_lte(nrow(result[["matches"]]), nrow(
    suppressMessages(hc_cluster_compare(adata_clustered(), clustered_b()))[["matches"]]
  ))
})

test_that("cluster_matching never produces many-to-many matches", {
  comparison <- tibble::tibble(
    cluster_A = c("1", "1", "2", "2", "3"),
    cluster_B = c("a", "b", "a", "b", "c"),
    n_overlap = c(100, 90, 80, 70, 60),
    adj_p_val = rep(0.001, 5)
  )

  result <- cluster_matching(comparison, total_genes = 400)
  matched <- result[["matched_clusters"]]

  # A pair is rejected only when both of its clusters are already matched.
  expect_false(any(duplicated(paste(matched[["cluster_A"]], matched[["cluster_B"]]))))
  expect_equal(result[["match_score"]], round(sum(matched[["n_overlap"]]) / 400 * 100, 2))
  expect_lte(nrow(matched), nrow(comparison))
})

test_that("cluster_matching drops non-significant pairs", {
  comparison <- tibble::tibble(
    cluster_A = c("1", "2"),
    cluster_B = c("a", "b"),
    n_overlap = c(100, 90),
    adj_p_val = c(0.001, 0.9)
  )

  matched <- cluster_matching(comparison, total_genes = 200)[["matched_clusters"]]
  expect_equal(nrow(matched), 1L)
  expect_equal(matched[["cluster_A"]], "1")
})

test_that("visualize_comparison_net rejects an empty match table", {
  empty <- tibble::tibble(
    cluster_A = character(), cluster_B = character(),
    n_genes_A = numeric(), n_genes_B = numeric(),
    percentage_overlap = numeric()
  )
  expect_error(visualize_comparison_net(empty), "No significant cluster matches")
})

test_that("hc_plot_umap draws points, hulls and both", {
  with_hulls <- hc_cluster_hulls(adata_umap(), n = 120)

  expect_error(hc_plot_umap(adata_clustered()), "run `hc_umap\\(\\)`")
  expect_error(hc_plot_umap(adata_umap(), plot = "hulls"), "run `hc_cluster_hulls\\(\\)`")

  for (style in c("points", "hulls", "both")) {
    plot <- hc_plot_umap(with_hulls, plot = style)
    expect_s3_class(plot, "ggplot")
    expect_no_error(ggplot2::ggplot_build(plot))
  }
})

test_that("hc_plot_expression returns one heatmap pair per cluster", {
  expect_error(
    hc_plot_expression(adata_snn()),
    "Call `hc_cluster_consensus\\(\\)`"
  )

  plots <- hc_plot_expression(adata_clustered(), show_sample_labels = FALSE)

  expect_setequal(names(plots), c("zscore", "scaled"))
  clusters <- as.character(sort(unique(adata_clustered()[["obs"]][["cluster"]])))
  expect_setequal(names(plots[["zscore"]]), clusters)
  expect_setequal(names(plots[["scaled"]]), clusters)
  expect_no_error(print(plots[["zscore"]][[1]]))
})

test_that("theme_hc returns a ggplot2 theme", {
  expect_s3_class(theme_hc(), "theme")
  plot <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) +
    ggplot2::geom_point() +
    theme_hc()
  expect_no_error(ggplot2::ggplot_build(plot))
})

test_that("visualize_ari returns the histogram and the scores", {
  scores <- tibble::tibble(seed1 = 1L, seed2 = 2L, ARI = 0.5)
  result <- visualize_ari(scores)

  expect_setequal(names(result), c("cluster_stability_hist", "ari_scores"))
  expect_identical(result[["ari_scores"]], scores)
  expect_s3_class(result[["cluster_stability_hist"]], "ggplot")
})
