test_that("rescale_umap centres both axes and preserves the aspect ratio", {
  embedding <- cbind(c(0, 10, 5), c(0, 2, 1))
  rescaled <- rescale_umap(embedding)

  expect_equal(dim(rescaled), dim(embedding))
  expect_equal(range(rescaled[, 1]), c(-1, 1))
  # Each axis is centred on its own midpoint.
  expect_equal(mean(range(rescaled[, 1])), 0)
  expect_equal(mean(range(rescaled[, 2])), 0)
  # The narrower axis stays narrower: the scaling is shared, not per-axis.
  expect_lt(diff(range(rescaled[, 2])), diff(range(rescaled[, 1])))
})

test_that("hc_umap requires an SNN graph", {
  expect_error(hc_umap(adata_pca()), "Call `hc_snn\\(\\)`")
})

test_that("hc_umap stores raw and rescaled embeddings and updates obs", {
  result <- adata_umap()

  raw <- result[["obsm"]][["X_umap_raw"]]
  rescaled <- result[["obsm"]][["X_umap"]]

  expect_equal(dim(raw), c(example_adata[["n_obs"]], 2L))
  expect_equal(dim(rescaled), c(example_adata[["n_obs"]], 2L))
  expect_equal(colnames(rescaled), c("UMAP_1", "UMAP_2"))
  expect_equal(rownames(rescaled), result[["obs_names"]])
  expect_false(anyNA(rescaled))
  expect_lte(max(abs(rescaled)), 1)

  expect_true(all(c("UMAP1", "UMAP2") %in% names(result[["obs"]])))
  expect_false(anyNA(result[["obs"]][["UMAP1"]]))
  expect_equal(nrow(result[["obs"]]), example_adata[["n_obs"]])
})

test_that("hc_umap is reproducible for a fixed seed", {
  first <- suppressMessages(hc_umap(adata_clustered(), seed = 7, verbose = FALSE))
  again <- suppressMessages(hc_umap(adata_clustered(), seed = 7, verbose = FALSE))
  expect_equal(first[["obsm"]][["X_umap_raw"]], again[["obsm"]][["X_umap_raw"]])
})

test_that("hc_umap does not duplicate UMAP columns when re-run", {
  twice <- suppressMessages(hc_umap(adata_umap(), verbose = FALSE))
  expect_equal(sum(names(twice[["obs"]]) == "UMAP1"), 1L)
  expect_false(any(grepl("UMAP1\\.(x|y)$", names(twice[["obs"]]))))
})

test_that("get_density returns a tidy grid of the requested size", {
  set.seed(15)
  density <- get_density(
    stats::rnorm(50),
    stats::rnorm(50),
    h = 0.5,
    n = 20,
    lims = c(-3, 3, -3, 3)
  )

  expect_equal(nrow(density), 400L)
  expect_setequal(names(density), c("x", "y", "z", "x_coord", "y_coord"))
  expect_equal(sort(unique(density[["x"]])), 1:20)
  expect_true(all(density[["z"]] >= 0))

  # Grid indices and coordinates must stay in step.
  grid <- seq(-3, 3, length.out = 20)
  expect_equal(density[["x_coord"]], grid[density[["x"]]])
  expect_equal(density[["y_coord"]], grid[density[["y"]]])
})

test_that("hc_cluster_hulls requires UMAP coordinates and clusters", {
  expect_error(hc_cluster_hulls(adata_clustered()), "Call `hc_umap\\(\\)`")

  without_clusters <- suppressMessages(hc_umap(
    suppressMessages(hc_snn(
      hc_distance(adata_pca(), components = 20),
      neighbors = 15,
      verbose = FALSE
    )),
    verbose = FALSE
  ))
  expect_error(
    hc_cluster_hulls(without_clusters),
    "Call `hc_cluster_consensus\\(\\)`"
  )
})

test_that("hc_cluster_hulls produces closed polygons per cluster", {
  result <- hc_cluster_hulls(adata_umap(), n = 120)
  hulls <- result[["uns"]][["UMAP_hulls"]]

  expect_setequal(
    names(hulls),
    c("hulls", "density", "landmass_pixels", "center_density")
  )

  polygons <- hulls[["hulls"]]
  expect_true(all(c("X", "Y", "cluster", "polygon_id") %in% names(polygons)))
  expect_gt(dplyr::n_distinct(polygons[["polygon_id"]]), 1L)
  expect_false(anyNA(polygons[["X"]]))

  # Degenerate one- and two-pixel landmasses are dropped, so every polygon is
  # a closed ring of at least three distinct vertices.
  vertex_counts <- table(polygons[["polygon_id"]])
  expect_true(all(vertex_counts >= 4L))

  # Hull coordinates stay inside the plotting range.
  limit <- 1.05 * max(abs(adata_umap()[["obsm"]][["X_umap"]]))
  expect_lte(max(abs(polygons[["X"]])), limit)
  expect_lte(max(abs(polygons[["Y"]])), limit)

  # Every cluster with a hull is a real cluster.
  expect_true(all(
    polygons[["cluster"]] %in% adata_umap()[["obs"]][["cluster"]]
  ))

  # UMAP columns are consumed by this step.
  expect_false(any(c("UMAP1", "UMAP2") %in% names(result[["obs"]])))
})

test_that("hc_cluster_hulls reads the observation names column, not a fixed name", {
  # Regression: the element id used to be hard-coded as `ensembl_id`.
  renamed <- adata_umap()$clone(deep = TRUE)
  names(renamed[["obs"]])[names(renamed[["obs"]]) == "ensembl_id"] <- "gene_id"
  renamed[["obs_names_col"]] <- "gene_id"

  expect_no_error(hc_cluster_hulls(renamed, n = 80))
})
