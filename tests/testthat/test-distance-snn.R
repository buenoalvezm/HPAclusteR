test_that("hc_distance requires a PCA first", {
  expect_error(hc_distance(example_adata), "Call `hc_pca\\(\\)`")
})

test_that("hc_distance computes distances over the requested components", {
  result <- hc_distance(adata_pca(), components = 20)
  distances <- result[["uns"]][["distance"]]

  expect_s3_class(distances, "dist")
  expect_equal(attr(distances, "Size"), example_adata[["n_obs"]])
  expect_equal(
    as.vector(distances),
    as.vector(hc_dist(adata_pca()[["obsm"]][["X_pca"]][, 1:20], method = "spearman"))
  )
})

test_that("hc_distance caps components and warns", {
  expect_warning(
    result <- hc_distance(adata_pca(), components = 1000),
    "exceeds number of PCA components"
  )
  expect_equal(
    as.vector(result[["uns"]][["distance"]]),
    as.vector(hc_dist(adata_pca()[["obsm"]][["X_pca"]], method = "spearman"))
  )
})

test_that("hc_distance picks the 80% variance component when unspecified", {
  expect_warning(
    result <- hc_distance(adata_pca()),
    "first component with over 80%"
  )
  expect_s3_class(result[["uns"]][["distance"]], "dist")
})

test_that("hc_distance supports every documented metric", {
  for (method in c("euclidean", "manhattan", "maximum", "pearson", "spearman")) {
    result <- hc_distance(adata_pca(), components = 5, method = method)
    expect_s3_class(result[["uns"]][["distance"]], "dist")
  }
})

test_that("build_snn produces a symmetric, pruned Jaccard graph", {
  set.seed(12)
  points <- matrix(stats::rnorm(60), nrow = 20)
  rownames(points) <- paste0("g", seq_len(20))
  distances <- stats::dist(points)

  k <- 5L
  graphs <- build_snn(distances, k = k, prune = 0)

  expect_setequal(names(graphs), c("nn", "snn"))
  expect_equal(dim(graphs[["snn"]]), c(20L, 20L))
  expect_equal(rownames(graphs[["snn"]]), rownames(points))
  expect_true(Matrix::isSymmetric(graphs[["snn"]]))

  # Each row of the nearest neighbour graph holds exactly k ones.
  expect_equal(unname(Matrix::rowSums(graphs[["nn"]])), rep(k, 20))

  # Self-similarity is 1: a set shares all k neighbours with itself.
  expect_equal(unname(Matrix::diag(graphs[["snn"]])), rep(1, 20))

  # Weights are the Jaccard index, so they lie in (0, 1].
  weights <- graphs[["snn"]]@x
  expect_true(all(weights > 0 & weights <= 1))
})

test_that("build_snn prunes weak edges", {
  set.seed(13)
  points <- matrix(stats::rnorm(80), nrow = 40)
  rownames(points) <- paste0("g", seq_len(40))
  distances <- stats::dist(points)

  loose <- build_snn(distances, k = 8L, prune = 0)
  strict <- build_snn(distances, k = 8L, prune = 0.5)

  expect_lt(length(strict[["snn"]]@x), length(loose[["snn"]]@x))
  expect_true(all(strict[["snn"]]@x >= 0.5))
})

test_that("build_snn reproduces the Jaccard formula exactly", {
  set.seed(14)
  points <- matrix(stats::rnorm(30), nrow = 10)
  rownames(points) <- paste0("g", seq_len(10))
  distance_matrix <- as.matrix(stats::dist(points))

  k <- 4L
  graphs <- build_snn(stats::as.dist(distance_matrix), k = k, prune = 0)

  neighbours <- lapply(seq_len(10), function(i) order(distance_matrix[i, ])[1:k])
  overlap <- length(intersect(neighbours[[1]], neighbours[[2]]))

  expect_equal(
    as.matrix(graphs[["snn"]])[1, 2],
    overlap / (2 * k - overlap)
  )
})

test_that("build_snn caps k at the number of observations", {
  points <- matrix(stats::rnorm(20), nrow = 10)
  rownames(points) <- paste0("g", seq_len(10))
  expect_warning(
    build_snn(stats::dist(points), k = 50L, prune = 0),
    "larger than the number of observations"
  )
})

test_that("hc_snn requires a distance matrix and stores both graphs", {
  expect_error(hc_snn(adata_pca()), "Call `hc_distance\\(\\)`")

  result <- adata_snn()
  expect_setequal(names(result[["uns"]][["neighbors"]]), c("nn", "snn"))
  expect_equal(
    rownames(result[["uns"]][["neighbors"]][["snn"]]),
    result[["obs_names"]]
  )
})
