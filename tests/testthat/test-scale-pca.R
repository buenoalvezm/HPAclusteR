test_that("prep_matrix centres and scales as documented", {
  x <- matrix(c(1, 2, 3, 10, 20, 30), ncol = 2)

  centred <- prep_matrix(x, scale = "none")
  expect_equal(colMeans(centred), c(0, 0))

  unit_variance <- prep_matrix(x, scale = "uv")
  expect_equal(apply(unit_variance, 2, stats::sd), c(1, 1))

  pareto <- prep_matrix(x, scale = "pareto")
  expect_equal(
    pareto,
    sweep(centred, 2, sqrt(apply(centred, 2, stats::sd)), "/")
  )

  vector_scaled <- prep_matrix(x, scale = "vector")
  expect_equal(
    vector_scaled,
    sweep(centred, 2, sqrt(colSums(centred^2)), "/")
  )
})

test_that("prep_matrix leaves constant columns alone and warns", {
  x <- cbind(c(1, 1, 1), c(1, 2, 3))
  expect_warning(result <- prep_matrix(x, scale = "uv"), "below eps")
  expect_equal(result[, 1], c(0, 0, 0))
  expect_false(any(is.nan(result)))
})

test_that("prep_matrix ignores missing values when computing statistics", {
  x <- cbind(c(1, 2, NA, 4), c(1, 2, 3, 4))
  result <- prep_matrix(x, scale = "uv")
  expect_true(is.na(result[3, 1]))
  expect_equal(mean(result[-3, 1]), 0)
})

test_that("hc_scale returns genes in rows for every method", {
  set.seed(2)
  # Samples in rows, genes in columns, as hc_pca passes it.
  data <- matrix(stats::runif(40, 1, 10), nrow = 5, ncol = 8)

  methods <- c(
    "sample", "gene", "pareto", "vector",
    "min-max", "max_sequencial", "max"
  )
  for (method in methods) {
    scaled <- hc_scale(data, method = method)
    expect_equal(dim(scaled), c(8L, 5L), info = method)
  }
})

test_that("hc_scale computes the documented transformations", {
  data <- matrix(c(1, 2, 4, 8), nrow = 2)

  expect_equal(hc_scale(data, "max"), t(data) / apply(t(data), 1, max))
  expect_equal(
    hc_scale(data, "max_sequencial"),
    t(data) / apply(t(data), 1, max)
  )

  min_max <- hc_scale(data, "min-max")
  expect_equal(apply(min_max, 1, min), rep(0, nrow(min_max)))
  expect_equal(apply(min_max, 1, max), rep(1, nrow(min_max)))
})

test_that("hc_scale rejects unknown methods", {
  expect_error(hc_scale(matrix(1:4, 2), "nope"), "Scaling method not defined")
})

test_that("pca_svd agrees with prcomp and reports cumulative variance", {
  set.seed(8)
  x <- matrix(stats::rnorm(200), nrow = 40)

  fit <- pca_svd(x, n_pcs = 3)
  reference <- stats::prcomp(x, center = TRUE, scale. = FALSE)

  expect_equal(fit[["scores"]], reference[["x"]][, 1:3])
  expect_equal(fit[["sdev"]], apply(fit[["scores"]], 2, stats::sd))
  expect_equal(
    fit[["r2cum"]],
    (cumsum(reference[["sdev"]]^2) / sum(reference[["sdev"]]^2))[1:3]
  )
  expect_true(all(diff(fit[["r2cum"]]) > 0))
  expect_lte(max(fit[["r2cum"]]), 1)
  expect_equal(fit[["method"]], "svd")
})

test_that("pca_svd caps components at the rank of the input", {
  x <- matrix(stats::rnorm(30), nrow = 10)
  expect_equal(ncol(pca_svd(x, n_pcs = 50)[["scores"]]), 3L)
})

test_that("impute_column_means fills gaps with the column mean", {
  x <- cbind(c(1, 2, NA, 5), c(NA, NA, NA, NA))
  result <- impute_column_means(x)

  expect_false(anyNA(result))
  expect_equal(result[3, 1], mean(c(1, 2, 5)))
  expect_equal(result[, 2], rep(0, 4))
})

test_that("impute_column_means is a no-op on complete matrices", {
  x <- matrix(1:6, nrow = 3)
  expect_identical(impute_column_means(x), x)
})
