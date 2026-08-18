test_that("hc_pca stores scores and explained variance on the object", {
  result <- adata_pca()

  expect_s3_class(result, "AnnDatR")
  expect_equal(ncol(result[["obsm"]][["X_pca"]]), 40L)
  expect_equal(nrow(result[["obsm"]][["X_pca"]]), example_adata[["n_obs"]])
  expect_setequal(
    names(result[["uns"]][["pca"]]),
    c("scores", "loadings", "sdev", "r2cum", "method")
  )
  expect_equal(result[["uns"]][["pca"]][["method"]], "svd")
  expect_identical(
    result[["obsm"]][["X_pca"]],
    result[["uns"]][["pca"]][["scores"]]
  )
})

test_that("hc_pca leaves the input object untouched", {
  before <- names(example_adata[["uns"]])
  hc_pca(example_adata, components = 5, verbose = FALSE)
  expect_equal(names(example_adata[["uns"]]), before)
})

test_that("hc_pca reports missing values and imputes them by default", {
  expect_message(
    hc_pca(example_adata, components = 5),
    "values .* are missing"
  )
  result <- suppressMessages(hc_pca(example_adata, components = 5))
  expect_false(anyNA(result[["obsm"]][["X_pca"]]))
})

test_that("hc_pca honours na_action", {
  expect_error(
    suppressMessages(hc_pca(example_adata, components = 5, na_action = "error")),
    "missing values found"
  )

  omitted <- suppressMessages(
    hc_pca(example_adata, components = 5, na_action = "omit")
  )
  expect_lt(nrow(omitted[["obsm"]][["X_pca"]]), example_adata[["n_obs"]])
  expect_false(anyNA(omitted[["obsm"]][["X_pca"]]))
})

test_that("hc_pca applies the requested transformation", {
  expect_error(
    hc_pca(example_adata, components = 5, transform = "cube"),
    "'arg' should be one of"
  )

  for (transform in c("log1p", "sqrt", "none")) {
    result <- suppressMessages(
      hc_pca(example_adata, components = 5, transform = transform, verbose = FALSE)
    )
    expect_equal(ncol(result[["obsm"]][["X_pca"]]), 5L, info = transform)
  }
})

test_that("hc_pca caps components at the number of observations", {
  expect_warning(
    suppressMessages(hc_pca(toy_adata(), components = 500)),
    "exceeds number of observations"
  )
})

test_that("hc_pca errors on an unknown layer", {
  expect_error(
    hc_pca(example_adata, components = 5, layer = "missing"),
    "not found in AnnDatR\\$layers"
  )
})

test_that("hc_kaisers_rule needs a PCA and returns a usable component count", {
  expect_error(hc_kaisers_rule(example_adata), "Call `hc_pca\\(\\)`")

  n_components <- suppressMessages(hc_kaisers_rule(adata_pca()))
  expect_true(is.numeric(n_components))
  expect_gte(n_components, 1L)
  expect_lte(n_components, ncol(adata_pca()[["obsm"]][["X_pca"]]))
})

test_that("hc_kaisers_rule falls back to the 80% variance rule", {
  fit <- adata_pca()
  # Force the Kaiser component to sit below the 80% threshold.
  fit[["uns"]][["pca"]][["sdev"]] <- c(2, 0.5, 0.4, 0.3)
  fit[["uns"]][["pca"]][["r2cum"]] <- c(0.3, 0.5, 0.85, 0.95)

  expect_message(
    result <- hc_kaisers_rule(fit, with_alternative = TRUE),
    "under 80%"
  )
  expect_equal(result, 3L)

  expect_message(
    fallback <- hc_kaisers_rule(fit, with_alternative = FALSE),
    "less than 80%"
  )
  expect_equal(fallback, 2L)
})

test_that("hc_kaisers_rule errors when no eigenvalue drops below 1", {
  fit <- adata_pca()
  fit[["uns"]][["pca"]][["sdev"]] <- c(3, 2, 1.5)
  fit[["uns"]][["pca"]][["r2cum"]] <- c(0.5, 0.8, 0.99)
  expect_error(hc_kaisers_rule(fit), "No eigenvalue is lower than 1")
})
