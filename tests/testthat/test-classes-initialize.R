test_that("check_annotation_alignment always returns a single logical", {
  expect_true(check_annotation_alignment(example_adata))

  # Regression: this used to fall off the end and return NULL, so the caller's
  # `if (!check_annotation_alignment(...))` failed with "argument is of length
  # zero" instead of reporting the misalignment.
  broken <- example_adata$clone(deep = TRUE)
  broken[["obs"]] <- broken[["obs"]][rev(seq_len(nrow(broken[["obs"]]))), ]

  result <- check_annotation_alignment(broken)
  expect_type(result, "logical")
  expect_length(result, 1L)
  expect_false(result)
})

test_that("AnnDatR validation rejects misaligned objects with a clear message", {
  expect_error(
    AnnDatR$new(
      X = example_adata[["X"]],
      obs = example_adata[["obs"]][rev(seq_len(example_adata[["n_obs"]])), ],
      var = example_adata[["var"]],
      var_names = "sample_id"
    ),
    "alignment failed"
  )
})

test_that("AnnDatR validation catches dimension mismatches", {
  expect_error(
    AnnDatR$new(
      X = example_adata[["X"]],
      obs = example_adata[["obs"]][1:10, ],
      var = example_adata[["var"]],
      var_names = "sample_id"
    ),
    "columns"
  )
})

test_that("AnnDatR requires either a prefix or all three tables", {
  expect_error(AnnDatR$new(), "Provide either")
})

test_that("AnnDatR active bindings report the right shape", {
  # Built fresh, not read from data/: R6 active bindings are lost on
  # serialisation, so the shipped dataset carries plain frozen fields.
  adata <- toy_adata(n_genes = 12L, n_samples = 6L)

  expect_equal(adata[["n_obs"]], nrow(adata[["obs"]]))
  expect_equal(adata[["n_vars"]], nrow(adata[["var"]]))
  expect_equal(adata[["shape"]], c(adata[["n_obs"]], adata[["n_vars"]]))
  expect_equal(
    adata[["obs_names"]],
    adata[["obs"]][[adata[["obs_names_col"]]]]
  )
  expect_error(adata$obs_names <- "x", "not yet supported")
  expect_error(adata$var_names <- "x", "not yet supported")
})

test_that("active bindings track changes to obs and var", {
  adata <- toy_adata(n_genes = 12L, n_samples = 6L)
  adata[["obs"]] <- adata[["obs"]][1:5, ]

  expect_equal(adata[["n_obs"]], 5L)
  expect_length(adata[["obs_names"]], 5L)
})

test_that("clustering reads gene names from obs, not the active binding", {
  # Regression: reading `obs_names` returned stale values for objects restored
  # from disk, because R6 active bindings do not survive serialisation.
  restored <- example_adata$clone(deep = TRUE)
  expect_equal(
    restored[["obs"]][[restored[["obs_names_col"]]]],
    adata_clustered()[["uns"]][["consensus_clustering"]][["gene"]]
  )
})

test_that("hc_initialize transposes genes into observations", {
  adata <- toy_adata(n_genes = 12L, n_samples = 6L)

  expect_s3_class(adata, "AnnDatR")
  expect_equal(adata[["n_obs"]], 12L)
  expect_equal(adata[["n_vars"]], 6L)
  expect_equal(adata[["obs_names_col"]], "ensembl_id")
  expect_equal(adata[["var_names_col"]], "sample_id")

  # X is samples x (genes + 1) after transposition.
  expect_equal(dim(adata[["X"]]), c(6L, 13L))
  expect_equal(colnames(adata[["X"]])[-1], adata[["obs_names"]])
  expect_equal(adata[["X"]][[1]], adata[["var_names"]])
})

test_that("hc_initialize preserves the expression values", {
  set.seed(30)
  X <- tibble::tibble(
    ensembl_id = c("g1", "g2", "g3"),
    s1 = c(1, 2, 3),
    s2 = c(4, 5, 6)
  )
  obs <- tibble::tibble(ensembl_id = c("g1", "g2", "g3"))
  var <- tibble::tibble(sample_id = c("s1", "s2"))

  adata <- hc_initialize(X, obs = obs, var = var)

  expect_equal(adata[["X"]][["g1"]], c(1, 4))
  expect_equal(adata[["X"]][["g3"]], c(3, 6))
})

test_that("hc_initialize reorders metadata to match the matrix", {
  X <- tibble::tibble(
    ensembl_id = c("g1", "g2", "g3"),
    s1 = c(1, 2, 3),
    s2 = c(4, 5, 6)
  )
  # Metadata deliberately supplied in a different order.
  obs <- tibble::tibble(ensembl_id = c("g3", "g1", "g2"), tag = c("c", "a", "b"))
  var <- tibble::tibble(sample_id = c("s2", "s1"), tag = c("B", "A"))

  adata <- hc_initialize(X, obs = obs, var = var)

  expect_equal(adata[["obs"]][["tag"]], c("a", "b", "c"))
  expect_equal(adata[["var"]][["tag"]], c("A", "B"))
})

test_that("the toy fixture runs through the clustering pipeline", {
  adata <- toy_adata(n_genes = 60L, n_samples = 9L)

  fit <- suppressWarnings(suppressMessages(
    hc_pca(adata, components = 5, verbose = FALSE)
  ))
  fit <- hc_distance(fit, components = 3)
  fit <- suppressWarnings(hc_snn(fit, neighbors = 5, verbose = FALSE))

  expect_equal(dim(fit[["uns"]][["neighbors"]][["snn"]]), c(60L, 60L))
})
