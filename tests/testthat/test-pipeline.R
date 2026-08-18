# End-to-end runs. These are the tests that would have caught the regressions
# the unit tests above pin down individually.

test_that("the documented step-by-step pipeline runs start to finish", {
  adata <- suppressMessages(hc_pca(example_adata, components = 40, verbose = FALSE))
  adata <- hc_distance(adata, components = hc_kaisers_rule(adata, verbose = FALSE))
  adata <- hc_snn(adata, neighbors = 15, verbose = FALSE)
  adata <- suppressMessages(
    hc_cluster_consensus(adata, resolution = 8, n_seeds = 20, verbose = FALSE)
  )
  adata <- suppressMessages(hc_umap(adata, verbose = FALSE))
  adata <- hc_cluster_hulls(adata, n = 150)

  expect_s3_class(adata, "AnnDatR")
  expect_true(all(
    c("pca", "distance", "neighbors", "consensus_clustering", "UMAP_hulls") %in%
      names(adata[["uns"]])
  ))
  expect_true(all(c("X_pca", "X_umap", "X_umap_raw") %in% names(adata[["obsm"]])))

  # Every gene keeps a cluster all the way through.
  expect_equal(nrow(adata[["obs"]]), example_adata[["n_obs"]])
  expect_false(anyNA(adata[["obs"]][["cluster"]]))

  expect_s3_class(hc_plot_umap(adata, plot = "both"), "ggplot")
})

test_that("hc_auto_cluster reproduces the step-by-step pipeline", {
  set.seed(1)
  auto <- suppressMessages(hc_auto_cluster(
    example_adata,
    pca_components = 40,
    distance_components = 20,
    snn_neighbors = 15,
    cluster_resolution = 8,
    cluster_seeds = 20,
    verbose = FALSE
  ))

  manual <- suppressMessages(hc_pca(example_adata, components = 40, verbose = FALSE))
  manual <- hc_distance(manual, components = 20)
  manual <- hc_snn(manual, neighbors = 15, verbose = FALSE)
  manual <- suppressMessages(
    hc_cluster_consensus(manual, resolution = 8, n_seeds = 20, verbose = FALSE)
  )

  expect_equal(
    auto[["uns"]][["consensus_clustering"]],
    manual[["uns"]][["consensus_clustering"]]
  )
  expect_true("UMAP_hulls" %in% names(auto[["uns"]]))
})

test_that("the pipeline never mutates its input object", {
  snapshot <- list(
    uns = names(example_adata[["uns"]]),
    obsm = names(example_adata[["obsm"]]),
    obs = names(example_adata[["obs"]]),
    n_obs = example_adata[["n_obs"]]
  )

  invisible(suppressMessages(hc_auto_cluster(
    example_adata,
    pca_components = 20,
    distance_components = 10,
    snn_neighbors = 15,
    cluster_resolution = 8,
    cluster_seeds = 20,
    verbose = FALSE
  )))

  expect_equal(names(example_adata[["uns"]]), snapshot[["uns"]])
  expect_equal(names(example_adata[["obsm"]]), snapshot[["obsm"]])
  expect_equal(names(example_adata[["obs"]]), snapshot[["obs"]])
  expect_equal(example_adata[["n_obs"]], snapshot[["n_obs"]])
})

test_that("clustering, classification and comparison agree on the gene set", {
  adata <- adata_clustered()

  classification <- hc_classify(adata, "tissue_name")
  classified_genes <- unlist(lapply(classification[["classification"]], `[[`, "ENSG"))

  expect_setequal(classified_genes, adata[["obs_names"]])
})

test_that("installing the package pulls in neither Python nor Seurat", {
  # The whole point of the rewrite. Checked against the declared dependency
  # closure rather than loadedNamespaces(), because the optional Bioconductor
  # annotation stack drags in unrelated packages during a test run.
  banned <- c(
    "Seurat", "reticulate", "pcaMethods", "factoextra", "concaveman",
    "fpc", "MASS", "matrixStats", "GGally", "network", "sna"
  )

  declared <- utils::packageDescription("HPAclusteR", fields = "Imports")
  declared <- trimws(sub("\\s*\\(.*\\)$", "", unlist(strsplit(declared, ","))))
  declared <- declared[nzchar(declared)]

  expect_false(any(banned %in% declared))

  # And nothing they depend on brings those back in either. Walk the installed
  # packages rather than querying a repository, so the test stays offline.
  hard_deps <- function(package) {
    description <- tryCatch(
      utils::packageDescription(package, fields = c("Depends", "Imports", "LinkingTo")),
      error = function(e) NULL
    )
    fields <- unlist(description, use.names = FALSE)
    fields <- fields[!is.na(fields)]
    if (length(fields) == 0L) {
      return(character())
    }
    names <- trimws(sub("\\s*\\(.*\\)$", "", unlist(strsplit(fields, ","))))
    setdiff(names[nzchar(names)], c("R", ""))
  }

  closure <- character()
  pending <- declared
  while (length(pending) > 0L) {
    current <- pending[1L]
    pending <- pending[-1L]
    if (current %in% closure) {
      next
    }
    closure <- c(closure, current)
    pending <- c(pending, setdiff(hard_deps(current), closure))
  }

  expect_false(any(c("Seurat", "reticulate", "concaveman", "V8", "sf") %in% closure))
})
