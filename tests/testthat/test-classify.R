test_that("hc_classify_fun assigns the documented specificity categories", {
  # Three genes with deliberately chosen profiles across six tissues.
  X <- tibble::tibble(
    sample_id = paste0("s", 1:6),
    enriched = c(100, 1, 1, 1, 1, 1),
    low_spec = c(10, 10, 10, 10, 10, 10),
    silent = c(0, 0, 0, 0, 0, 0)
  )
  obs <- tibble::tibble(ensembl_id = c("enriched", "low_spec", "silent"))
  var <- tibble::tibble(
    sample_id = paste0("s", 1:6),
    tissue_name = paste0("t", 1:6)
  )
  adata <- AnnDatR$new(X = X, obs = obs, var = var, var_names = "sample_id")

  result <- hc_classify_fun(adata, "tissue_name")

  categories <- stats::setNames(result[["spec_category"]], result[["ENSG"]])
  expect_equal(unname(categories[["enriched"]]), "Enriched")
  expect_equal(unname(categories[["low_spec"]]), "Low specificity")
  expect_equal(unname(categories[["silent"]]), "Not detected")

  specific <- stats::setNames(result[["spec_sample_categories"]], result[["ENSG"]])
  expect_equal(unname(specific[["enriched"]]), "t1")
})

test_that("hc_classify_fun assigns the documented distribution categories", {
  # Eight samples, so that two detected is 25% -- under the 31% "many" cut-off.
  X <- tibble::tibble(
    sample_id = paste0("s", 1:8),
    all_of_them = rep(10, 8),
    just_one = c(50, rep(0, 7)),
    a_few = c(50, 40, rep(0, 6)),
    many = c(rep(50, 3), rep(0, 5)),
    none = rep(0, 8)
  )
  obs <- tibble::tibble(
    ensembl_id = c("all_of_them", "just_one", "a_few", "many", "none")
  )
  var <- tibble::tibble(
    sample_id = paste0("s", 1:8),
    tissue_name = paste0("t", 1:8)
  )
  adata <- AnnDatR$new(X = X, obs = obs, var = var, var_names = "sample_id")

  result <- hc_classify_fun(adata, "tissue_name")
  distribution <- stats::setNames(result[["dist_category"]], result[["ENSG"]])

  expect_equal(unname(distribution[["all_of_them"]]), "Detected in all")
  expect_equal(unname(distribution[["just_one"]]), "Detected in single")
  expect_equal(unname(distribution[["a_few"]]), "Detected in some")
  expect_equal(unname(distribution[["many"]]), "Detected in many")
  expect_equal(unname(distribution[["none"]]), "Not detected")
})

test_that("hc_classify_fun computes tau between 0 and 1", {
  X <- tibble::tibble(
    sample_id = paste0("s", 1:4),
    specific = c(100, 0, 0, 0),
    uniform = c(10, 10, 10, 10)
  )
  obs <- tibble::tibble(ensembl_id = c("specific", "uniform"))
  var <- tibble::tibble(sample_id = paste0("s", 1:4), tissue_name = paste0("t", 1:4))
  adata <- AnnDatR$new(X = X, obs = obs, var = var, var_names = "sample_id")

  result <- hc_classify_fun(adata, "tissue_name")
  tau <- stats::setNames(result[["tau"]], result[["ENSG"]])

  expect_true(all(result[["tau"]] >= 0 & result[["tau"]] <= 1))
  expect_equal(unname(tau[["uniform"]]), 0)
  expect_gt(tau[["specific"]], 0.9)
})

test_that("hc_classify_fun returns one row per gene with the expected columns", {
  adata <- toy_adata(n_genes = 20L, n_samples = 9L)
  result <- hc_classify_fun(adata, "tissue_name")

  expect_equal(nrow(result), 20L)
  expect_setequal(
    names(result),
    c("ENSG", "spec_category", "spec_sample_categories", "tau", "dist_category")
  )
  expect_false(anyNA(result[["spec_category"]]))
  expect_false(anyNA(result[["dist_category"]]))
})

test_that("hc_classify splits results by cluster", {
  result <- hc_classify(adata_clustered(), "tissue_name")

  expect_setequal(names(result), c("classification", "pie_charts"))
  expect_equal(names(result[["classification"]]), names(result[["pie_charts"]]))

  clusters <- unique(adata_clustered()[["obs"]][["cluster"]])
  expect_setequal(names(result[["classification"]]), as.character(clusters))

  # Every gene appears in exactly one cluster's classification.
  classified <- unlist(lapply(result[["classification"]], `[[`, "ENSG"))
  expect_equal(anyDuplicated(classified), 0L)

  expect_s3_class(result[["classification"]][[1]], "tbl_df")
  expect_s3_class(result[["pie_charts"]][[1]], "ggplot")
})

test_that("plot_specificity_distribution renders with missing categories", {
  minimal <- tibble::tibble(
    ENSG = c("a", "b"),
    spec_category = c("Enriched", "Enriched"),
    spec_sample_categories = c("t1", "t1"),
    tau = c(0.9, 0.8),
    dist_category = c("Detected in all", "Detected in all")
  )

  plot <- plot_specificity_distribution(minimal)
  expect_s3_class(plot, "ggplot")
  expect_no_error(ggplot2::ggplot_build(plot))
})
