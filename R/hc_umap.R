#' Check if version v1 is greater than or equal to version v2
#'
#' @param v1 version string 1
#' @param v2 version string 2
#'
#' @returns TRUE if v1 >= v2, FALSE otherwise
#' @keywords internal
version_gte <- function(v1, v2) {
  v1 <- as.integer(strsplit(v1, "\\.")[[1]])
  v2 <- as.integer(strsplit(v2, "\\.")[[1]])
  length(v1) <- length(v2) <- max(length(v1), length(v2))
  any((v1 > v2)[cumsum(v1 != v2) == 1])
}

#' Create UMAP embeddings from SNN graph
#'
#' `hc_umap()` generates UMAP embeddings based on the Shared Nearest Neighbors (SNN) graph stored in the AnnDatR object.
#'
#' @param AnnDatR An AnnDatR object containing the data with SNN graph results.
#' @param n_epochs Number of epochs for UMAP optimization. If NULL, it will be set to 200 for datasets with more than 10,000 cells and 500 otherwise (default is NULL).
#' @param seed Random seed for reproducibility (default is 42).
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#'
#' @returns UMAP embeddings stored within the AnnDatR object.
#'
#' @export
#' @examples
#' # Generate UMAP embeddings on an AnnDatR object after SNN graph computation
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_umap(adata_res)
#' head(adata_res$obsm$X_umap)
hc_umap <- function(AnnDatR, n_epochs = NULL, seed = 42, verbose = TRUE) {
  if (is.null(AnnDatR[["uns"]][["neighbors"]])) {
    stop(
      "AnnDatR$uns$neighbors not found. Call `hc_snn()` before `hc_umap()`."
    )
  }

  n_cells <- nrow(AnnDatR[["uns"]][["neighbors"]][["snn"]])
  if (is.null(n_epochs)) {
    n_epochs <- ifelse(n_cells > 10000, 200, 500)
  }

  # Check if umap-learn is available
  if (!reticulate::py_module_available("umap")) {
    message(
      "The 'umap-learn' Python package is not installed. Installing it now..."
    )
    reticulate::py_install("umap-learn")
  }
  # Check numpy version and downgrade if needed
  np <- reticulate::import("numpy", delay_load = TRUE)
  np_version <- as.character(np$`__version__`)
  if (version_gte(np_version, "2.3.0")) {
    message("Downgrading numpy to 2.3.0 for compatibility...")
    reticulate::py_install("numpy==2.3.0", pip = TRUE)
  }

  # Import the umap module and fix issue with version string
  umap <- reticulate::import("umap", delay_load = TRUE)
  umap$`__version__` <- gsub("[^0-9.]", "", umap$`__version__`) # Remove non-numeric parts

  set.seed(seed)

  # exact piping and extraction as in original:
  umap_raw_mat <- AnnDatR[["uns"]][["neighbors"]][["snn"]] |>
    Seurat::RunUMAP(
      umap.method = "umap-learn",
      n.epochs = n_epochs,
      seed.use = seed,
      assay = "RNA",
      verbose = verbose
    ) |>
    (\(x) x@cell.embeddings)()

  if (is.null(rownames(umap_raw_mat))) {
    rownames(umap_raw_mat) <- rownames(AnnDatR[["uns"]][["neighbors"]][["snn"]])
  }

  AnnDatR_out <- AnnDatR$clone(deep = TRUE)

  AnnDatR_out[["obsm"]][["X_umap_raw"]] <- umap_raw_mat
  AnnDatR_out[["obsm"]][["X_umap"]] <-
    umap_raw_mat |>
    tibble::as_tibble(rownames = "gene") |>
    tidyr::gather(
      !!rlang::sym("UMAP"),
      !!rlang::sym("UMAP_value"),
      !!rlang::sym("UMAP_1"),
      !!rlang::sym("UMAP_2")
    ) |>
    dplyr::mutate(
      UMAP_value = (!!rlang::sym("UMAP_value") -
        min(!!rlang::sym("UMAP_value"))) /
        (max(!!rlang::sym("UMAP_value")) - min(!!rlang::sym("UMAP_value"))) *
        2 -
        1
    ) |>
    dplyr::group_by(!!rlang::sym("UMAP")) |>
    dplyr::mutate(
      UMAP_value = !!rlang::sym("UMAP_value") -
        mean(range(!!rlang::sym("UMAP_value")))
    ) |>
    dplyr::ungroup() |>
    tidyr::spread(!!rlang::sym("UMAP"), !!rlang::sym("UMAP_value")) |>
    tibble::column_to_rownames("gene") |>
    as.matrix()

  AnnDatR_out[["obs"]] <- AnnDatR_out[["obs"]] |>
    dplyr::left_join(
      AnnDatR_out[["obsm"]][["X_umap"]] |>
        tibble::as_tibble() |>
        (\(x) {
          colnames(x) <- paste0("UMAP", 1:ncol(x))
          x
        })() |>
        dplyr::mutate(ensembl_id = rownames(AnnDatR_out[["obsm"]][["X_umap"]])),
      by = dplyr::join_by(!!rlang::sym("ensembl_id"))
    )

  return(AnnDatR_out)
}
