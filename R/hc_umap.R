#' Rescale UMAP coordinates to a centred, symmetric range
#'
#' Both axes are mapped jointly onto `[-1, 1]` so that the aspect ratio of the
#' embedding is preserved, then each axis is centred on its own midpoint.
#'
#' @param embedding Two-column numeric matrix of UMAP coordinates.
#'
#' @returns A matrix with the same dimensions and row names.
#' @keywords internal
#' @noRd
rescale_umap <- function(embedding) {
  overall_range <- range(embedding, na.rm = TRUE)
  rescaled <- (embedding - overall_range[1]) /
    (overall_range[2] - overall_range[1]) *
    2 -
    1

  # Centre each axis on the midpoint of its own range.
  axis_midpoints <- apply(rescaled, 2L, function(axis) mean(range(axis)))
  sweep(rescaled, 2L, axis_midpoints, "-")
}

#' Create UMAP embeddings from SNN graph
#'
#' `hc_umap()` generates UMAP embeddings based on the Shared Nearest Neighbors
#' (SNN) graph stored in the AnnDatR object.
#'
#' @param AnnDatR An AnnDatR object containing the data with SNN graph results.
#' @param n_epochs Number of epochs for UMAP optimization. If `NULL`, it will be
#'   set to 200 for datasets with more than 10,000 observations and 500
#'   otherwise (default is `NULL`).
#' @param min_dist Minimum distance between embedded points (default is 0.3).
#' @param spread Effective scale of the embedded points (default is 1).
#' @param seed Random seed for reproducibility (default is 42).
#' @param verbose Logical indicating whether to print progress messages
#'   (default is `TRUE`).
#'
#' @details
#' The embedding is produced by [uwot::optimize_graph_layout()], a native R and
#' C++ implementation of UMAP's simplicial set embedding. Earlier versions
#' called the Python `umap-learn` package through `reticulate` by way of
#' `Seurat::RunUMAP()`, which required a working Python installation with
#' pinned `numpy` and `umap-learn` versions and was the main source of
#' installation and continuous integration failures.
#'
#' The algorithm and its parameters are the same, but the two implementations
#' use different random number generators and optimisers, so embeddings are
#' equivalent in structure rather than numerically identical to those produced
#' by earlier versions.
#'
#' @returns An AnnDatR object with the raw embedding in `obsm$X_umap_raw`, the
#'   rescaled embedding in `obsm$X_umap`, and `UMAP1`/`UMAP2` columns added to
#'   `obs`.
#'
#' @export
#' @examples
#' # Generate UMAP embeddings on an AnnDatR object after SNN graph computation
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_umap(adata_res)
#' head(adata_res$obsm$X_umap)
hc_umap <- function(
  AnnDatR,
  n_epochs = NULL,
  min_dist = 0.3,
  spread = 1,
  seed = 42,
  verbose = TRUE
) {
  check_installed("uwot", "to compute UMAP embeddings")

  if (is.null(AnnDatR[["uns"]][["neighbors"]])) {
    stop(
      "AnnDatR$uns$neighbors not found. Call `hc_snn()` before `hc_umap()`.",
      call. = FALSE
    )
  }

  snn_graph <- AnnDatR[["uns"]][["neighbors"]][["snn"]]
  n_observations <- nrow(snn_graph)
  n_epochs <- n_epochs %||% if (n_observations > 10000) 200L else 500L

  graph <- snn_graph
  Matrix::diag(graph) <- 0
  graph <- Matrix::drop0(graph)

  set.seed(seed)
  embedding <- uwot::optimize_graph_layout(
    graph = graph,
    n_components = 2L,
    n_epochs = n_epochs,
    learning_rate = 1,
    init = "spectral",
    min_dist = min_dist,
    spread = spread,
    repulsion_strength = 1,
    negative_sample_rate = 5,
    verbose = verbose
  )

  embedding <- as.matrix(embedding)

  # Mean-centre each axis, as `Seurat::RunUMAP()` did before returning. uwot
  # already returns a centred embedding, so this is normally a no-op, but it
  # keeps the downstream rescaling independent of that implementation detail.
  embedding <- sweep(embedding, 2L, colMeans(embedding), "-")

  rownames(embedding) <- rownames(snn_graph)
  colnames(embedding) <- c("UMAP_1", "UMAP_2")

  AnnDatR_out <- AnnDatR$clone(deep = TRUE)
  AnnDatR_out[["obsm"]][["X_umap_raw"]] <- embedding
  AnnDatR_out[["obsm"]][["X_umap"]] <- rescale_umap(embedding)

  gene_column <- AnnDatR_out[["obs_names_col"]]
  rescaled <- AnnDatR_out[["obsm"]][["X_umap"]]
  umap_coordinates <- tibble::tibble(
    gene = rownames(rescaled),
    UMAP1 = rescaled[, 1],
    UMAP2 = rescaled[, 2]
  )
  names(umap_coordinates)[1] <- gene_column

  AnnDatR_out[["obs"]] <- AnnDatR_out[["obs"]] |>
    dplyr::select(-dplyr::any_of(c("UMAP1", "UMAP2"))) |>
    dplyr::left_join(umap_coordinates, by = gene_column)

  AnnDatR_out
}
