#' Build a shared nearest neighbour graph from a distance matrix
#'
#' Pure-R implementation of the shared nearest neighbour construction that
#' `Seurat::FindNeighbors()` applies to a precomputed distance matrix. For each
#' observation the `k` nearest neighbours are taken (the observation itself is
#' always the first, at distance zero). Edge weights are the Jaccard index of
#' the two neighbour sets, `overlap / (2 * k - overlap)`, and edges weighing
#' less than `prune` are dropped.
#'
#' @param distances A `dist` object or square distance matrix.
#' @param k Number of nearest neighbours.
#' @param prune Edges below this weight are removed.
#'
#' @returns A list with the binary nearest neighbour graph `nn` and the
#'   weighted shared nearest neighbour graph `snn`, both sparse matrices.
#' @keywords internal
#' @noRd
build_snn <- function(distances, k, prune) {
  distance_matrix <- as.matrix(distances)
  n_obs <- nrow(distance_matrix)
  observation_names <- rownames(distance_matrix)

  if (k > n_obs) {
    warning(
      "`neighbors` is larger than the number of observations. Setting it to ",
      n_obs - 1L,
      ".",
      call. = FALSE
    )
    k <- n_obs - 1L
  }

  # Row i lists the k closest observations to i, itself included.
  nn_ranked <- t(apply(distance_matrix, 1L, function(row) order(row)[seq_len(k)]))

  nn_graph <- Matrix::sparseMatrix(
    i = rep(seq_len(n_obs), each = k),
    j = as.vector(t(nn_ranked)),
    x = 1,
    dims = c(n_obs, n_obs),
    dimnames = list(observation_names, observation_names)
  )

  # The number of shared neighbours between every pair is the inner product of
  # their neighbour indicator rows.
  snn_graph <- Matrix::tcrossprod(nn_graph)
  snn_graph <- methods::as(snn_graph, "generalMatrix")
  snn_graph@x <- snn_graph@x / (2 * k - snn_graph@x)
  snn_graph@x[snn_graph@x < prune] <- 0
  snn_graph <- Matrix::drop0(snn_graph)
  dimnames(snn_graph) <- list(observation_names, observation_names)

  list(nn = nn_graph, snn = snn_graph)
}

#' Compute Shared Nearest Neighbors (SNN) Graph
#'
#' `hc_snn()` computes a Shared Nearest Neighbors (SNN) graph from the distance
#' matrix stored in the AnnDatR object.
#'
#' @param AnnDatR An AnnDatR object containing the data with distance matrix
#'   results.
#' @param neighbors Number of nearest neighbors to consider (default is 20).
#' @param prune Pruning threshold for the SNN graph (default is 1/15). Edges
#'   with a Jaccard weight below this value are removed.
#' @param verbose Logical indicating whether to print progress messages
#'   (default is `TRUE`).
#'
#' @details
#' Neighbours are read directly from the distance matrix computed by
#' [hc_distance()], so the metric is the one chosen there. Earlier versions
#' accepted a `similarity` argument that was passed to `Seurat` as an Annoy
#' metric; it never had any effect on this code path, because Annoy is only
#' used when neighbours are searched in a coordinate space rather than taken
#' from a precomputed distance matrix. The argument has been removed.
#'
#' @returns An AnnDatR object with the nearest neighbour and shared nearest
#'   neighbour graphs stored in `uns$neighbors`.
#'
#' @export
#' @examples
#' # Compute SNN graph on an AnnDatR object after distance calculation
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res$uns$neighbors$snn[1:5, 1:5]
hc_snn <- function(
  AnnDatR,
  neighbors = 20,
  prune = 1 / 15,
  verbose = TRUE
) {
  distances <- AnnDatR[["uns"]][["distance"]]
  if (is.null(distances)) {
    stop(
      "AnnDatR$uns$distance not found. Call `hc_distance()` before `hc_snn()`.",
      call. = FALSE
    )
  }

  if (verbose) {
    message("Building SNN graph from the provided distance matrix.")
  }

  AnnDatR_out <- AnnDatR$clone(deep = TRUE)
  AnnDatR_out[["uns"]][["neighbors"]] <- build_snn(
    distances,
    k = neighbors,
    prune = prune
  )

  AnnDatR_out
}
