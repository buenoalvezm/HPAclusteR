#' Compute Shared Nearest Neighbors (SNN) Graph
#'
#' `hc_snn()` computes a Shared Nearest Neighbors (SNN) graph from the distance matrix stored in the AnnDatR object.
#'
#' @param AnnDatR An AnnDatR object containing the data with distance matrix results.
#' @param neighbors Number of nearest neighbors to consider (default is 20).
#' @param prune Pruning threshold for the SNN graph (default is 1/15).
#' @param similarity Distance metric to use for nearest neighbor search (default is "euclidean"). Other options include "cosine", "manhattan", and "hamming".
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#'
#' @returns SNN graph stored within the AnnDatR object.
#'
#' @export
#' @examples
#' # Compute SNN graph on an AnnDatR object after distance calculation
#' adata_pca <- hc_pca(adata_t, components = 40)
#' adata_dist <- hc_distance(adata_pca, components = 20)
#' adata_snn <- hc_snn(adata_dist, neighbors = 15)
#' adata_snn$uns$neighbors$snn
hc_snn <- function(
  AnnDatR,
  neighbors = 20,
  prune = 1 / 15,
  similarity = "euclidean",
  verbose = TRUE
) {
  if (is.null(AnnDatR[["uns"]][["distance"]])) {
    stop(
      "AnnDatR$uns$distance not found. Call `hc_distance()` before `hc_snn()`."
    )
  }
  dist_results <- AnnDatR[["uns"]][["distance"]]

  AnnDatR[["uns"]][["neighbors"]] <- dist_results |>
    as.matrix() |>
    stats::as.dist() |>
    Seurat::FindNeighbors(
      k.param = neighbors,
      compute.SNN = TRUE,
      prune.SNN = prune,
      nn.method = "annoy",
      annoy.metric = similarity,
      nn.eps = 0,
      verbose = verbose,
      force.recalc = FALSE
    )

  return(AnnDatR)
}
