#' One-call clustering pipeline
#'
#' `hc_auto_cluster` performs a full clustering pipeline on an AnnDataR object,
#' including PCA, distance calculation, SNN graph construction, consensus clustering,
#' UMAP embedding, and cluster hull generation. It prints a UMAP plot of the
#' resulting clusters.
#'
#' @param AnnDatR An AnnDataR object.
#' @param layer Optional layer name within the AnnDatR object to use instead of the default data.
#' @param pca_components Number of principal components to compute (default is 50).
#' @param pca_transform Data transformation method to apply before PCA. Options are "log1p", "sqrt", or "none" (default is "log1p").
#' @param distance_components Number of principal components to be used. If NULL, uses Kaiser’s rule to determine the number of components.
#' @param distance_method Distance metric to use (default is 'spearman'). This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
#' @param snn_neighbors Number of nearest neighbors to consider (default is 20).
#' @param snn_prune Pruning threshold for the SNN graph (default is 1/15).
#' @param cluster_resolution Resolution parameter for clustering (default is 6).
#' @param cluster_method Clustering method to use: "louvain" (default) or "leiden".
#' @param cluster_seeds Number of random seeds to aggregate over in the
#'   consensus step (default is 100).
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#'
#' @returns An AnnDataR object with clustering and intermediate results added. The UMAP plot of the clusters is also printed.
#'
#' @export
#' @examples
#' adata_res <- hc_auto_cluster(
#'   example_adata,
#'   cluster_resolution = 8,
#'   cluster_seeds = 20
#' )
#' head(adata_res$uns$consensus_clustering)
hc_auto_cluster <- function(
  AnnDatR,
  layer = NULL,
  pca_components = 50,
  pca_transform = "log1p",
  distance_components = NULL,
  distance_method = "spearman",
  snn_neighbors = 20,
  snn_prune = 1 / 15,
  cluster_resolution = 6,
  cluster_method = "louvain",
  cluster_seeds = 100,
  verbose = TRUE
) {
  AnnDatR_res <- AnnDatR |>
    hc_pca(
      components = pca_components,
      transform = pca_transform,
      scale_by = "sample",
      layer = layer,
      verbose = verbose
    )

  AnnDatR_res <- AnnDatR_res |>
    hc_distance(
      components = distance_components %||%
        hc_kaisers_rule(AnnDatR_res, verbose = verbose),
      method = distance_method
    ) |>
    hc_snn(
      neighbors = snn_neighbors,
      prune = snn_prune,
      verbose = verbose
    ) |>
    hc_cluster_consensus(
      resolution = cluster_resolution,
      method = cluster_method,
      n_seeds = cluster_seeds,
      verbose = verbose
    ) |>
    hc_umap(verbose = verbose) |>
    hc_cluster_hulls()

  print(hc_plot_umap(AnnDatR_res))

  AnnDatR_res
}
