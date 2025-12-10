#' Compute distance matrix from PCA scores
#'
#' `hc_do_distance()` computes a distance matrix based on PCA scores stored in the AnnDatR object.
#'
#' @param AnnDatR An AnnDatR object containing the data with PCA results.
#' @param components Number of principal components to be used. If NULL, it will be set to the first component with over 80% explained variance.
#' @param method Distance metric to use (default is 'spearman'). This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman" or "kendall".
#'
#' @returns Distance matrix stored within the AnnDatR object.
#'
#' @export
#' @examples
#' # Compute distance matrix on an AnnDatR object after PCA
#' adata_pca <- hc_do_pca(adata_t, components = 40)
#' adata_dist <- hc_do_distance(adata_pca, components = 20)
#' head(adata_dist$uns$distance)
hc_do_distance <- function(
  AnnDatR,
  components = NULL,
  method = "spearman"
) {
  if (is.null(AnnDatR[["uns"]][["pca"]])) {
    stop(
      "AnnDatR$uns$pca not found. Call `hc_do_pca()` before `hc_do_distance()`."
    )
  }
  pca_results <- AnnDatR[["uns"]][["pca"]]

  if (is.null(components)) {
    warning(
      "Number of components not specified. Setting components to first component with over 80% explained variance."
    )
    components <- which(pca_results@R2cum > 0.8)[1]
  } else if (components > pca_results@nPcs) {
    warning(
      "Number of components exceeds number of PCA components. Setting components to ",
      pca_results@nPcs
    )
    components <- pca_results@nPcs
  }

  if (is.null(AnnDatR[["uns"]][["distance"]])) {
    AnnDatR[["uns"]][["distance"]] = list()
  }

  distance <- pca_results@scores[, 1:components] |>
    factoextra::get_dist(method = method)

  AnnDatR[["uns"]][["distance"]] <- distance

  return(AnnDatR)
}
