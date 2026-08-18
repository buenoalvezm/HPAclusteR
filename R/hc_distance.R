#' Compute distance matrix from PCA scores
#'
#' `hc_distance()` computes a distance matrix based on PCA scores stored in the
#' AnnDatR object.
#'
#' @param AnnDatR An AnnDatR object containing the data with PCA results.
#' @param components Number of principal components to be used. If `NULL`, it
#'   will be set to the first component with over 80% explained variance.
#' @param method Distance metric to use (default is `"spearman"`). This must be
#'   one of `"euclidean"`, `"maximum"`, `"manhattan"`, `"canberra"`,
#'   `"binary"`, `"minkowski"`, `"pearson"`, `"spearman"` or `"kendall"`.
#'
#' @returns An AnnDatR object with the distance matrix stored in `uns$distance`.
#'
#' @export
#' @examples
#' # Compute distance matrix on an AnnDatR object after PCA
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' head(adata_res$uns$distance)
hc_distance <- function(
  AnnDatR,
  components = NULL,
  method = "spearman"
) {
  pca_results <- AnnDatR[["uns"]][["pca"]]
  if (is.null(pca_results)) {
    stop(
      "AnnDatR$uns$pca not found. Call `hc_pca()` before `hc_distance()`.",
      call. = FALSE
    )
  }

  n_available <- ncol(pca_results[["scores"]])

  if (is.null(components)) {
    components <- which(pca_results[["r2cum"]] > 0.8)[1]
    if (is.na(components)) {
      stop(
        paste0(
          "No principal component reaches 80% explained variance. ",
          "Specify `components` explicitly, or compute more components in `hc_pca()`."
        ),
        call. = FALSE
      )
    }
    warning(
      "Number of components not specified. Setting components to ",
      components,
      ", the first component with over 80% explained variance.",
      call. = FALSE
    )
  } else if (components > n_available) {
    warning(
      "Number of components exceeds number of PCA components. Setting components to ",
      n_available,
      call. = FALSE
    )
    components <- n_available
  }

  distance <- hc_dist(
    pca_results[["scores"]][, seq_len(components), drop = FALSE],
    method = method
  )

  AnnDatR_out <- AnnDatR$clone(deep = TRUE)
  AnnDatR_out[["uns"]][["distance"]] <- distance

  AnnDatR_out
}
