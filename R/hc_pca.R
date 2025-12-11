#' Perform PCA on AnnDatR object
#'
#' `hc_pca()` performs Principal Component Analysis (PCA) on the data contained within an AnnDatR object.
#'
#' @param AnnDatR An AnnDatR object containing the data to perform PCA on.
#' @param components Number of principal components to compute (default is 50).
#' @param transform Data transformation method to apply before PCA. Options are "log1p", "sqrt", or "none" (default is "log1p").
#' @param scale_by Scaling method to apply to the data before PCA. Options include "sample", "gene", "pareto", "vector", "min-max", "max_sequencial", or "max" (default is "sample").
#' @param layer Optional layer name within the AnnDatR object to use for PCA. If NULL, the main data matrix is used (default is NULL).
#'
#' @returns PCA object and scores are stored within the AnnDatR object.
#'
#' @export
#' @examples
#' # Run PCA on an AnnDatR object
#' adata_res <- hc_pca(adata_t, components = 40)
#'
#' # PCA results
#' adata_res$obsm$X_pca
#'
#' # Explained variance
#' adata_res$uns$pca
hc_pca <- function(
  AnnDatR,
  components = 50,
  transform = "log1p",
  scale_by = "sample",
  layer = NULL
) {
  # Check obj structure
  if (!is.null(layer)) {
    wide_data <- AnnDatR$layers$layer
  } else {
    wide_data <- AnnDatR$X
  }

  # Check number of components
  if (components > dim(wide_data)[1]) {
    warning(
      "Number of components exceeds number of observations. Setting components to ",
      dim(wide_data)[1]
    )
    components <- dim(wide_data)[1]
  }

  # Apply the desired transformation
  if (transform == "log1p") {
    transformed_data <- wide_data |>
      dplyr::mutate_if(is.numeric, function(x) {
        log1p(x)
      }) |>
      tibble::column_to_rownames(colnames(wide_data)[1])
  } else if (transform == "sqrt") {
    transformed_data <- wide_data |>
      dplyr::mutate_if(is.numeric, function(x) {
        sqrt(x)
      }) |>
      tibble::column_to_rownames(colnames(wide_data)[1])
  } else if (transform == "none") {
    transformed_data <- wide_data |>
      tibble::column_to_rownames(colnames(wide_data)[1])
  } else {
    stop(
      "Transformation not defined, only 'log1p', 'sqrt', or 'none' options allowed"
    )
  }

  # Scale the data using the scaling function
  scaled_data <- hc_scale(transformed_data, method = scale_by)

  # Perform PCA
  pca_results <- pcaMethods::pca(scaled_data, nPcs = components)

  # Store the PCA results in AnnDatR
  AnnDatR[["uns"]][['pca']] <- pca_results
  AnnDatR[["obsm"]][[paste0("X_", "pca")]] <- pcaMethods::scores(AnnDatR[[
    "uns"
  ]][['pca']])

  return(AnnDatR)
}
