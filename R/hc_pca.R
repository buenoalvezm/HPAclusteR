#' Principal component analysis by singular value decomposition
#'
#' Base-R replacement for `pcaMethods::pca(method = "svd")`. Like `pcaMethods`,
#' the matrix is mean-centred (but not rescaled) before decomposition.
#'
#' @param x Numeric matrix without missing values.
#' @param n_pcs Number of components to return.
#'
#' @returns A list with `scores`, `loadings`, `sdev`, `r2cum` and `method`.
#' @keywords internal
#' @noRd
pca_svd <- function(x, n_pcs) {
  decomposition <- stats::prcomp(x, center = TRUE, scale. = FALSE)
  n_pcs <- min(n_pcs, ncol(decomposition[["rotation"]]))
  keep <- seq_len(n_pcs)

  variances <- decomposition[["sdev"]]^2
  scores <- decomposition[["x"]][, keep, drop = FALSE]

  list(
    scores = scores,
    loadings = decomposition[["rotation"]][, keep, drop = FALSE],
    sdev = apply(scores, 2L, stats::sd),
    r2cum = cumsum(variances)[keep] / sum(variances),
    method = "svd"
  )
}

#' Impute missing values with column means
#'
#' @param x Numeric matrix.
#'
#' @returns The matrix with `NA` replaced by the mean of their column. Columns
#'   that are entirely missing are set to zero.
#' @keywords internal
#' @noRd
impute_column_means <- function(x) {
  missing <- which(is.na(x), arr.ind = TRUE)
  if (nrow(missing) == 0L) {
    return(x)
  }
  column_means <- colMeans(x, na.rm = TRUE)
  column_means[is.nan(column_means)] <- 0
  x[missing] <- column_means[missing[, "col"]]
  x
}

#' Perform PCA on an AnnDatR object
#'
#' `hc_pca()` performs Principal Component Analysis (PCA) on the data contained
#' within an AnnDatR object.
#'
#' @param AnnDatR An AnnDatR object containing the data to perform PCA on.
#' @param components Number of principal components to compute (default is 50).
#' @param transform Data transformation method to apply before PCA. Options are
#'   `"log1p"`, `"sqrt"` or `"none"` (default is `"log1p"`).
#' @param scale_by Scaling method to apply to the data before PCA. Options
#'   include `"sample"`, `"gene"`, `"pareto"`, `"vector"`, `"min-max"`,
#'   `"max_sequencial"` or `"max"` (default is `"sample"`).
#' @param layer Optional layer name within the AnnDatR object to use for PCA.
#'   If `NULL`, the main data matrix is used (default is `NULL`).
#' @param na_action How to deal with missing values in the expression matrix.
#'   One of `"impute"` (default, replace missing values by their column mean and
#'   run a fast SVD), `"nipals"` (iterative PCA that handles missing values
#'   directly; exact but much slower, and requires the `pcaMethods` package),
#'   `"omit"` (drop genes with any missing value) or `"error"` (stop).
#' @param verbose Logical indicating whether to report on missing values
#'   (default is `TRUE`).
#'
#' @details
#' Missing values are the main performance trap in this step. Earlier versions
#' of this package delegated to `pcaMethods::pca()`, which silently fell back to
#' the iterative NIPALS algorithm whenever it encountered an `NA`. NIPALS is
#' orders of magnitude slower than an SVD on a complete matrix, and nothing in
#' the output indicated why the call had become slow.
#'
#' `hc_pca()` now reports how many values are missing and, by default, imputes
#' them with their column mean so that the fast SVD path is taken. Set
#' `na_action = "nipals"` to reproduce the previous behaviour exactly.
#'
#' @returns An AnnDatR object with the PCA fit stored in `uns$pca` (a list with
#'   `scores`, `loadings`, `sdev`, `r2cum` and `method`) and the scores in
#'   `obsm$X_pca`.
#'
#' @export
#' @examples
#' # Run PCA on an AnnDatR object
#' adata_res <- hc_pca(example_adata, components = 40)
#'
#' # PCA scores
#' head(adata_res$obsm$X_pca[, 1:5])
#'
#' # Cumulative explained variance
#' head(adata_res$uns$pca$r2cum)
hc_pca <- function(
  AnnDatR,
  components = 50,
  transform = c("log1p", "sqrt", "none"),
  scale_by = "sample",
  layer = NULL,
  na_action = c("impute", "nipals", "omit", "error"),
  verbose = TRUE
) {
  transform <- match.arg(transform)
  na_action <- match.arg(na_action)

  wide_data <- if (is.null(layer)) {
    AnnDatR[["X"]]
  } else {
    if (is.null(AnnDatR[["layers"]][[layer]])) {
      stop(
        sprintf("Layer '%s' not found in AnnDatR$layers.", layer),
        call. = FALSE
      )
    }
    AnnDatR[["layers"]][[layer]]
  }

  if (components > nrow(wide_data)) {
    warning(
      "Number of components exceeds number of observations. Setting components to ",
      nrow(wide_data),
      call. = FALSE
    )
    components <- nrow(wide_data)
  }

  id_column <- colnames(wide_data)[1]
  numeric_data <- as.matrix(wide_data[, -1, drop = FALSE])
  rownames(numeric_data) <- wide_data[[id_column]]

  numeric_data <- switch(
    transform,
    log1p = log1p(numeric_data),
    sqrt = sqrt(numeric_data),
    none = numeric_data
  )

  # Genes in rows, samples in columns after scaling.
  scaled_data <- hc_scale(numeric_data, method = scale_by)

  n_missing <- sum(is.na(scaled_data))
  if (n_missing > 0L && verbose) {
    message(sprintf(
      "%d of %d values (%.2f%%) are missing; using na_action = '%s'.",
      n_missing,
      length(scaled_data),
      100 * n_missing / length(scaled_data),
      na_action
    ))
  }

  if (n_missing > 0L && na_action == "error") {
    stop(
      sprintf(
        paste0(
          "%d missing values found in the expression matrix. ",
          "Set `na_action` to \"impute\", \"nipals\" or \"omit\" to handle them."
        ),
        n_missing
      ),
      call. = FALSE
    )
  }

  if (n_missing > 0L && na_action == "nipals") {
    check_installed(
      "pcaMethods",
      "for `na_action = \"nipals\"`",
      bioc = TRUE
    )
    if (verbose) {
      message(
        "Running NIPALS; this is much slower than the default `na_action = \"impute\"`."
      )
    }
    components <- min(components, nrow(scaled_data), ncol(scaled_data))
    fit <- pcaMethods::pca(scaled_data, nPcs = components, method = "nipals")
    pca_results <- list(
      scores = pcaMethods::scores(fit),
      loadings = pcaMethods::loadings(fit),
      sdev = pcaMethods::sDev(fit),
      r2cum = fit@R2cum,
      method = "nipals"
    )
  } else {
    if (n_missing > 0L) {
      scaled_data <- switch(
        na_action,
        impute = impute_column_means(scaled_data),
        omit = {
          complete <- stats::complete.cases(scaled_data)
          if (verbose) {
            message(sprintf(
              "Dropping %d of %d genes with missing values.",
              sum(!complete),
              length(complete)
            ))
          }
          scaled_data[complete, , drop = FALSE]
        }
      )
    }
    components <- min(components, nrow(scaled_data), ncol(scaled_data))
    pca_results <- pca_svd(scaled_data, n_pcs = components)
  }

  AnnDatR_out <- AnnDatR$clone(deep = TRUE)
  AnnDatR_out[["uns"]][["pca"]] <- pca_results
  AnnDatR_out[["obsm"]][["X_pca"]] <- pca_results[["scores"]]

  AnnDatR_out
}
