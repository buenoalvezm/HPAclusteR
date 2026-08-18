#' Default value for `NULL`
#'
#' @param x Value to test.
#' @param y Fallback used when `x` is `NULL`.
#'
#' @returns `x` unless it is `NULL`, in which case `y`.
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Require a suggested package
#'
#' Fails with a consistent, actionable message when an optional dependency is
#' missing.
#'
#' @param pkg Package name.
#' @param what Short description of what the package is needed for.
#' @param bioc Logical, whether the package lives on Bioconductor.
#'
#' @returns `TRUE` invisibly if the package is available, otherwise an error.
#' @keywords internal
#' @noRd
check_installed <- function(pkg, what = NULL, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }

  install_call <- if (bioc) {
    sprintf('BiocManager::install("%s")', pkg)
  } else {
    sprintf('install.packages("%s")', pkg)
  }

  stop(
    sprintf(
      "The '%s' package is required%s. Install it with `%s`.",
      pkg,
      if (is.null(what)) "" else paste0(" ", what),
      install_call
    ),
    call. = FALSE
  )
}

#' Compute a distance matrix
#'
#' Base-R replacement for `factoextra::get_dist()`. Correlation-based metrics
#' are converted to distances as `1 - r`; all other metrics are forwarded to
#' [stats::dist()].
#'
#' @param x Numeric matrix or data frame with observations in rows.
#' @param method One of `"euclidean"`, `"maximum"`, `"manhattan"`,
#'   `"canberra"`, `"binary"`, `"minkowski"`, `"pearson"`, `"spearman"` or
#'   `"kendall"`.
#' @param ... Further arguments passed to [stats::dist()], e.g. `p` for
#'   Minkowski distance.
#'
#' @returns An object of class `dist`.
#' @keywords internal
#' @noRd
hc_dist <- function(x, method = "euclidean", ...) {
  correlation_methods <- c("pearson", "spearman", "kendall")
  stats_methods <- c(
    "euclidean",
    "maximum",
    "manhattan",
    "canberra",
    "binary",
    "minkowski"
  )
  method <- match.arg(method, c(stats_methods, correlation_methods))

  x <- as.matrix(x)

  if (method %in% correlation_methods) {
    correlations <- stats::cor(t(x), method = method)
    return(stats::as.dist(1 - correlations))
  }

  stats::dist(x, method = method, ...)
}

#' Row-wise maxima
#'
#' Base-R replacement for `matrixStats::rowMaxs()`.
#'
#' @param x Numeric matrix.
#' @param na.rm Logical, whether to drop missing values.
#'
#' @returns Numeric vector of row maxima.
#' @keywords internal
#' @noRd
row_maxs <- function(x, na.rm = FALSE) {
  do.call(pmax, c(lapply(seq_len(ncol(x)), function(j) x[, j]), na.rm = na.rm))
}
