#' Centre and scale the columns of a matrix
#'
#' Base-R replacement for `pcaMethods::prep()`. Columns are mean-centred and
#' then divided by a per-column scale factor. Scale factors below `eps` are
#' replaced by 1 so that constant columns are left untouched rather than
#' producing `NaN`.
#'
#' @param x Numeric matrix.
#' @param scale One of `"none"`, `"uv"`, `"pareto"` or `"vector"`.
#' @param center Logical, whether to mean-centre the columns.
#' @param eps Lower bound below which a scale factor is treated as zero.
#'
#' @returns The centred and scaled matrix.
#' @keywords internal
#' @noRd
prep_matrix <- function(
  x,
  scale = c("none", "uv", "pareto", "vector"),
  center = TRUE,
  eps = 1e-12
) {
  scale <- match.arg(scale)
  x <- as.matrix(x)

  if (center) {
    x <- sweep(x, 2L, colMeans(x, na.rm = TRUE), "-")
  }

  scale_factors <- switch(
    scale,
    none = rep(1, ncol(x)),
    uv = apply(x, 2L, stats::sd, na.rm = TRUE),
    pareto = sqrt(apply(x, 2L, stats::sd, na.rm = TRUE)),
    vector = apply(x, 2L, function(column) sqrt(sum(column^2, na.rm = TRUE)))
  )

  below_eps <- scale_factors < eps
  if (any(below_eps)) {
    warning(
      sprintf(
        "Variance is below eps for %d variable(s). Not scaling them.",
        sum(below_eps)
      ),
      call. = FALSE
    )
    scale_factors[below_eps] <- 1
  }

  sweep(x, 2L, scale_factors, "/")
}

#' Scale a data matrix
#'
#' @param data Data to be scaled, with samples in rows and genes in columns.
#' @param method Scaling method: `"sample"` (default), `"gene"`, `"pareto"`,
#'   `"vector"`, `"min-max"`, `"max_sequencial"` or `"max"`.
#'
#' @returns Scaled data matrix with genes in rows and samples in columns.
#' @keywords internal
hc_scale <- function(data, method = "sample") {
  methods <- c(
    "sample",
    "gene",
    "pareto",
    "vector",
    "min-max",
    "max_sequencial",
    "max"
  )
  if (!method %in% methods) {
    stop(
      sprintf(
        "Scaling method not defined, only %s are allowed.",
        paste(sprintf('"%s"', methods), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  data <- as.matrix(data)

  switch(
    method,
    # Scale each gene (column) to unit variance, then transpose.
    sample = t(prep_matrix(
      data,
      scale = "uv",
      eps = .Machine[["double.eps"]]
    )),
    # Transpose first, so each sample (column) is scaled to unit variance.
    gene = prep_matrix(
      t(data),
      scale = "uv",
      eps = .Machine[["double.eps"]]
    ),
    pareto = prep_matrix(t(data), scale = "pareto"),
    vector = prep_matrix(t(data), scale = "vector"),
    `min-max` = {
      transposed <- t(data)
      row_min <- apply(transposed, 1L, min)
      row_max <- apply(transposed, 1L, max)
      (transposed - row_min) / (row_max - row_min)
    },
    max_sequencial = {
      transposed <- t(data)
      transposed / apply(transposed, 1L, max)
    },
    max = {
      transposed <- t(data)
      transposed / row_maxs(transposed)
    }
  )
}
