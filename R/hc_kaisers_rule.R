#' Determine number of components using Kaiser's rule
#'
#' `hc_kaisers_rule()` determines the number of principal components to retain
#' based on Kaiser's rule, which suggests keeping components with eigenvalues
#' greater than 1. If the cumulative explained variance at this point is less
#' than 80%, it can optionally suggest an alternative number of components that
#' achieves at least 80% explained variance.
#'
#' @param AnnDatR An AnnDatR object containing the data with PCA results.
#' @param with_alternative Logical indicating whether to suggest an alternative
#'   number of components if the explained variance at Kaiser's rule is below
#'   80% (default is `TRUE`).
#' @param verbose Logical indicating whether to report which rule was applied
#'   (default is `TRUE`).
#'
#' @returns Number of principal components to retain based on Kaiser's rule.
#'
#' @export
#' @examples
#' # Determine number of components using Kaiser's rule
#' adata_res <- hc_pca(example_adata, components = 40)
#' n_components <- hc_kaisers_rule(adata_res)
#' print(n_components)
hc_kaisers_rule <- function(AnnDatR, with_alternative = TRUE, verbose = TRUE) {
  pca_results <- AnnDatR[["uns"]][["pca"]]
  if (is.null(pca_results)) {
    stop(
      "AnnDatR$uns$pca not found. Call `hc_pca()` before `hc_kaisers_rule()`.",
      call. = FALSE
    )
  }

  # Eigenvalues of the correlation matrix are the squared component standard
  # deviations; Kaiser's rule keeps those above 1.
  eigenvalues <- pca_results[["sdev"]]^2
  n_comp <- which(eigenvalues < 1)[1]

  if (is.na(n_comp)) {
    stop(
      "No eigenvalue is lower than 1; Kaiser's rule cannot select a number of components.",
      call. = FALSE
    )
  }

  r2cum <- pca_results[["r2cum"]]

  if (r2cum[n_comp] >= 0.8) {
    if (verbose) {
      message("Kaiser's rule is above 80% variation. Success.")
    }
    return(n_comp)
  }

  if (!with_alternative) {
    if (verbose) {
      message("Suggested number of components explains less than 80% variation.")
    }
    return(n_comp)
  }

  alternative <- which(r2cum > 0.8)[1]
  if (is.na(alternative)) {
    stop(
      "No principal component achieves the cumulative R2 threshold of 0.8.",
      call. = FALSE
    )
  }

  if (verbose) {
    message(sprintf(
      paste0(
        "Explained variance at Kaiser's rule (component = %d) is under 80%%; ",
        "selected %d, which explains at least 80%% of the variation."
      ),
      n_comp,
      alternative
    ))
  }

  alternative
}
