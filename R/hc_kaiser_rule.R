#' Determine number of components using Kaiser's rule
#'
#' `hc_kaisers_rule()` determines the number of principal components to retain based on Kaiser's rule, which suggests keeping components with eigenvalues greater than 1. If the cumulative explained variance at this point is less than 80%, it can optionally suggest an alternative number of components that achieves at least 80% explained variance.
#'
#' @param AnnDatR An AnnDatR object containing the data with PCA results.
#' @param with_alternative Logical indicating whether to suggest an alternative number of components if the explained variance at Kaiser's rule is below 80% (default is TRUE).
#'
#' @returns Number of principal components to retain based on Kaiser's rule.
#'
#' @export
#' @examples
#' # Determine number of components using Kaiser's rule
#' adata_res <- hc_pca(adata_t, components = 40)
#' n_components <- hc_kaisers_rule(adata_res)
#' print(n_components)
hc_kaisers_rule <- function(AnnDatR, with_alternative = TRUE) {
  if (is.null(AnnDatR[["uns"]][["pca"]])) {
    stop(
      "AnnDatR$uns$pca not found. Call `hc_pca()` before `hc_kaisers_rule()`."
    )
  }
  pca_results <- AnnDatR[["uns"]][["pca"]]

  # Extract the squared standard deviations (eigenvalues)
  squared_devs <- pcaMethods::sDev(pca_results)^2

  # Find the first component where the squared standard deviation (eigenvalue) is less than 1
  n_comp <- which(squared_devs < 1)[1]

  # Check if n_comp is NA or exceeds the number of components
  if (is.na(n_comp) || n_comp > length(squared_devs)) {
    stop("No eigenvalue is lower than 1")
  }

  # Check cumulative R2 value at the identified component
  if (pca_results@R2cum[n_comp] < 0.8) {
    if (with_alternative) {
      print(
        paste0(
          "Explained variance at Kaiser rule (component = ",
          n_comp,
          ") is under 80%, selected ",
          which(pca_results@R2cum > 0.8)[1],
          " which explains at least 80% Variation"
        )
      )
      n_comp <- which(pca_results@R2cum > 0.8)[1]
      # If no component satisfies the 0.8 threshold
      if (is.na(n_comp)) {
        stop(
          "No principal component achieves the cumulative R2 threshold of 0.8"
        )
      }
    } else {
      print("Suggested number of components explains less than 80% variation")
    }
  } else {
    print("Kaiser\'s rule is above 80% variation. Success")
  }

  return(n_comp)
}
