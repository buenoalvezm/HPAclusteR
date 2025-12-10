#' Title
#'
#' @param pca_results
#' @param with_alternative
#'
#' @returns
#'
#' @export
#' @examples
hc_do_kaisers_rule <- function(pca_results, with_alternative = TRUE) {
  if (class(pca_results) == 'prcomp') {
    return(kaisers_PCA_rule_prcomp(pca_results, with_alternative))
  }
  # Extract the squared standard deviations (eigenvalues)
  squared_devs <- pcaMethods::sDev(pca_results)^2

  # Find the first component where the squared standard deviation (eigenvalue) is less than 1
  n_comp <- which(squared_devs < 1)[1]

  # Check if n_comp is NA or exceeds the number of components
  if (is.na(n_comp) || n_comp > length(squared_devs)) {
    stop("No eigenvalue is lower than 1")
  }

  # Check cumulative R² value at the identified component
  if (pca_results@R2cum[n_comp] < 0.8) {
    if (with_alternative) {
      print(
        'Explained variance at Kaiser rule is under 80%, suggesting at least 80% Variation'
      )
      n_comp <- which(pca_results@R2cum > 0.8)[1]
      # If no component satisfies the 0.8 threshold
      if (is.na(n_comp)) {
        stop(
          "No principal component achieves the cumulative R² threshold of 0.8"
        )
      }
    } else {
      print('Suggested n_comp explains less than 80% variation')
    }
  } else {
    print('Kaiser\'s rule is also above 80% variation')
  }

  return(n_comp)
}
