#' Calculate Tau Score
#'
#' Computes the tissue specificity measure Tau (Yanai et al., 2005) for each gene
#' based on expression values across samples (e.g. tissues, single cell types, etc.). Tau ranges from 0 (ubiquitous expression)
#' to 1 (tissue-specific expression).
#'
#' @param df A tibble containing gene expression values (genes as rows, samples (e.g. tissues, single cell types, etc.) as columns).
#' @param gene_col Name of the column in `df` that contains gene identifiers.
#'
#' @return A tibble with two columns:
#'   \describe{
#'     \item{gene}{Gene identifier.}
#'     \item{tau_score}{Tau score for each gene.}
#'   }
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tidyr)
#' example_data <- tibble::tibble(
#'   gene_id = c("gene1", "gene2"),
#'   tissueA = c(10, 2),
#'   tissueB = c(8, 0),
#'   tissueC = c(7, 0)
#' )
#' calculate_tau_score(example_data, gene_col = "gene_id")

calculate_tau_score <- function(df, gene_col) {
  # Check if X is a tibble
  if (!tibble::is_tibble(df)) {
    stop("Input must be a tibble.")
  }

  # Check if the gene_col exists
  if (!gene_col %in% names(df)) {
    stop(paste(
      "Gene identifier column '",
      gene_col,
      "' not found in tibble."
    ))
  }

  # Pivot the data to long format
  long_data <- df %>%
    tidyr::pivot_longer(
      cols = -{{ gene_col }},
      names_to = "tissue",
      values_to = "expression"
    ) %>%
    dplyr::rename(gene = {{ gene_col }})

  # Calculate Tau score per gene
  tau_scores <- long_data %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      tau_score = {
        expression_finite <- expression[!is.na(expression)]
        n_finite <- length(expression_finite)

        if (n_finite <= 1) {
          NA_real_
        } else {
          max_val <- max(expression_finite, na.rm = TRUE)
          if (max_val == 0) {
            NA_real_
          } else {
            sum_val <- sum(1 - (expression_finite / max_val), na.rm = TRUE)
            sum_val / (n_finite - 1)
          }
        }
      },
      .groups = "drop"
    )

  return(tau_scores)
}
