calculate_tau_score_ <- function(X, gene_col) {
  # Check if X is a tibble
  if (!is_tibble(X)) {
    stop("Input X must be a tibble.")
  }

  # Check if the gene_col exists
  if (!gene_col %in% names(X)) {
    stop(paste("Gene identifier column '", gene_col, "' not found in X."))
  }

  # Pivot the data to long format
  long_data <- X %>%
    pivot_longer(
      cols = -{{ gene_col }},
      names_to = "tissue",
      values_to = "expression"
    ) %>%
    rename(gene = {{ gene_col }})

  # Calculate Tau score per gene
  tau_scores <- long_data %>%
    group_by(gene) %>%
    summarise(
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
