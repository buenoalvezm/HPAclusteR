#' Calculate Adjusted Rand Index (ARI) between two clusterings
#'
#' @param labels1 First set of cluster labels
#' @param labels2 Second set of cluster labels
#'
#' @returns Adjusted Rand Index (ARI) value
#' @keywords internal
calculate_ari <- function(labels1, labels2) {
  # The Adjusted Rand Index (ARI) algorithm evaluates the similarity between two clusterings, adjusted for chance.
  # Steps:
  # 1. Construct a contingency table from the two cluster assignments.
  # 2. Compute the total number of pairs of elements in the dataset.
  # 3. Calculate the number of pairs within clusters for both clusterings and their overlap.
  # 4. Compute the expected index (chance agreement) and the maximum index (perfect agreement).
  # 5. Use the ARI formula: (observed - expected) / (maximum - expected) to adjust for chance.
  contingency <- table(labels1, labels2)
  sum_rows <- rowSums(contingency)
  sum_cols <- colSums(contingency)
  total_pairs <- choose(sum(contingency), 2)
  sum_comb_c <- sum(choose(contingency, 2))
  sum_comb_r <- sum(choose(sum_rows, 2))
  sum_comb_s <- sum(choose(sum_cols, 2))
  expected_index <- (sum_comb_r * sum_comb_s) / total_pairs
  max_index <- (sum_comb_r + sum_comb_s) / 2
  ari <- (sum_comb_c - expected_index) / (max_index - expected_index)
  return(ari)
}

#' Visualize ARI scores with histogram
#'
#' @param ARI_scores A tibble containing ARI scores
#'
#' @returns A list containing the histogram plot and ARI scores
#' @keywords internal
visualize_ari <- function(ARI_scores) {
  ari_ranges <- tibble::tibble(
    xmin = c(-1, 0, 0.4, 0.7),
    xmax = c(0, 0.4, 0.7, 1),
    fill = factor(
      c("red", "yellow", "lightgreen", "green"),
      levels = c("red", "yellow", "lightgreen", "green")
    )
  )

  ari_hist <- ggplot2::ggplot(
    ARI_scores,
    ggplot2::aes(x = !!rlang::sym("ARI"))
  ) +
    ggplot2::geom_rect(
      data = ari_ranges,
      ggplot2::aes(
        xmin = !!rlang::sym("xmin"),
        xmax = !!rlang::sym("xmax"),
        ymin = -Inf,
        ymax = Inf,
        fill = !!rlang::sym("fill")
      ),
      inherit.aes = FALSE,
      alpha = 0.2
    ) +
    ggplot2::geom_histogram(
      binwidth = 0.05,
      color = "black",
      fill = "white"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "red" = "red",
        "yellow" = "yellow",
        "lightgreen" = "lightgreen",
        "green" = "green"
      ),
      guide = "legend",
      labels = c(
        "Poor (< 0)",
        "Weak (0-0.5)",
        "Moderate (0.5-0.8)",
        "Strong (> 0.8)"
      )
    ) +
    ggplot2::labs(
      x = "Adjusted Rand Index (ARI)",
      y = "Count",
      fill = "Stability Range"
    ) +
    theme_hc()

  return(list(cluster_stability_hist = ari_hist, ari_scores = ARI_scores))
}

#' Calculate ARI scores for multiple clustering seeds
#'
#' `hc_cluster_stability()` computes the Adjusted Rand Index (ARI) scores between clustering results obtained from multiple random seeds.
#'
#' @param AnnDatR AnnDatR An AnnDatR object containing clustering results.
#'
#' @returns A tibble containing pairs of seeds and their corresponding ARI scores.
#'
#' @export
#' @examples
#' # Perform clustering with multiple seeds and calculate ARI scores
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_consensus_cluster(adata_res, resolution = 6.3)
#' hc_cluster_stability(adata_res)
hc_cluster_stability <- function(AnnDatR) {
  if (is.null(AnnDatR[["uns"]][["cluster_data"]])) {
    stop(
      "AnnDatR$uns$cluster_data not found. Call `hc_consensus_cluster()` before `hc_cluster_stability()`."
    )
  }

  n_seeds <- ncol(AnnDatR[["uns"]][["cluster_data"]]) - 2

  ARI_scores <- tidyr::expand_grid(
    seed1 = 1:n_seeds,
    seed2 = 1:n_seeds
  ) |>
    dplyr::filter(!!rlang::sym("seed1") < !!rlang::sym("seed2")) |>
    dplyr::group_by(!!rlang::sym("seed1"), !!rlang::sym("seed2")) |>
    dplyr::mutate(
      ARI = calculate_ari(
        AnnDatR[["uns"]][["cluster_data"]] |>
          dplyr::select(
            !!rlang::sym("gene"),
            paste("seed", !!rlang::sym("seed1"), sep = "_")
          ) |>
          tibble::deframe(),
        AnnDatR[["uns"]][["cluster_data"]] |>
          dplyr::select(
            !!rlang::sym("gene"),
            paste("seed", !!rlang::sym("seed2"), sep = "_")
          ) |>
          tibble::deframe()
      )
    )

  ARI_res <- visualize_ari(ARI_scores)

  return(ARI_res)
}
