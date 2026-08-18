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
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 8, n_seeds = 20)
#' hc_cluster_stability(adata_res)
hc_cluster_stability <- function(AnnDatR) {
  cluster_data <- AnnDatR[["uns"]][["cluster_data"]]
  if (is.null(cluster_data)) {
    stop(
      paste0(
        "AnnDatR$uns$cluster_data not found. ",
        "Call `hc_cluster_consensus()` before `hc_cluster_stability()`."
      ),
      call. = FALSE
    )
  }

  # Pull the per-seed labels out once, rather than re-selecting and deframing
  # the tibble inside every one of the n_seeds * (n_seeds - 1) / 2 comparisons.
  seed_columns <- grep("^seed_", names(cluster_data), value = TRUE)
  seed_numbers <- as.integer(sub("^seed_", "", seed_columns))
  seed_columns <- seed_columns[order(seed_numbers)]
  seed_numbers <- sort(seed_numbers)

  labels <- as.matrix(cluster_data[, seed_columns, drop = FALSE])
  n_seeds <- length(seed_columns)

  if (n_seeds < 2L) {
    stop(
      "At least two clustering seeds are required to assess stability.",
      call. = FALSE
    )
  }

  pairs <- utils::combn(n_seeds, 2L)

  ARI_scores <- tibble::tibble(
    seed1 = seed_numbers[pairs[1, ]],
    seed2 = seed_numbers[pairs[2, ]],
    ARI = vapply(
      seq_len(ncol(pairs)),
      function(pair) {
        calculate_ari(labels[, pairs[1, pair]], labels[, pairs[2, pair]])
      },
      numeric(1)
    )
  )

  visualize_ari(ARI_scores)
}
