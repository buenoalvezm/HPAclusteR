#' Cluster matching algorithm
#'
#' `cluster_matching()` implements an algorithm to match clusters from two datasets based on significant overlaps, ensuring that there are no many-to-many matches, while still possible to have one-to-many or many-to-one matches.
#'
#' @param cluster_comparison A tibble containing pairwise cluster comparison results, including columns: cluster_A, cluster_B, n_overlap, adj_p_val.
#' @param total_genes Total number of genes considered in the comparison.
#'
#' @returns A list containing:
#' - matched_clusters: A tibble of matched clusters adhering to the matching rules.
#' - match_score: A numeric score representing the total normalized overlap percentage.
#' @keywords internal
cluster_matching <- function(
  cluster_comparison,
  total_genes
) {
  # Step 1: Filter significant matches based on threshold
  matches <- cluster_comparison |>
    dplyr::filter(!!rlang::sym("adj_p_val") < 0.05) |>
    dplyr::select(dplyr::any_of(c("cluster_A", "cluster_B", "n_overlap")))

  # Step 2: Sort matches by overlap (descending) to prioritize stronger matches
  matches <- matches |>
    dplyr::arrange(dplyr::desc(!!rlang::sym("n_overlap")))

  # Step 3: Initialize empty results and track matched clusters
  matched_clusters <- dplyr::tibble(
    cluster_A = character(),
    cluster_B = character(),
    n_overlap = numeric()
  )
  matched_1 <- character() # Track matched clusters in cluster_A
  matched_2 <- character() # Track matched clusters in cluster_B

  # Step 4: Iterate through matches and add valid ones
  for (i in seq_len(nrow(matches))) {
    row <- matches[i, ]
    # Check if the match violates the many-to-many rule
    if (
      !(row[["cluster_A"]] %in% matched_1 && row[["cluster_B"]] %in% matched_2)
    ) {
      matched_clusters <- dplyr::bind_rows(matched_clusters, row)
      matched_1 <- c(matched_1, row[["cluster_A"]])
      matched_2 <- c(matched_2, row[["cluster_B"]])
    }
  }

  # Step 5: Calculate total weight (normalized similarity)
  total_weight <- sum(matched_clusters[["n_overlap"]]) / total_genes * 100

  return(list(
    matched_clusters = matched_clusters,
    match_score = round(total_weight, 2)
  ))
}

#' Compare clusters using hypergeometric test
#'
#' `hc_cluster_compare()` compares clusters from two AnnDatR objects using a hypergeometric test to assess the significance of overlap between clusters.
#'
#' @param AnnDatR_A First AnnDatR object
#' @param AnnDatR_B Second AnnDatR object
#' @param graph_type Type of graph to return: "bipartite" for simple overlap graph, "optimized" for optimized matching graph (default: "bipartite")
#'
#' @details
#' If `graph_type` is set to "optimized", the function will run a cluster matching algorithm to ensure no many-to-many matches between clusters, while still allowing one-to-many or many-to-one matches.
#' The algorithm will:
#' 1) Filter significant matches based on adjusted p-value threshold (0.05)
#' 2) Sort matches by overlap size in descending order
#' 3) Iteratively add matches to the result set, ensuring no many-to-many matches occur.
#' 4) Calculate a match score representing the total normalized overlap percentage.
#'
#' @return A list containing:
#' - matches: A tibble of cluster comparisons with overlap statistics and p-values.
#' - heatmap: A ggplot2 heatmap visualizing cluster overlaps.
#' - match_score (if `graph_type` is "optimized"): A numeric score representing the total normalized overlap percentage.
#' @export
#'
#' @examples
#' # Run clustering A
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res_a <- hc_cluster_consensus(adata_res, resolution = 6)
#'
#' # Run clustering B
#' adata_res_b <- hc_cluster_consensus(adata_res, resolution = 7)
#'
#' # Compare clusters between the two results
#' hc_cluster_compare(adata_res_a, adata_res_b)
hc_cluster_compare <- function(
  AnnDatR_A,
  AnnDatR_B,
  graph_type = c("bipartite", "optimized")
) {
  graph_type <- match.arg(graph_type)

  # Check if consensus clustering exists in both objects
  if (is.null(AnnDatR_A[["uns"]][["consensus_clustering"]])) {
    stop(
      "Consensus clustering not found in the first AnnDatR object. Run `hc_cluster_consensus()` first."
    )
  }
  if (is.null(AnnDatR_B[["uns"]][["consensus_clustering"]])) {
    stop(
      "Consensus clustering not found in the second AnnDatR object. Run `hc_cluster_consensus()` first."
    )
  }

  # Extract consensus clustering
  consensus1 <- AnnDatR_A[["uns"]][["consensus_clustering"]]
  consensus2 <- AnnDatR_B[["uns"]][["consensus_clustering"]]

  # Align genes between the two datasets
  common_genes <- dplyr::inner_join(
    consensus1,
    consensus2,
    by = "gene",
    suffix = c("_A", "_B")
  )

  # Calculate the percentage of genes that are not in both datasets
  total_genes_1 <- nrow(consensus1)
  total_genes_2 <- nrow(consensus2)
  common_gene_count <- nrow(common_genes)
  unmatched_percentage <- 100 *
    (1 - (common_gene_count / min(total_genes_1, total_genes_2)))

  # Notify the user about unmatched genes
  if (unmatched_percentage > 10) {
    warning(sprintf(
      "More than 10%% of genes do not match between the two datasets (%.2f%% unmatched).",
      unmatched_percentage
    ))
  } else {
    message(sprintf(
      "%.2f%% of genes do not match between the two datasets.",
      unmatched_percentage
    ))
  }

  # Perform pairwise cluster comparisons
  results <- common_genes |>
    dplyr::group_by(!!rlang::sym("cluster_A"), !!rlang::sym("cluster_B")) |>
    dplyr::summarize(
      n_genes_A = dplyr::n_distinct(consensus1[["gene"]][
        consensus1[["cluster"]] == unique(!!rlang::sym("cluster_A"))
      ]),
      n_genes_B = dplyr::n_distinct(consensus2[["gene"]][
        consensus2[["cluster"]] == unique(!!rlang::sym("cluster_B"))
      ]),
      n_overlap = dplyr::n(),
      percentage_overlap = round(
        100 *
          !!rlang::sym("n_overlap") /
            min(!!rlang::sym("n_genes_A"), !!rlang::sym("n_genes_B")),
        2
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      # Perform one-sided hypergeometric test
      p_val = stats::phyper(
        !!rlang::sym("n_overlap") - 1, # Overlap minus 1 (one-sided test)
        !!rlang::sym("n_genes_A"), # Total genes in cluster_A
        !!rlang::sym("total_genes_1") - !!rlang::sym("n_genes_A"), # Genes not in cluster_A
        !!rlang::sym("n_genes_B"), # Total genes in cluster_B
        lower.tail = FALSE # One-sided test
      ),
      adj_p_val = stats::p.adjust(!!rlang::sym("p_val"), method = "BH")
    )

  heatmap <- visualize_comparison_heatmap(results)

  if (graph_type == "optimized") {
    matches_list <- cluster_matching(results, common_gene_count)
    matches <- dplyr::inner_join(
      results,
      matches_list[["matched_clusters"]] |>
        dplyr::select(-dplyr::any_of(c("n_overlap"))),
      by = c("cluster_A", "cluster_B")
    )
    network <- visualize_comparison_net(matches)
    return(list(
      matches = matches,
      heatmap = heatmap,
      network = network,
      match_score = matches_list[["match_score"]]
    ))
  } else {
    matches <- results |>
      dplyr::filter(!!rlang::sym("adj_p_val") < 0.05)
    network <- visualize_comparison_net(matches)
    return(list(
      matches = matches,
      heatmap = heatmap,
      network = network
    ))
  }
}
