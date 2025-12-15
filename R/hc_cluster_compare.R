#' Visualize cluster comparison results with a bubble heatmap
#'
#' @param summary_data A tibble containing cluster comparison results
#'
#' @returns A buble heatmap plot
#' @keywords internal
visualize_comparison <- function(summary_data) {
  # Remove rows where adj_p_val is above 0.05 before clustering
  summary_data <- summary_data |>
    dplyr::filter(!!rlang::sym("adj_p_val") < 0.05)

  heatmap_matrix <- summary_data |>
    dplyr::select(dplyr::any_of(c("cluster_1", "cluster_2", "n_overlap"))) |>
    tidyr::pivot_wider(
      names_from = !!rlang::sym("cluster_2"),
      values_from = !!rlang::sym("n_overlap"),
      values_fill = 0
    ) |>
    tibble::column_to_rownames("cluster_1") |>
    as.matrix()

  # Perform hierarchical clustering on rows and columns
  row_dendrogram <- stats::hclust(stats::dist(heatmap_matrix))
  col_dendrogram <- stats::hclust(stats::dist(t(heatmap_matrix)))

  # Reorder rows and columns based on clustering
  row_order <- rownames(heatmap_matrix)[row_dendrogram[["order"]]]
  col_order <- colnames(heatmap_matrix)[col_dendrogram[["order"]]]

  summary_data <- summary_data |>
    dplyr::mutate(
      cluster_1 = factor(!!rlang::sym("cluster_1"), levels = row_order),
      cluster_2 = factor(!!rlang::sym("cluster_2"), levels = col_order),
      color_value = -log10(!!rlang::sym("adj_p_val")) # Use -log10(adj_p_val) for coloring
    )

  # Plot heatmap with ggplot2
  heatmap_plot <- ggplot2::ggplot(
    summary_data,
    ggplot2::aes(x = !!rlang::sym("cluster_2"), y = !!rlang::sym("cluster_1"))
  ) +
    ggplot2::geom_point(
      ggplot2::aes(
        color = !!rlang::sym("color_value"),
        size = !!rlang::sym("n_overlap")
      ),
    ) +
    ggplot2::scale_color_gradient(
      na.value = NA, # Invisible for adj_p_val >= 0.05
      name = "-log10(adj_p_val)"
    ) +
    ggplot2::scale_size(range = c(1, 8), name = "Overlap") +
    ggplot2::labs(
      x = "Cluster 2",
      y = "Cluster 1"
    ) +
    theme_hc() +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(), # Remove axis lines
      axis.ticks = ggplot2::element_blank(), # Remove axis ticks
      panel.grid.major = ggplot2::element_line(color = "lightgrey", size = 0.2) # Light grey grid
    )

  return(heatmap_plot)
}

#' Compare clusters using hypergeometric test
#'
#' `hc_cluster_compare()` compares clusters from two AnnDatR objects using a hypergeometric test to assess the significance of overlap between clusters.
#'
#' @param AnnDatR1 First AnnDatR object
#' @param AnnDatR2 Second AnnDatR object
#'
#' @returns A tibble summarizing the pairwise cluster comparisons, including overlap, percentage overlap, p-values, and adjusted p-values.
#' @export
#'
#' @examples
#' # Run clustering 1
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res_1 <- hc_cluster_consensus(adata_res, resolution = 6)
#'
#' # Run clustering 2
#' adata_res_2 <- hc_cluster_consensus(adata_res, resolution = 7)
#'
#' # Compare clusters between the two results
#' hc_cluster_compare(adata_res_1, adata_res_2)
hc_cluster_compare <- function(AnnDatR1, AnnDatR2) {
  # Check if consensus clustering exists in both objects
  if (is.null(AnnDatR1[["uns"]][["consensus_clustering"]])) {
    stop(
      "Consensus clustering not found in the first AnnDatR object. Run `hc_cluster_consensus()` first."
    )
  }
  if (is.null(AnnDatR2[["uns"]][["consensus_clustering"]])) {
    stop(
      "Consensus clustering not found in the second AnnDatR object. Run `hc_cluster_consensus()` first."
    )
  }

  # Extract consensus clustering
  consensus1 <- AnnDatR1[["uns"]][["consensus_clustering"]]
  consensus2 <- AnnDatR2[["uns"]][["consensus_clustering"]]

  # Align genes between the two datasets
  common_genes <- dplyr::inner_join(
    consensus1,
    consensus2,
    by = "gene",
    suffix = c("_1", "_2")
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
    dplyr::group_by(!!rlang::sym("cluster_1"), !!rlang::sym("cluster_2")) |>
    dplyr::summarize(
      n_genes_1 = dplyr::n_distinct(consensus1[["gene"]][
        consensus1[["cluster"]] == unique(!!rlang::sym("cluster_1"))
      ]),
      n_genes_2 = dplyr::n_distinct(consensus2[["gene"]][
        consensus2[["cluster"]] == unique(!!rlang::sym("cluster_2"))
      ]),
      n_overlap = dplyr::n(),
      percentage_overlap = round(
        100 *
          !!rlang::sym("n_overlap") /
            min(!!rlang::sym("n_genes_1"), !!rlang::sym("n_genes_2")),
        2
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      # Perform one-sided hypergeometric test
      p_val = stats::phyper(
        !!rlang::sym("n_overlap") - 1, # Overlap minus 1 (one-sided test)
        !!rlang::sym("n_genes_1"), # Total genes in cluster_1
        !!rlang::sym("total_genes_1") - !!rlang::sym("n_genes_1"), # Genes not in cluster_1
        !!rlang::sym("n_genes_2"), # Total genes in cluster_2
        lower.tail = FALSE # One-sided test
      ),
      adj_p_val = stats::p.adjust(!!rlang::sym("p_val"), method = "BH")
    )

  heatmap <- visualize_comparison(results)
}
