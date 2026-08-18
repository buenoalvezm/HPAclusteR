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

#' Visualize cluster comparison results with a bubble heatmap
#'
#' @param summary_data A tibble containing cluster comparison results
#'
#' @returns A buble heatmap plot
#' @keywords internal
visualize_comparison_heatmap <- function(summary_data) {
  # Remove rows where adj_p_val is above 0.05 before clustering
  summary_data <- summary_data |>
    dplyr::filter(!!rlang::sym("adj_p_val") < 0.05)

  heatmap_matrix <- summary_data |>
    dplyr::select(dplyr::any_of(c("cluster_A", "cluster_B", "n_overlap"))) |>
    tidyr::pivot_wider(
      names_from = !!rlang::sym("cluster_B"),
      values_from = !!rlang::sym("n_overlap"),
      values_fill = 0
    ) |>
    tibble::column_to_rownames("cluster_A") |>
    as.matrix()

  # Perform hierarchical clustering on rows and columns
  row_dendrogram <- stats::hclust(stats::dist(heatmap_matrix))
  col_dendrogram <- stats::hclust(stats::dist(t(heatmap_matrix)))

  row_order <- rownames(heatmap_matrix)[row_dendrogram[["order"]]]
  col_order <- colnames(heatmap_matrix)[col_dendrogram[["order"]]]

  summary_data <- summary_data |>
    dplyr::mutate(
      cluster_A = factor(!!rlang::sym("cluster_A"), levels = row_order),
      cluster_B = factor(!!rlang::sym("cluster_B"), levels = col_order),
      color_value = -log10(!!rlang::sym("adj_p_val"))
    )

  heatmap_plot <- ggplot2::ggplot(
    summary_data,
    ggplot2::aes(x = !!rlang::sym("cluster_B"), y = !!rlang::sym("cluster_A"))
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
      x = "Cluster B",
      y = "Cluster A"
    ) +
    theme_hc() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        colour = "lightgrey",
        linewidth = 0.2
      )
    )

  return(heatmap_plot)
}

#' Visualize cluster comparison results with a network graph
#'
#' @param matched_clusters A tibble containing matched clusters between two datasets
#'
#' @returns A network graph plot
#' @keywords internal
visualize_comparison_net <- function(matched_clusters) {
  check_installed("igraph", "to lay out the cluster comparison network")

  if (nrow(matched_clusters) == 0L) {
    stop(
      "No significant cluster matches to draw a network from.",
      call. = FALSE
    )
  }

  # Prefix node names so that clusters from the two datasets stay distinct.
  edges <- matched_clusters |>
    dplyr::mutate(
      from = paste0("A_", !!rlang::sym("cluster_A")),
      to = paste0("B_", !!rlang::sym("cluster_B"))
    ) |>
    dplyr::select(dplyr::all_of(c("from", "to", "percentage_overlap")))

  nodes <- dplyr::bind_rows(
    matched_clusters |>
      dplyr::transmute(
        name = paste0("A_", !!rlang::sym("cluster_A")),
        label = as.character(!!rlang::sym("cluster_A")),
        size = !!rlang::sym("n_genes_A"),
        dataset = "Dataset A"
      ),
    matched_clusters |>
      dplyr::transmute(
        name = paste0("B_", !!rlang::sym("cluster_B")),
        label = as.character(!!rlang::sym("cluster_B")),
        size = !!rlang::sym("n_genes_B"),
        dataset = "Dataset B"
      )
  ) |>
    dplyr::distinct()

  graph <- igraph::graph_from_data_frame(
    edges,
    directed = FALSE,
    vertices = nodes
  )

  set.seed(42)
  layout <- igraph::layout_with_fr(graph)
  colnames(layout) <- c("x", "y")

  node_positions <- dplyr::bind_cols(
    tibble::tibble(name = igraph::vertex_attr(graph, "name")),
    tibble::as_tibble(layout)
  ) |>
    dplyr::left_join(nodes, by = "name")

  edge_positions <- edges |>
    dplyr::left_join(
      node_positions |>
        dplyr::select(from = "name", x = "x", y = "y"),
      by = "from"
    ) |>
    dplyr::left_join(
      node_positions |>
        dplyr::select(to = "name", xend = "x", yend = "y"),
      by = "to"
    )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_positions,
      ggplot2::aes(
        x = !!rlang::sym("x"),
        y = !!rlang::sym("y"),
        xend = !!rlang::sym("xend"),
        yend = !!rlang::sym("yend"),
        alpha = !!rlang::sym("percentage_overlap")
      ),
      colour = "grey40",
      linewidth = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::scale_alpha_continuous(range = c(0.15, 0.9), limits = c(0, 100)) +
    ggplot2::geom_point(
      data = node_positions,
      ggplot2::aes(
        x = !!rlang::sym("x"),
        y = !!rlang::sym("y"),
        size = !!rlang::sym("size"),
        fill = !!rlang::sym("dataset"),
        shape = !!rlang::sym("dataset")
      ),
      colour = "black",
      stroke = 0.3
    ) +
    ggplot2::scale_shape_manual(
      values = c("Dataset A" = 21, "Dataset B" = 24)
    ) +
    ggplot2::scale_fill_manual(
      values = c("Dataset A" = "#66C2A5", "Dataset B" = "#FC8D62")
    ) +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggrepel::geom_text_repel(
      data = node_positions,
      ggplot2::aes(
        x = !!rlang::sym("x"),
        y = !!rlang::sym("y"),
        label = !!rlang::sym("label")
      ),
      size = 3,
      seed = 42,
      max.overlaps = Inf
    ) +
    ggplot2::guides(size = "none", shape = "none", fill = "none") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      caption = paste0(
        "Node colour/shape = Dataset, Node size = Number of genes in each cluster\n",
        "Edge opacity = Percentage overlap between two clusters"
      )
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        face = "italic"
      )
    )
}
