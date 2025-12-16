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
      panel.grid.major = ggplot2::element_line(color = "lightgrey", size = 0.2)
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
  set.seed(42)

  # Add prefixes to ensure unique node names
  edges <- matched_clusters |>
    dplyr::mutate(
      from = paste0("A_", !!rlang::sym("cluster_A")),
      to = paste0("B_", !!rlang::sym("cluster_B"))
    ) |>
    dplyr::select(dplyr::any_of(c("from", "to", "percentage_overlap")))

  # Combine node sizes from both datasets
  node_sizes <- matched_clusters |>
    dplyr::select(dplyr::any_of(c(
      "cluster_A",
      "cluster_B",
      "n_genes_A",
      "n_genes_B"
    ))) |>
    tidyr::pivot_longer(
      cols = c(!!rlang::sym("cluster_A"), !!rlang::sym("cluster_B")),
      names_to = "dataset",
      values_to = "cluster"
    ) |>
    dplyr::mutate(
      size = ifelse(
        !!rlang::sym("dataset") == "cluster_A",
        !!rlang::sym("n_genes_A"),
        !!rlang::sym("n_genes_B")
      ),
      cluster = paste0(
        ifelse(!!rlang::sym("dataset") == "cluster_A", "A_", "B_"),
        !!rlang::sym("cluster")
      ),
      shape = ifelse(
        !!rlang::sym("dataset") == "cluster_A",
        "circle",
        "triangle"
      ) # Assign shapes
    ) |>
    dplyr::distinct(
      !!rlang::sym("cluster"),
      !!rlang::sym("size"),
      !!rlang::sym("shape")
    )

  net <- network::network(edges, directed = FALSE)

  # Add node attributes
  network::set.vertex.attribute(
    net,
    "dataset",
    ifelse(
      network::get.vertex.attribute(net, "vertex.names") %in% edges$from,
      "Dataset A",
      "Dataset B"
    )
  )
  network::set.vertex.attribute(
    net,
    "size",
    node_sizes[["size"]][match(
      network::get.vertex.attribute(net, "vertex.names"),
      node_sizes[["cluster"]]
    )]
  )
  network::set.vertex.attribute(
    net,
    "shape",
    node_sizes[["shape"]][match(
      network::get.vertex.attribute(net, "vertex.names"),
      node_sizes[["cluster"]]
    )]
  )

  # Add edge attributes
  network::set.edge.attribute(
    net,
    "color",
    ifelse(
      edges[["percentage_overlap"]] > 90,
      "#9f9f9f", # Fully solid grey
      ifelse(
        edges[["percentage_overlap"]] > 80,
        "#a8a8a8", # 90% solid
        ifelse(
          edges[["percentage_overlap"]] > 70,
          "#b2b2b2", # 80% solid
          ifelse(
            edges[["percentage_overlap"]] > 60,
            "#bbbbbb", # 70% solid
            ifelse(
              edges[["percentage_overlap"]] > 50,
              "#c4c4c4", # 60% solid
              ifelse(
                edges[["percentage_overlap"]] > 40,
                "#cecece", # 50% solid
                ifelse(
                  edges[["percentage_overlap"]] > 30,
                  "#d8d8d8", # 40% solid
                  ifelse(
                    edges[["percentage_overlap"]] > 20,
                    "#e1e1e1", # 30% solid
                    "#ebebeb" # 20% solid
                  )
                )
              )
            )
          )
        )
      )
    )
  )

  if (!requireNamespace("sna", quietly = TRUE)) {
    message("Package 'sna' not found. Installing now...")
    install.packages("sna")
  }

  network <- GGally::ggnet2(
    net,
    color = "dataset",
    palette = "Set2",
    size = "size",
    shape = "shape",
    edge.color = "color",
    edge.size = 0.5,
    label = TRUE,
    label.size = 3,
    vjust = 2
  ) +
    ggplot2::guides(size = "none", shape = "none", color = "none") +
    ggplot2::labs(
      caption = "Node color/shape = Dataset, Node size = Number of genes in each cluster\nEdge color = Percentage overlap between two clusters"
    ) +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(
        hjust = 0.5,
        size = 10,
        face = "italic"
      )
    )

  return(network)
}
