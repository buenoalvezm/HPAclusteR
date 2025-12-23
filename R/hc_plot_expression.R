#' Plot per-cluster heatmaps for scaled and z-score expression
#'
#' `hc_plot_expression()` creates heatmaps of gene expression for each cluster identified in the AnnDatR object. It generates two types of plots per cluster: one showing z-score normalized expression and another showing scaled relative expression. Each heatmap is accompanied by a confidence bar indicating the membership confidence of genes in the respective cluster.
#'
#' @param AnnDatR AnnDatR object (genes x samples)
#' @param value_limits List with zscore and scaled limits (default: zscore = c(-2,2), scaled = c(0,1))
#' @param show_sample_labels Logical, whether to show sample labels on x-axis (default: TRUE)
#'
#' @return List with two named lists: $zscore and $scaled, each a list of ggplot objects per cluster
#' @export
#'
#' @examples
#' # Run clustering pipeline
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 7)
#'
#' # Plot expression heatmaps
#' expression_plots <- hc_plot_expression(adata_res, show_sample_labels = FALSE)
#' expression_plots$zscore[["1"]]  # Z-score plot for cluster 1
#' expression_plots$scaled[["1"]]  # Scaled plot for cluster 1
hc_plot_expression <- function(
  AnnDatR,
  value_limits = list(zscore = c(-2, 2), scaled = c(0, 1)),
  show_sample_labels = TRUE
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "The 'patchwork' package is required for this function. Please install it using install.packages('patchwork')."
    )
  }
  if (is.null(AnnDatR[["obs"]][["cluster"]])) {
    stop(
      "AnnDatR$obs$cluster not found. Call `hc_cluster_consensus()` before `hc_plot_expression()`."
    )
  }
  # Prepare long expression data
  sample_col <- AnnDatR[["var_names_col"]]
  gene_col <- AnnDatR[["obs_names_col"]]
  X <- AnnDatR[["X"]]
  obs <- AnnDatR[["obs"]]
  var <- AnnDatR[["var"]]

  expr_long <- X |>
    tidyr::pivot_longer(
      -!!rlang::sym("sample_col"),
      names_to = "gene",
      values_to = "value"
    ) |>
    dplyr::rename(sample = !!rlang::sym("sample_col"))

  # Add cluster info and filter NAs
  expr_long <- expr_long |>
    dplyr::left_join(
      obs |>
        dplyr::select(
          gene = !!rlang::sym("gene_col"),
          cluster = "cluster"
        ),
      by = "gene"
    ) |>
    dplyr::filter(!is.na(!!rlang::sym("cluster")))

  # Calculate scaled and z-score per gene
  expr_long <- expr_long |>
    dplyr::group_by(!!rlang::sym("gene")) |>
    dplyr::mutate(
      value_scaled = !!rlang::sym("value") /
        max(!!rlang::sym("value"), na.rm = TRUE),
      value_zscore = (!!rlang::sym("value") -
        mean(!!rlang::sym("value"), na.rm = TRUE)) /
        stats::sd(!!rlang::sym("value"), na.rm = TRUE)
    ) |>
    dplyr::ungroup()

  clusters <- sort(unique(expr_long[["cluster"]]))

  # Membership confidence
  memb_conf <- AnnDatR[["uns"]][["membership_matrix"]] |>
    dplyr::rename(
      cluster_membership_confidence = !!rlang::sym("membership"),
      gene = !!rlang::sym("gene")
    )

  # Helper for ordering
  get_order <- function(mat) {
    if (nrow(mat) < 2 || ncol(mat) < 2) {
      return(list(row = rownames(mat), col = colnames(mat)))
    }
    rowclust <- stats::hclust(stats::dist(mat))
    colclust <- stats::hclust(stats::dist(t(mat)))
    list(
      row = rownames(mat)[rowclust[["order"]]],
      col = colnames(mat)[colclust[["order"]]]
    )
  }

  axis_x <- if (show_sample_labels) {
    ggplot2::element_text(
      angle = -90,
      hjust = 0,
      vjust = 0.5
    )
  } else {
    ggplot2::element_blank()
  }

  # Z-score plots
  plots <- lapply(clusters, function(clust) {
    dat <- expr_long |> dplyr::filter(!!rlang::sym("cluster") == clust)
    mat <- dat |>
      dplyr::select(dplyr::any_of(c("gene", "sample", "value_zscore"))) |>
      tidyr::pivot_wider(
        names_from = !!rlang::sym("sample"),
        values_from = !!rlang::sym("value_zscore")
      ) |>
      tibble::column_to_rownames("gene") |>
      as.matrix()

    ord <- get_order(mat)
    dat <- dat |>
      dplyr::mutate(
        gene = factor(!!rlang::sym("gene"), rev(ord[["row"]])),
        sample = factor(!!rlang::sym("sample"), ord[["col"]])
      )
    conf_data <- memb_conf |>
      dplyr::filter(!!rlang::sym("cluster") == clust) |>
      dplyr::filter(!!rlang::sym("gene") %in% levels(dat[["gene"]])) |>
      dplyr::mutate(gene = factor(!!rlang::sym("gene"), levels(dat[["gene"]])))

    plot_conf <- ggplot2::ggplot(
      conf_data,
      ggplot2::aes(
        x = "Confidence",
        y = !!rlang::sym("gene"),
        fill = !!rlang::sym("cluster_membership_confidence")
      )
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(
        low = "#D1EEEA",
        high = "#2A5674",
        limits = c(0.2, 1),
        breaks = c(0.5, 0.75, 1),
        name = "Confidence"
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = -90, face = "bold"),
        legend.position = "left"
      )

    plot_heat_z <- ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x = !!rlang::sym("sample"),
        y = !!rlang::sym("gene"),
        fill = !!rlang::sym("value_zscore")
      )
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "#4575B4",
        high = "#D73027",
        mid = "#FFFFFF",
        name = "Z-score",
        midpoint = 0
      ) +
      theme_hc() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.x = axis_x,
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        title = ggplot2::element_text(face = "bold"),
        legend.title = ggplot2::element_text(face = "plain"),
        legend.position = "top"
      ) +
      ggplot2::ggtitle(paste("Cluster", clust, "- Z-score"))

    plot_heat_scaled <- ggplot2::ggplot(
      dat,
      ggplot2::aes(
        x = !!rlang::sym("sample"),
        y = !!rlang::sym("gene"),
        fill = !!rlang::sym("value_scaled")
      )
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "#4575B4",
        high = "#D73027",
        mid = "#FFFFFF",
        name = "Relative Expression",
        limits = value_limits[["scaled"]],
        breaks = c(0, 0.5, 1)
      ) +
      theme_hc() +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text.x = axis_x,
        axis.text.y = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        title = ggplot2::element_text(face = "bold"),
        legend.title = ggplot2::element_text(face = "plain"),
        legend.position = "top"
      ) +
      ggplot2::ggtitle(paste("Cluster", clust, "- Scaled"))

    list(
      zscore = plot_conf +
        plot_heat_z +
        patchwork::plot_layout(widths = c(0.05, 1)),
      scaled = plot_conf +
        plot_heat_scaled +
        patchwork::plot_layout(widths = c(0.05, 1))
    )
  })
  names(plots) <- as.character(clusters)

  list(
    zscore = lapply(plots, `[[`, "zscore"),
    scaled = lapply(plots, `[[`, "scaled")
  )
}
