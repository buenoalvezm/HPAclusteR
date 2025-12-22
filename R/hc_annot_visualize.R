#' Plot treemap for each cluster in reduced GO terms
#'
#' @param reduced_terms Data frame, e.g. simp_res$reducedTerms
#' @param ontology_pal Named vector of ontology colors (default: c(MF="#FF6F00", CC="#C71B00", BP="#018EA0"))
#' @return Named list of ggplot treemaps, one per cluster
plot_enrichment_treemap <- function(
  reduced_terms,
  ontology_pal = c("MF" = "#FF6F00", "CC" = "#C71B00", "BP" = "#018EA0")
) {
  if (!requireNamespace("treemapify", quietly = TRUE)) {
    stop(
      "The 'treemapify' package is required for this function. Please install it using install.packages('treemapify')."
    )
  }
  clusters <- sort(unique(reduced_terms[["cluster"]]))
  plots <- lapply(clusters, function(cluster_id) {
    plot_data <- reduced_terms |>
      dplyr::filter(!!rlang::sym("cluster") == cluster_id) |>
      dplyr::mutate(
        ontology = factor(
          !!rlang::sym("ontology"),
          levels = names(ontology_pal)
        )
      )
    ggplot2::ggplot(
      plot_data,
      ggplot2::aes(
        area = !!rlang::sym("score"),
        subgroup = !!rlang::sym("parentTerm"),
        fill = !!rlang::sym("ontology")
      )
    ) +
      treemapify::geom_treemap(color = "black", show.legend = TRUE) +
      treemapify::geom_treemap_subgroup_border(color = "black") +
      treemapify::geom_treemap_text(
        ggplot2::aes(label = !!rlang::sym("term")),
        colour = "black",
        place = "centre",
        grow = TRUE,
        min.size = 0,
        reflow = TRUE
      ) +
      treemapify::geom_treemap_subgroup_text(
        place = "centre",
        grow = TRUE,
        reflow = TRUE,
        alpha = 1,
        colour = "white",
        fontface = "bold",
        min.size = 0
      ) +
      ggplot2::scale_fill_manual(values = ontology_pal, na.value = "grey80") +
      ggplot2::theme_void() +
      ggplot2::ggtitle("Gene Ontology treemap") +
      ggplot2::theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = NULL))
  })
  names(plots) <- as.character(clusters)
  plots
}

#' Plot enrichment bubblemap for annotation results
#'
#' @param enrichment_results Data frame with annotation results (e.g. database_enrichment, kegg_enrichment, or go_enrichment)
#' @param odds_ratio_cap Maximum value for odds ratio to cap bubble size (default: 50)
#' @param size_range Numeric vector of length 2 for point size range (default: c(1, 4))
#' @param fill_palette Color palette for odds ratio (default: c("white", "orangered"))
#' @param facet_var Optional variable to facet by (default: NULL)
#' @return ggplot object (bubblemap)
#' @export
plot_enrichment_bubblemap <- function(
  enrichment_results,
  odds_ratio_cap = 50,
  size_range = c(1, 4),
  fill_palette = c("white", "orangered"),
  facet_var = NULL
) {
  parse_ratio <- function(ratio_str) {
    sapply(strsplit(ratio_str, "/"), function(x) {
      as.numeric(x[1]) / as.numeric(x[2])
    })
  }

  plot_data <- enrichment_results |>
    dplyr::mutate(
      GeneFrac = parse_ratio(!!rlang::sym("GeneRatio")),
      BgFrac = parse_ratio(!!rlang::sym("BgRatio")),
      odds_ratio = !!rlang::sym("GeneFrac") / !!rlang::sym("BgFrac"),
      capped_odds_ratio = pmin(
        !!rlang::sym("odds_ratio"),
        !!rlang::sym("odds_ratio_cap")
      ),
      cluster = as.factor(!!rlang::sym("Cluster ID")),
      Description = as.character(!!rlang::sym("Term"))
    )

  # Keep only the top 50 terms (by mean odds ratio across clusters)
  top_terms <- plot_data |>
    dplyr::group_by(!!rlang::sym("Description")) |>
    dplyr::summarise(
      mean_odds = mean(!!rlang::sym("odds_ratio"), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(!!rlang::sym("mean_odds"))) |>
    dplyr::slice_head(n = 150) |>
    dplyr::pull(!!rlang::sym("Description"))

  plot_data <- plot_data |>
    dplyr::filter(!!rlang::sym("Description") %in% top_terms)

  odds_mat <- plot_data |>
    dplyr::group_by(!!rlang::sym("cluster"), !!rlang::sym("Description")) |>
    dplyr::summarise(
      odds_ratio = max(!!rlang::sym("odds_ratio"), na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(
      names_from = !!rlang::sym("Description"),
      values_from = !!rlang::sym("odds_ratio"),
      values_fill = 0
    ) |>
    tibble::column_to_rownames("cluster") |>
    as.matrix()

  # Cluster rows (clusters) and columns (terms)
  row_order <- rownames(odds_mat)[stats::hclust(stats::dist(odds_mat))[[
    "order"
  ]]]
  col_order <- colnames(odds_mat)[stats::hclust(stats::dist(t(odds_mat)))[[
    "order"
  ]]]

  plot_data <- plot_data |>
    dplyr::mutate(
      cluster = factor(!!rlang::sym("cluster"), levels = row_order),
      Description = factor(!!rlang::sym("Description"), levels = col_order),
      Description = ifelse(
        nchar(as.character(!!rlang::sym("Description"))) > 50,
        paste0(substr(as.character(!!rlang::sym("Description")), 1, 50), "..."),
        as.character(!!rlang::sym("Description"))
      )
    )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = !!rlang::sym("Description"),
      y = !!rlang::sym("cluster"),
      size = !!rlang::sym("capped_odds_ratio"),
      fill = !!rlang::sym("capped_odds_ratio")
    )
  ) +
    ggplot2::geom_point(shape = 21, color = "black") +
    ggplot2::scale_size_continuous(
      name = NULL,
      range = size_range,
      limits = c(0, odds_ratio_cap),
      breaks = rev(pretty(c(0, odds_ratio_cap)))
    ) +
    ggplot2::scale_fill_gradient(
      name = "Odds ratio",
      low = fill_palette[1],
      high = fill_palette[2],
      limits = c(0, odds_ratio_cap)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = -90, hjust = 0, vjust = 0.5),
      axis.title = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (!is.null(facet_var) && facet_var %in% names(plot_data)) {
    p <- p + ggplot2::facet_wrap(facet_var, scales = "free_y")
  }

  return(p)
}
