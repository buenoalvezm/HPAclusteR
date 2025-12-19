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

# #' Plot bubble plot for enrichment results
# #'
# #' @param enrichment_df Enrichment result dataframe (must have Cluster ID, Term, odds_ratio, `Adjusted P-value`)
# #' @param database Database to plot (optional)
# #' @return ggplot object
# plot_enrichment_bubble <- function(enrichment_df, database = NULL) {
#   df <- enrichment_df
#   if (!is.null(database)) {
#     df <- df[df$Database == database, ]
#   }
#   ggplot2::ggplot(
#     df,
#     ggplot2::aes(
#       x = Term,
#       y = `Cluster ID`,
#       size = odds_ratio,
#       fill = `Adjusted P-value`
#     )
#   ) +
#     ggplot2::geom_point(shape = 21, color = "black") +
#     ggplot2::scale_size(range = c(2, 10)) +
#     ggplot2::theme_bw() +
#     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
#     ggplot2::labs(
#       title = paste("Bubble plot:", database),
#       y = "Cluster",
#       x = "Term"
#     )
# }

# #' Plot heatmap for enrichment results
# #'
# #' @param enrichment_df Enrichment result dataframe (must have Cluster ID, Term, odds_ratio)
# #' @param database Database to plot (optional)
# #' @return ggplot object
# plot_enrichment_heatmap <- function(enrichment_df, database = NULL) {
#   df <- enrichment_df
#   if (!is.null(database)) {
#     df <- df[df$Database == database, ]
#   }
#   ggplot2::ggplot(
#     df,
#     ggplot2::aes(x = Term, y = factor(`Cluster ID`), fill = odds_ratio)
#   ) +
#     ggplot2::geom_tile() +
#     ggplot2::theme_minimal() +
#     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
#     ggplot2::labs(
#       title = paste("Heatmap:", database),
#       y = "Cluster",
#       x = "Term"
#     )
# }
