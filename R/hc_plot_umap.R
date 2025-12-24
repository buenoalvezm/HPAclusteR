#' Plot UMAP with clusters and hulls
#'
#' `hc_plot_umap()` creates a UMAP plot from an AnnData object, displaying points, cluster hulls, or both.
#'
#' @param AnnDatR An AnnData object containing UMAP coordinates and cluster hulls.
#' @param plot Type of plot to create: "points" for scatter plot, "hulls" for cluster hulls, or "both" for overlaying points on hulls. Default is "points".
#'
#' @returns A ggplot2 object representing the UMAP visualization.
#'
#' @export
#' @examples
#' # Plot UMAP points from an AnnDatR object
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 7)
#' adata_res <- hc_umap(adata_res)
#' adata_res <- hc_cluster_hulls(adata_res)
#' hc_plot_umap(adata_res)
hc_plot_umap <- function(AnnDatR, plot = c("points", "hulls", "both")) {
  plot <- match.arg(plot)

  if (is.null(AnnDatR[["obsm"]][["X_umap"]])) {
    stop(
      "UMAP coordinates not found in AnnData object. Please run `hc_umap()` first."
    )
  }
  if (plot != "points" && is.null(AnnDatR[["uns"]][["UMAP_hulls"]])) {
    stop(
      "UMAP hulls not found in AnnData object. Please run `hc_cluster_hulls()` first."
    )
  }

  plot_data <- AnnDatR[["obs"]] |>
    dplyr::left_join(
      AnnDatR[["obsm"]][["X_umap"]] |>
        as.data.frame() |>
        tibble::rownames_to_column("ensembl_id") |>
        tibble::as_tibble(),
      by = dplyr::join_by(!!rlang::sym("ensembl_id"))
    )

  pal <- plot_data |>
    dplyr::select(dplyr::any_of(c("cluster", "cluster_colors"))) |>
    dplyr::arrange(!!rlang::sym("cluster")) |>
    dplyr::distinct() |>
    tibble::deframe()

  # Initialize the ggplot object
  umap_plot <- ggplot2::ggplot()

  # Add points if `plot` is "points" or "both"
  if (plot %in% c("points", "both")) {
    umap_plot <- umap_plot +
      ggplot2::geom_point(
        data = plot_data,
        ggplot2::aes(
          !!rlang::sym("UMAP_1"),
          !!rlang::sym("UMAP_2"),
          color = !!rlang::sym("cluster")
        ),
        inherit.aes = FALSE,
        size = 0.5,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_manual(values = pal)
  }

  # Add hulls if `plot` is "hulls" or "both"
  if (plot %in% c("hulls", "both")) {
    umap_plot <- umap_plot +
      ggplot2::geom_polygon(
        data = AnnDatR[["uns"]][["UMAP_hulls"]][["hulls"]],
        ggplot2::aes(
          !!rlang::sym("X"),
          !!rlang::sym("Y"),
          group = paste(
            !!rlang::sym("cluster"),
            !!rlang::sym("sub_cluster"),
            !!rlang::sym("landmass")
          ),
          fill = as.factor(!!rlang::sym("cluster")),
          alpha = 0.3
        ),
        show.legend = FALSE,
        color = "black"
      ) +
      ggplot2::scale_fill_manual(values = pal)
  }

  # Add center of density only if `plot` is "hulls"
  if (plot == "hulls") {
    umap_plot <- umap_plot +
      ggplot2::geom_point(
        data = AnnDatR[["uns"]][["UMAP_hulls"]][["center_density"]],
        ggplot2::aes(!!rlang::sym("x"), !!rlang::sym("y")),
        inherit.aes = FALSE,
        size = 0.2,
        show.legend = FALSE
      )
  }

  # Finalize the plot
  umap_plot <- umap_plot +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "UMAP1", y = "UMAP2") +
    theme_hc()

  return(umap_plot)
}
