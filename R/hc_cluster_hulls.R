#' Estimate 2D density
#'
#' @param x x coordinates
#' @param y y coordinates
#' @param h bandwidth
#' @param n number of grid points
#' @param lims limits
#'
#' @returns Tidy tibble with x, y, z (density), x_coord, y_coord
#' @keywords internal
get_density <-
  function(x, y, h = 0.5, n = 100, lims = c(range(x), range(y))) {
    g_density <-
      MASS::kde2d(x, y, h = h, n = n, lims = lims)

    g_density[["z"]] |>
      tibble::as_tibble(rownames = "x") |>
      tidyr::gather(!!rlang::sym("y"), !!rlang::sym("z"), -!!rlang::sym("x")) |>
      dplyr::mutate(y = gsub("V", "", !!rlang::sym("y"))) |>
      dplyr::mutate_at(
        dplyr::vars(!!rlang::sym("x"), !!rlang::sym("y")),
        as.integer
      ) |>
      dplyr::left_join(
        tibble::enframe(g_density[["x"]], "x", "x_coord"),
        by = "x"
      ) |>
      dplyr::left_join(
        tibble::enframe(g_density[["y"]], "y", "y_coord"),
        by = "y"
      )
  }

#' Calculate UMAP cluster hulls
#'
#' `hc_cluster_hulls()` computes density-based cluster hulls on UMAP embeddings stored within an AnnDatR object. It identifies subclusters using DBSCAN, estimates density, filters based on cumulative density and fraction limits, and generates concave hulls for each cluster.
#'
#' @param AnnDatR An AnnDatR object containing UMAP embeddings and cluster assignments.
#' @param n Number of grid points for density estimation (default is 300).
#' @param cum_z_lim Cumulative density limit for filtering (default is 0.95).
#' @param frac_lim Minimum fraction limit for subcluster inclusion (default is 0.05).
#' @param plot_range Range of the plot in the format c(xmin, xmax, ymin, ymax) (default is 1.05 times the min and max of UMAP1 and UMAP2).
#' @param poly_concavity Controls the concavity of the polygons (default is 1). Higher values result in less detailed polygons.
#' @param poly_smoothing Controls the smoothing of the polygons (default is 1). Higher values result in less detailed polygons.
#' @param relative_bandwidth Relative bandwidth for density estimation as a fraction of the plot diagonal (default is 1/200).
#'
#' @returns AnnDatR object with UMAP cluster hulls within the AnnDatR object.
#'
#' @export
#' @examples
#' # Calculate UMAP cluster hulls on an AnnDatR object after UMAP and clustering
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 6.3)
#' adata_res <- hc_umap(adata_res)
#' adata_res <- hc_cluster_hulls(adata_res)
#' head(adata_res$uns$UMAP_hulls$hulls)
hc_cluster_hulls <-
  function(
    AnnDatR,
    n = 300,
    cum_z_lim = 0.95,
    frac_lim = 0.05,
    plot_range = rep(
      c(
        min(c(AnnDatR[["obs"]][["UMAP1"]], AnnDatR[["obs"]][["UMAP2"]])),
        max(c(AnnDatR[["obs"]][["UMAP1"]], AnnDatR[["obs"]][["UMAP2"]]))
      ),
      2
    ) *
      1.05,
    poly_concavity = 1,
    poly_smoothing = 1,
    relative_bandwidth = 1 / 200
  ) {
    if (!requireNamespace("fpc", quietly = TRUE)) {
      stop(
        "The 'fpc' package is required for this function. Please install it using install.packages('fpc')."
      )
    }
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop(
        "The 'MASS' package is required for this function. Please install it using install.packages('MASS')."
      )
    }
    if (!requireNamespace("concaveman", quietly = TRUE)) {
      stop(
        "The 'concaveman' package is required for this function. Please install it using install.packages('concaveman')."
      )
    }
    if (
      is.null(AnnDatR[["obs"]][["UMAP1"]]) &&
        is.null(AnnDatR[["obs"]][["UMAP2"]])
    ) {
      stop(
        "UMAP1 or 2 was not found in AnnDatR$obs. Call `hc_umap()` before `hc_cluster_hulls()`."
      )
    }
    if (is.null(AnnDatR[["obs"]][["cluster"]])) {
      stop(
        "AnnDatR$obs$cluster not found. Call `hc_cluster_consensus()` before `hc_cluster_hulls()`."
      )
    }
    V1 <- AnnDatR[["obs"]][["UMAP1"]]
    V2 <- AnnDatR[["obs"]][["UMAP2"]]
    element_id <- AnnDatR[["obs"]][["ensembl_id"]]
    cluster_membership <- AnnDatR[["obs"]][["cluster"]]

    # Combine input data
    cluster_data <-
      tibble::tibble(V1, V2, element_id, cluster = cluster_membership)

    # Calculate plot range
    plot_range_tb <-
      (\(x) {
        names(x) <- c("xmin", "xmax", "ymin", "ymax")
        x
      })(plot_range) |>
      tibble::enframe() |>
      tidyr::spread(!!rlang::sym("name"), !!rlang::sym("value"))

    # Calculate diagonal length
    plot_diagonal <-
      sqrt(
        (plot_range_tb[["xmax"]] - plot_range_tb[["xmin"]])^2 +
          (plot_range_tb[["ymax"]] - plot_range_tb[["ymin"]])^2
      )

    # Set bandwidth to fraction of diagonal
    plot_bandwidth <-
      relative_bandwidth * plot_diagonal

    # Find subclusters
    subclusters <-
      cluster_data |>
      dplyr::group_by(!!rlang::sym("cluster")) |>
      dplyr::mutate(
        n_cluster_genes = dplyr::n_distinct(!!rlang::sym("element_id")),
        sub_cluster = {
          dbscan_result <- fpc::dbscan(data.frame(V1, V2), eps = plot_bandwidth)
          dbscan_result[["cluster"]]
        }
      ) |>
      dplyr::ungroup() |>
      dplyr::group_by(!!rlang::sym("cluster"), !!rlang::sym("sub_cluster")) |>
      dplyr::mutate(
        n_sub_genes = dplyr::n_distinct(!!rlang::sym("element_id"))
      ) |>
      dplyr::ungroup()

    subclusters_classes <-
      subclusters |>
      dplyr::select(dplyr::any_of(c(
        "cluster",
        "sub_cluster",
        "n_cluster_genes",
        "n_sub_genes"
      ))) |>
      dplyr::distinct() |>
      dplyr::group_by(!!rlang::sym("cluster")) |>
      dplyr::mutate(
        n_sub_genes = ifelse(
          !!rlang::sym("sub_cluster") == 0,
          0,
          !!rlang::sym("n_sub_genes")
        ),
        sub_type = dplyr::case_when(
          !!rlang::sym("sub_cluster") == 0 ~ "outlier",
          !!rlang::sym("n_sub_genes") /
            !!rlang::sym("n_cluster_genes") < frac_lim ~ "outlier",
          rank(-n_sub_genes, ties.method = "first") == 1 ~ "primary",
          T ~ "secondary"
        )
      ) |>
      dplyr::select(dplyr::any_of(c(
        "cluster",
        "sub_cluster",
        "sub_type"
      )))

    subclusters_classed <-
      subclusters |>
      dplyr::left_join(subclusters_classes, by = c("cluster", "sub_cluster"))

    # Calculate plot density
    plot_density <-
      subclusters_classed |>
      dplyr::filter(!!rlang::sym("sub_type") != "outlier") |>
      dplyr::group_by(
        !!rlang::sym("cluster"),
        !!rlang::sym("sub_cluster"),
        !!rlang::sym("sub_type")
      ) |>
      tidyr::nest() |>
      dplyr::mutate(
        density_data = lapply(
          !!rlang::sym("data"),
          function(df) {
            get_density(
              df$V1,
              df$V2,
              h = plot_bandwidth,
              n = n,
              lims = plot_range
            )
          }
        )
      ) |>
      tidyr::unnest(!!rlang::sym("density_data")) |>
      dplyr::ungroup() |>
      dplyr::filter(!!rlang::sym("z") > 1e-200) |>
      dplyr::group_by(!!rlang::sym("cluster"), !!rlang::sym("sub_cluster")) |>
      dplyr::mutate(z = !!rlang::sym("z") / sum(!!rlang::sym("z"))) |>
      dplyr::arrange(
        !!rlang::sym("cluster"),
        !!rlang::sym("sub_cluster"),
        -!!rlang::sym("z")
      ) |>
      dplyr::mutate(cum_z = cumsum(!!rlang::sym("z"))) |>
      dplyr::ungroup()

    # Filter pixels such that 95% of density is included
    # Each point is then assigned to the cluster with highest density
    plot_density_filtered <-
      plot_density |>
      dplyr::filter(!!rlang::sym("cum_z") < cum_z_lim) |>
      dplyr::group_by(!!rlang::sym("x"), !!rlang::sym("y")) |>
      dplyr::top_n(1, !!rlang::sym("z")) |>
      dplyr::slice(1) |>
      dplyr::ungroup()

    # Calculate size of landmass
    plot_density_landmass <-
      plot_density_filtered |>
      dplyr::group_by(!!rlang::sym("cluster"), !!rlang::sym("sub_cluster")) |>
      dplyr::mutate(
        landmass = {
          dbscan_result <- fpc::dbscan(
            data.frame(
              !!rlang::sym("x_coord"),
              !!rlang::sym("y_coord")
            ),
            eps = plot_bandwidth
          )
          dbscan_result[["cluster"]]
        }
      ) |>
      dplyr::group_by(
        !!rlang::sym("cluster"),
        !!rlang::sym("sub_cluster"),
        !!rlang::sym("landmass")
      ) |>
      dplyr::mutate(n_landmass_points = length(!!rlang::sym("x"))) |>
      dplyr::ungroup() |>
      dplyr::group_by(!!rlang::sym("cluster")) |>
      dplyr::mutate(n_total_points = length(!!rlang::sym("x"))) |>
      dplyr::ungroup()

    # Classify landmasses
    plot_density_landmass_classes <-
      plot_density_landmass |>
      dplyr::select(
        dplyr::any_of(c(
          "cluster",
          "sub_cluster",
          "landmass",
          "n_landmass_points",
          "n_total_points"
        ))
      ) |>
      dplyr::distinct() |>
      dplyr::group_by(!!rlang::sym("cluster"), !!rlang::sym("sub_cluster")) |>
      dplyr::mutate(
        frac_landmass = !!rlang::sym("n_landmass_points") /
          !!rlang::sym("n_total_points"),
        landmass_type = dplyr::case_when(
          rank(-!!rlang::sym("n_landmass_points"), ties.method = "first") ==
            1 ~ "primary",
          T ~ "secondary"
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::select(
        dplyr::any_of(c(
          "cluster",
          "sub_cluster",
          "landmass",
          "landmass_type",
          "frac_landmass"
        ))
      )

    plot_density_landmass_classed <-
      plot_density_landmass |>
      dplyr::left_join(
        plot_density_landmass_classes,
        by = c("cluster", "sub_cluster", "landmass")
      ) |>
      dplyr::arrange(!!rlang::sym("cluster"))

    plot_density_mainland_filtered <-
      plot_density_landmass_classed |>
      dplyr::filter(!!rlang::sym("frac_landmass") > frac_lim) |>
      dplyr::filter(!!rlang::sym("cum_z") < cum_z_lim) |>
      dplyr::group_by(!!rlang::sym("x"), !!rlang::sym("y")) |>
      dplyr::top_n(1, !!rlang::sym("z")) |>
      dplyr::slice(1) |>
      dplyr::ungroup()

    ######

    # Create polygons
    # Poly smoothing: How small distances should be further detailed -
    # higher values --> less detailed
    # poly concavity: How convex polygons should be -
    # higher values --> less detailed

    plot_data_hulls <-
      plot_density_mainland_filtered |>
      dplyr::group_by(
        !!rlang::sym("cluster"),
        !!rlang::sym("sub_cluster"),
        !!rlang::sym("landmass"),
        !!rlang::sym("sub_type")
      ) |>
      tidyr::nest() |>
      dplyr::mutate(
        coords = lapply(
          !!rlang::sym("data"),
          function(.x) {
            as.matrix(data.frame(
              x = .x[["x_coord"]],
              y = .x[["y_coord"]]
            )) |>
              concaveman::concaveman(
                concavity = poly_concavity,
                length_threshold = plot_bandwidth * poly_smoothing
              ) |>
              tibble::as_tibble()
          }
        )
      ) |>
      tidyr::unnest(!!rlang::sym("coords")) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        polygon_id = paste(
          !!rlang::sym("cluster"),
          !!rlang::sym("sub_cluster"),
          !!rlang::sym("landmass"),
          sep = "_"
        )
      ) |>
      dplyr::select(-dplyr::any_of(c("data"))) |>
      dplyr::rename(X = "V1", Y = "V2")

    plot_density_center <-
      plot_density |>
      dplyr::filter(!!rlang::sym("sub_type") == "primary") |>
      dplyr::group_by(!!rlang::sym("cluster")) |>
      dplyr::top_n(1, !!rlang::sym("z")) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::select(
        !!rlang::sym("cluster"),
        x = !!rlang::sym("x_coord"),
        y = !!rlang::sym("y_coord")
      )

    AnnDatR_out <- AnnDatR$clone(deep = TRUE)
    AnnDatR_out[["uns"]][["UMAP_hulls"]] <- list(
      hulls = plot_data_hulls,
      density = plot_density,
      landmass_pixels = plot_density_mainland_filtered,
      center_density = plot_density_center
    )

    AnnDatR_out[["obs"]] <- AnnDatR_out[["obs"]] |>
      dplyr::select(-dplyr::any_of(c("UMAP1", "UMAP2")))

    return(AnnDatR_out)
  }
