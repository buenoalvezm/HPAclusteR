#' Perform sample category gene classification (HPA logic, sample-agnostic)
#'
#' @param AnnDatR AnnDatR object containing gene expression data and sample annotations
#' @param sample_categories Name of the column in AnnDatR$var containing sample categories
#' @param cutoff_detected Expression value threshold to consider a gene as detected (default: 1)
#' @param lowest_exp_limit Minimum expression value to consider for calculations (default: 0.1)
#' @param group_enriched_limit Maximum number of sample categories for group enriched classification (default: NULL, auto-calculated)
#' @param fold Fold change threshold for classification (default: 4)
#'
#' @return A tibble with gene classification results including ENSG, specificity category, specific sample categories, tau value, and distribution category
#' @keywords internal
hc_classify_fun <- function(
  AnnDatR,
  sample_categories,
  cutoff_detected = 1,
  lowest_exp_limit = 0.1,
  group_enriched_limit = NULL,
  fold = 4
) {
  sample_id <- AnnDatR[["var_names_col"]]
  df <- AnnDatR[["X"]] |>
    tidyr::pivot_longer(
      cols = -sample_id,
      names_to = "ENSG",
      values_to = "Value"
    ) |>
    dplyr::left_join(
      AnnDatR[["var"]] |>
        dplyr::select(dplyr::any_of(c(sample_id, sample_categories))),
      by = sample_id
    )

  n_categories <- df |>
    dplyr::pull(!!rlang::sym(sample_categories)) |>
    unique() |>
    length()

  ge_limit <- if (is.null(group_enriched_limit)) {
    if (n_categories < 50) 5L else max(2L, round(n_categories * 0.1))
  } else {
    group_enriched_limit
  }

  df_sorted <- df |>
    dplyr::arrange(!!rlang::sym("ENSG"), dplyr::desc(!!rlang::sym("Value")))

  genes <- df_sorted |>
    dplyr::group_by(!!rlang::sym("ENSG")) |>
    dplyr::summarise(
      categories = list(get(sample_categories)),
      exps = list(!!rlang::sym("Value")),
      exp_avg = mean(!!rlang::sym("Value"), na.rm = TRUE),
      .groups = "drop"
    )

  out_ens <- character()
  out_cat <- character()
  out_categories <- character()
  out_tau <- numeric()

  for (r in seq_len(nrow(genes))) {
    ens <- genes[["ENSG"]][[r]]
    exps <- genes[["exps"]][[r]]
    categories <- genes[["categories"]][[r]]
    exp_avg <- genes[["exp_avg"]][[r]]

    if (length(categories) < 2) {
      next
    }

    log_exps <- log(exps + 1)
    if (length(log_exps) <= 1 || max(log_exps, na.rm = TRUE) == 0) {
      tau_val <- NA_real_
    } else {
      mx <- max(log_exps, na.rm = TRUE)
      tau_val <- sum(1 - (log_exps / mx), na.rm = TRUE) / (length(log_exps) - 1)
    }

    max_exp <- if (exps[[1]] > lowest_exp_limit) exps[[1]] else lowest_exp_limit
    second_max <- if (exps[[2]] > lowest_exp_limit) {
      exps[[2]]
    } else {
      lowest_exp_limit
    }

    if (max_exp < cutoff_detected) {
      out_ens <- c(out_ens, ens)
      out_cat <- c(out_cat, "Not detected")
      out_categories <- c(out_categories, NA_character_)
      out_tau <- c(out_tau, tau_val)
      next
    }

    if ((max_exp / second_max) >= fold) {
      out_ens <- c(out_ens, ens)
      out_cat <- c(out_cat, "Enriched")
      out_categories <- c(out_categories, categories[[1]])
      out_tau <- c(out_tau, tau_val)
      next
    }

    lim <- exps[[1]] / fold
    if (
      length(exps) >= 2 &&
        exps[[2]] >= cutoff_detected &&
        (length(exps) >= (ge_limit + 1) &&
          exps[[ge_limit + 1]] < lim)
    ) {
      group_size <- 1L
      max_check_idx <- min(ge_limit + 1L, length(exps))
      for (k in 2:max_check_idx) {
        if (exps[[k]] > lim) {
          group_size <- group_size + 1L
        } else {
          break
        }
      }

      if (
        group_size > 1L &&
          group_size <= ge_limit &&
          exps[[group_size]] >= cutoff_detected
      ) {
        mean_group <- sum(exps[1:group_size]) / group_size
        if ((group_size + 1L) <= length(exps)) {
          mo <- exps[[group_size + 1L]]
        } else {
          mo <- NA_real_
        }
        max_other <- if (!is.na(mo) && mo >= lowest_exp_limit) {
          mo
        } else {
          lowest_exp_limit
        }

        if ((mean_group / max_other) >= fold) {
          out_ens <- c(out_ens, ens)
          out_cat <- c(out_cat, "Group enriched")
          out_categories <- c(
            out_categories,
            paste(unique(categories[1:group_size]), collapse = ";")
          )
          out_tau <- c(out_tau, tau_val)
          next
        }
      }
    }

    enhanced <- character()
    for (idx in seq_along(exps)) {
      fp <- exps[[idx]]
      if (!is.na(fp) && fp >= exp_avg * fold && fp >= cutoff_detected) {
        enhanced <- c(enhanced, categories[[idx]])
      }
    }
    if (length(enhanced) > 0) {
      out_ens <- c(out_ens, ens)
      out_cat <- c(out_cat, "Enhanced")
      out_categories <- c(
        out_categories,
        paste(unique(enhanced), collapse = ";")
      )
      out_tau <- c(out_tau, tau_val)
      next
    }

    out_ens <- c(out_ens, ens)
    out_cat <- c(out_cat, "Low specificity")
    out_categories <- c(out_categories, NA_character_)
    out_tau <- c(out_tau, tau_val)
  }

  # Distribution classification
  dist_category <- vector("character", length(out_ens))
  for (i in seq_along(out_ens)) {
    ens <- out_ens[i]
    gene_df <- df |> dplyr::filter(!!rlang::sym("ENSG") == ens)
    n_detected <- sum(gene_df[["Value"]] >= cutoff_detected, na.rm = TRUE)
    n_total <- sum(!is.na(gene_df[["Value"]]))
    frac <- n_detected / n_total * 100
    if (n_detected == 0) {
      dist_category[i] <- "Not detected"
    } else if (frac == 100) {
      dist_category[i] <- "Detected in all"
    } else if (frac >= 31) {
      dist_category[i] <- "Detected in many"
    } else if (n_detected > 1) {
      dist_category[i] <- "Detected in some"
    } else if (n_detected == 1) {
      dist_category[i] <- "Detected in single"
    }
  }

  tibble::tibble(
    ENSG = out_ens,
    spec_category = out_cat,
    spec_sample_categories = out_categories,
    tau = round(out_tau, 2),
    dist_category = dist_category
  )
}

#' Plot specificity and distribution classification results
#'
#' @param class_tbl A tibble containing gene classification results with 'spec_category' and 'dist_category' columns
#'
#' @returns A combined pie chart plot showing the distribution of specificity and distribution categories
#'
#' @keywords internal
plot_specificity_distribution <- function(class_tbl) {
  # Ensure all levels are present
  spec_levels <- c("Enriched", "Group enriched", "Enhanced", "Low specificity")
  dist_levels <- c(
    "Detected in single",
    "Detected in some",
    "Detected in many",
    "Detected in all"
  )
  spec_colors <- c(
    "Enriched" = "#e41a1c",
    "Group enriched" = "#ff9900",
    "Enhanced" = "#984ea3",
    "Low specificity" = "#666666"
  )
  dist_colors <- c(
    "Detected in single" = "#a6dba0",
    "Detected in some" = "#38a6a5",
    "Detected in many" = "#3066b1",
    "Detected in all" = "#253494"
  )

  # Specificity table
  spec_tbl <- class_tbl |>
    dplyr::count(!!rlang::sym("spec_category")) |>
    dplyr::mutate(
      spec_category = factor(
        !!rlang::sym("spec_category"),
        levels = spec_levels
      )
    ) |>
    tidyr::complete(
      spec_category = spec_levels,
      fill = list(n = 0)
    ) |>
    dplyr::mutate(
      fraction = !!rlang::sym("n") / sum(!!rlang::sym("n")),
      ymax = cumsum(!!rlang::sym("fraction")),
      ymin = dplyr::lag(!!rlang::sym("ymax"), default = 0),
      label_pos = (!!rlang::sym("ymax") + !!rlang::sym("ymin")) / 2
    )

  # Distribution table
  dist_tbl <- class_tbl |>
    dplyr::count(!!rlang::sym("dist_category")) |>
    dplyr::mutate(
      dist_category = factor(
        !!rlang::sym("dist_category"),
        levels = dist_levels
      )
    ) |>
    tidyr::complete(dist_category = dist_levels, fill = list(n = 0)) |>
    dplyr::mutate(
      fraction = !!rlang::sym("n") / sum(!!rlang::sym("n")),
      ymax = cumsum(!!rlang::sym("fraction")),
      ymin = dplyr::lag(!!rlang::sym("ymax"), default = 0),
      label_pos = (!!rlang::sym("ymax") + !!rlang::sym("ymin")) / 2
    )

  # Pie chart for specificity
  p1 <- ggplot2::ggplot(spec_tbl) +
    ggplot2::geom_rect(
      ggplot2::aes(
        ymin = !!rlang::sym("ymin"),
        ymax = !!rlang::sym("ymax"),
        xmin = 0,
        xmax = 1,
        fill = !!rlang::sym("spec_category")
      ),
      color = "white"
    ) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(
      values = spec_colors,
      breaks = spec_levels,
      na.value = "grey80"
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(fill = NULL) +
    ggrepel::geom_text_repel(
      ggplot2::aes(
        x = 1,
        y = !!rlang::sym("label_pos"),
        label = !!rlang::sym("n")
      ),
      size = 4,
      direction = "y",
      nudge_x = 0.2,
      segment.size = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "vertical"
    )

  # Pie chart for distribution
  p2 <- ggplot2::ggplot(dist_tbl) +
    ggplot2::geom_rect(
      ggplot2::aes(
        ymin = !!rlang::sym("ymin"),
        ymax = !!rlang::sym("ymax"),
        xmin = 0,
        xmax = 1,
        fill = !!rlang::sym("dist_category")
      ),
      color = "white"
    ) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(
      values = dist_colors,
      breaks = dist_levels,
      na.value = "grey80"
    ) +
    ggplot2::theme_void() +
    ggplot2::labs(fill = NULL) +
    ggrepel::geom_text_repel(
      ggplot2::aes(
        x = 1,
        y = !!rlang::sym("label_pos"),
        label = !!rlang::sym("n")
      ),
      size = 4,
      direction = "y",
      nudge_x = 0.2,
      segment.size = 0.4,
      show.legend = FALSE
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "vertical"
    )

  p1 +
    p2 +
    patchwork::plot_annotation(
      title = "Specificity and distribution",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
    )
}

#' Perform sample category gene classification (HPA logic, sample-agnostic)
#'
#' `hc_classify` classifies genes based on their expression patterns across different sample categories
#' using a logic similar to that employed by the Human Protein Atlas (HPA). The function processes
#' each cluster of genes identified in the AnnDatR object and categorizes them into specificity categories
#' such as "Enriched", "Group enriched", "Enhanced", and "Low specificity". It also generates pie charts to visualize
#' the distribution of these categories within each cluster.
#'
#' @param AnnDatR AnnDatR object containing gene expression data and sample annotations
#' @param sample_categories Name of the column in AnnDatR$var containing sample categories
#' @param cutoff_detected Expression value threshold to consider a gene as detected (default: 1)
#' @param lowest_exp_limit Minimum expression value to consider for calculations (default: 0.1)
#' @param group_enriched_limit Maximum number of sample categories for group enriched classification (default: NULL, auto-calculated)
#' @param fold Fold change threshold for classification (default: 4)
#'
#' @return A tibble with gene classification results including ENSG, specificity category, specific sample categories, tau value, and distribution category
#' @export
#'
#' @examples
#' # Run clustering pipeline
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 7)
#'
#' # Classify genes based on sample categories
#' gene_classification <- hc_classify(adata_res, "tissue_name")
#' head(gene_classification$classification[["1"]])
#' gene_classification$pie_charts[["1"]]
hc_classify <- function(
  AnnDatR,
  sample_categories,
  cutoff_detected = 1,
  lowest_exp_limit = 0.1,
  group_enriched_limit = NULL,
  fold = 4
) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "The 'ggrepel' package is required for this function. Please install it using install.packages('ggrepel')."
    )
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "The 'patchwork' package is required for this function. Please install it using install.packages('patchwork')."
    )
  }
  clusters <- unique(AnnDatR[["obs"]][["cluster"]])
  class_res <- lapply(clusters, function(clust) {
    genes_in_cluster <- AnnDatR[["obs"]] |>
      dplyr::filter(!!rlang::sym("cluster") == clust) |>
      dplyr::pull(AnnDatR[["obs_names_col"]])
    AnnDatR_subset <- AnnDatR$clone(deep = TRUE)
    AnnDatR_subset[["X"]] <- AnnDatR_subset[["X"]] |>
      dplyr::select(dplyr::any_of(c(
        AnnDatR[["var_names_col"]],
        genes_in_cluster
      )))

    res <- hc_classify_fun(
      AnnDatR_subset,
      sample_categories,
      cutoff_detected = cutoff_detected,
      lowest_exp_limit = lowest_exp_limit,
      group_enriched_limit = group_enriched_limit,
      fold = fold
    )
    res
  })
  names(class_res) <- clusters
  class_res <- class_res[order(names(class_res))]

  pie_charts <- lapply(class_res, plot_specificity_distribution)
  names(pie_charts) <- names(class_res)

  return(list(
    classification = class_res,
    pie_charts = pie_charts
  ))
}
