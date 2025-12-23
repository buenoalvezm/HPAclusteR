#' Clusters genes using specified method and parameters
#'
#' @param genes Gene names
#' @param neighbors Neighbor graph object
#' @param method Clustering method: "louvain" or "leiden"
#' @param resolution Resolution parameter for clustering
#' @param seed Random seed for clustering
#'
#' @returns Data frame with gene cluster assignments
#' @keywords internal
cluster_genes <- function(
  genes,
  neighbors,
  method = "louvain",
  resolution = 1,
  seed = seed
) {
  if (method == "louvain") {
    alg = 1
  } else if (method == "leiden") {
    alg = 4
  }

  suppressWarnings(
    {
      empty_matrix <- matrix(
        0,
        nrow = length(genes),
        ncol = length(genes),
        dimnames = list(genes, genes)
      )
      louv <-
        Seurat::CreateSeuratObject(assay = "Exp", counts = empty_matrix)

      # Create the Seurat object with the correct structure
      louv@graphs[["Exp_snn"]] <- neighbors[["snn"]]

      louv <- Seurat::FindClusters(
        louv,
        graph.name = "Exp_snn",
        resolution = resolution,
        algorithm = alg,
        random.seed = seed,
        verbose = FALSE
      )
    }
  )

  col <- paste("Exp_snn_res.", resolution, sep = "")
  res <- as.numeric(as.character(louv@meta.data[, col]))
  res <- data.frame(gene = genes, cluster = res)

  return(res)
}

#' Finds consensus clustering from multiple clusterings
#'
#' @param all_clusterings Data frame containing clustering results from multiple runs
#' @param n Number of clustering runs
#' @param seed Random seed for consensus clustering
#' @param verbose Logical indicating whether to print progress messages
#'
#' @returns Consensus clustering results
#' @keywords internal
find_consensus <- function(
  all_clusterings,
  n,
  seed = 42,
  verbose = TRUE
) {
  # Create ensemble of clusterings
  ensemble <- lapply(c(1:n), function(x) {
    all_clusterings |>
      dplyr::select(!!rlang::sym("gene"), paste("seed_", x, sep = "")) |>
      tibble::deframe() |>
      clue::as.cl_partition()
  }) |>
    clue::as.cl_ensemble()

  # Determine number of clusters for consensus
  num_clusters <- lapply(c(1:n), function(x) {
    all_clusterings |>
      dplyr::pull(paste("seed_", x, sep = "")) |>
      dplyr::n_distinct()
  })

  k <- as.numeric(num_clusters) |>
    stats::median() |>
    floor()

  set.seed(seed)

  # Compute consensus clustering
  cons_clustering <- clue::cl_consensus(
    ensemble,
    method = "SE",
    control = list(k = k, nruns = 1, verbose = verbose)
  )

  # Extract final clustering assignments
  final_clustering <-
    clue::cl_class_ids(cons_clustering) |>
    tibble::enframe(name = "gene", value = "cluster")

  # Handle empty clusters
  empty_clusters <-
    final_clustering |>
    dplyr::group_by(!!rlang::sym("cluster")) |>
    dplyr::mutate(size = dplyr::n_distinct(!!rlang::sym("gene"))) |>
    dplyr::filter(!!rlang::sym("size") < 5) |>
    dplyr::select(-!!rlang::sym("size")) |>
    dplyr::mutate(cluster = as.numeric(!!rlang::sym("cluster")))

  if (dim(empty_clusters)[1] != 0) {
    to_rename <- lapply(seq_len(nrow(empty_clusters)), function(i) {
      current <- empty_clusters[i, ]

      probabilities <- cons_clustering$.Data[,] |>
        as.data.frame() |>
        tibble::as_tibble(rownames = "gene") |>
        dplyr::filter(!!rlang::sym("gene") == current[["gene"]]) |>
        dplyr::select(-dplyr::any_of(c("gene"))) |>
        tidyr::gather(!!rlang::sym("cluster"), !!rlang::sym("probability")) |>
        dplyr::arrange(-!!rlang::sym("probability")) |>
        dplyr::mutate(
          cluster = as.numeric(sub("V", "", !!rlang::sym("cluster")))
        ) |>
        dplyr::filter(!!rlang::sym("cluster") != current[["cluster"]])

      if (nrow(probabilities) > 0) {
        data.frame(
          gene = current[["gene"]],
          new_cluster = probabilities[[1, 1]]
        )
      } else {
        data.frame(
          gene = current[["gene"]],
          new_cluster = NA
        )
      }
    }) |>
      dplyr::bind_rows()

    final_clustering_corrected <-
      final_clustering |>
      dplyr::mutate(cluster = as.numeric(!!rlang::sym("cluster"))) |>
      dplyr::left_join(to_rename) |>
      dplyr::mutate(
        new_cluster = dplyr::if_else(
          is.na(!!rlang::sym("new_cluster")),
          !!rlang::sym("cluster"),
          !!rlang::sym("new_cluster")
        )
      )

    mapping_table <-
      final_clustering_corrected |>
      dplyr::select(dplyr::any_of(c("new_cluster"))) |>
      dplyr::distinct() |>
      dplyr::arrange(!!rlang::sym("new_cluster")) |>
      tibble::rownames_to_column("renumbered_cluster")

    final_clustering <-
      final_clustering_corrected |>
      dplyr::left_join(mapping_table) |>
      dplyr::select(-dplyr::any_of(c("new_cluster"))) |>
      dplyr::rename(
        cons_cluster = !!rlang::sym("cluster"),
        cluster = !!rlang::sym("renumbered_cluster")
      ) |>
      dplyr::mutate(
        cluster = as.character(!!rlang::sym("cluster")),
        cons_cluster = as.character(
          !!rlang::sym("cons_cluster")
        )
      )
  } else {
    mapping_table <-
      final_clustering |>
      dplyr::select(new_cluster = !!rlang::sym("cluster")) |>
      dplyr::mutate(new_cluster = as.character(!!rlang::sym("new_cluster"))) |>
      dplyr::distinct() |>
      dplyr::arrange(as.numeric(!!rlang::sym("new_cluster"))) |>
      tibble::rownames_to_column("renumbered_cluster")

    final_clustering <-
      final_clustering |>
      dplyr::mutate(cons_cluster = as.character(!!rlang::sym("cluster"))) |>
      dplyr::select(-dplyr::any_of(c("cluster"))) |>
      dplyr::left_join(
        mapping_table,
        by = c("cons_cluster" = "new_cluster")
      ) |>
      dplyr::rename(cluster = !!rlang::sym("renumbered_cluster"))
  }

  # Extract cluster membership matrix
  cons_matrix <-
    cons_clustering$.Data[,] |>
    as.data.frame() |>
    tibble::as_tibble(rownames = "gene") |>
    tidyr::gather(!!rlang::sym("cluster"), !!rlang::sym("membership"), -1) |>
    dplyr::filter(!!rlang::sym("membership") > 0) |>
    dplyr::mutate(
      cluster = as.character(gsub("V", "", !!rlang::sym("cluster")))
    ) |>
    dplyr::rename(cons_cluster = !!rlang::sym("cluster")) |>
    dplyr::left_join(
      mapping_table |>
        dplyr::mutate(
          new_cluster = as.character(!!rlang::sym("new_cluster"))
        ) |>
        dplyr::rename(
          cons_cluster = !!rlang::sym("new_cluster"),
          cluster = !!rlang::sym("renumbered_cluster")
        )
    ) |>
    dplyr::filter(!is.na(!!rlang::sym("cluster")))

  return(list(
    consensus_clustering = final_clustering,
    membership_matrix = cons_matrix
  ))
}

#' Consensus clustering of genes
#'
#' Clusters genes using multiple runs of a clustering algorithm with different seeds and resolutions, then derives a consensus clustering.
#'
#' @param AnnDatR An AnnDatR object containing the data with SNN graph.
#' @param resolution Resolution parameter for clustering (default is 6).
#' @param method Clustering method to use: "louvain" (default) or "leiden".
#' @param n_seeds Number of different random seeds to use for clustering (default is 100).
#' @param seed Random seed for reproducibility (default is 42).
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#'
#' @returns Consensus clustering results stored within the AnnDatR object.
#'
#' @export
#' @examples
#' # Perform consensus clustering on an AnnDatR object after SNN graph computation
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 7)
#' head(adata_res$uns$consensus_clustering)
#' head(adata_res$obs)
hc_cluster_consensus <- function(
  AnnDatR,
  resolution = 6,
  method = "louvain",
  n_seeds = 100,
  seed = 42,
  verbose = TRUE
) {
  if (is.null(AnnDatR[["uns"]][["neighbors"]])) {
    stop(
      "AnnDatR$uns$neighbors not found. Call `hc_snn()` before `hc_cluster_consensus()`."
    )
  }
  seeds = 1:n_seeds

  cluster_data <-
    # Create all combinations of resolution and seed
    tidyr::crossing(resolution = c(resolution), seed = seeds) |>
    # For each combination, perform clustering with the respective parameters
    dplyr::rowwise() |>
    dplyr::mutate(
      result = list(
        cluster_genes(
          genes = AnnDatR[["obs_names"]],
          neighbors = AnnDatR[["uns"]][["neighbors"]],
          method = method,
          resolution = !!rlang::sym("resolution"),
          seed = !!rlang::sym("seed")
        )
      )
    ) |>
    dplyr::ungroup() |>
    tidyr::unnest(!!rlang::sym("result"))

  # Process clustering results
  cluster_data <- cluster_data |>
    dplyr::mutate(cluster = as.character(!!rlang::sym("cluster") + 1)) |>
    dplyr::group_by(!!rlang::sym("resolution"), !!rlang::sym("seed")) |>
    dplyr::mutate(k = dplyr::n_distinct(!!rlang::sym("cluster"))) |>
    dplyr::ungroup() |>
    dplyr::group_by(resolution) |>
    dplyr::mutate(resolution_k = ceiling(stats::median(!!rlang::sym("k")))) |>
    dplyr::ungroup() |>
    dplyr::filter(
      !!rlang::sym("resolution_k") >= 30,
      !!rlang::sym("resolution_k") <= 110
    ) |>
    dplyr::select(-dplyr::any_of(c("k", "resolution_k")))

  if (cluster_data |> nrow() == 0) {
    stop(
      "No clusterings found with the specified resolution parameter. Try increasing the `resolution` argument."
    )
  }
  cluster_data <- cluster_data |>
    tidyr::pivot_wider(
      names_from = seed,
      values_from = !!rlang::sym("cluster"),
      names_prefix = "seed_"
    )

  # Create consensus clustering
  cluster_consensus <- find_consensus(
    cluster_data,
    n = n_seeds,
    seed = seed,
    verbose = verbose
  )

  # Store results in AnnDatR
  AnnDatR_out <- AnnDatR$clone(deep = TRUE)
  AnnDatR_out[["uns"]][["cluster_data"]] <- cluster_data

  AnnDatR_out[["uns"]][["mapping_table"]] <-
    cluster_consensus[["consensus_clustering"]] |>
    dplyr::select(dplyr::any_of(c("cluster", "cons_cluster"))) |>
    dplyr::distinct() |>
    dplyr::arrange(as.numeric(!!rlang::sym("cluster")))

  AnnDatR_out[["uns"]][[
    "consensus_clustering"
  ]] <- cluster_consensus[["consensus_clustering"]] |>
    dplyr::select(-dplyr::any_of(c("cons_cluster")))
  AnnDatR_out[["uns"]][[
    "membership_matrix"
  ]] <- cluster_consensus[["membership_matrix"]] |>
    dplyr::select(-dplyr::any_of(c("cons_cluster")))

  AnnDatR_out[["obs"]] <- AnnDatR_out[["obs"]] |>
    dplyr::left_join(
      AnnDatR_out[["uns"]][[
        "consensus_clustering"
      ]],
      by = dplyr::join_by(!!AnnDatR_out[["obs_names_col"]] == "gene")
    )

  names <- AnnDatR_out[["obs"]][["cluster"]] |>
    unique() |>
    as.double() |>
    sort() |>
    as.character()

  base_colors <- c(
    "#1F83B4",
    "#BCBD22",
    "#FF7F0E",
    "#D62728",
    "#9467BD",
    "#8C564B",
    "#d7b1cbff",
    "#7F7F7F",
    "#1cd5ea",
    "#1f731f",
    "#ff4dac",
    "#370054ff",
    "#FFD700"
  )
  cluster_colors <- base_colors[rep(
    seq_along(base_colors),
    length.out = length(names)
  )]

  cluster_colors <- as.data.frame(cluster_colors, (names)) |>
    tibble::rownames_to_column("cluster")

  AnnDatR_out[["obs"]] <- AnnDatR_out[["obs"]] |>
    dplyr::left_join(cluster_colors, by = "cluster")

  return(AnnDatR_out)
}
