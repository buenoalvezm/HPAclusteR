#' Convert a shared nearest neighbour graph to an igraph object
#'
#' @param snn Sparse, symmetric, weighted adjacency matrix.
#'
#' @returns An undirected weighted `igraph` graph.
#' @keywords internal
#' @noRd
snn_to_igraph <- function(snn) {
  igraph::graph_from_adjacency_matrix(
    snn,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
}

#' Clusters genes using specified method and parameters
#'
#' Community detection on the shared nearest neighbour graph. Both algorithms
#' optimise modularity at the requested resolution and both depend on the
#' random number generator, so different seeds give different partitions --
#' which is what makes the consensus step meaningful.
#'
#' @param genes Gene names, in the row order of the graph.
#' @param graph An `igraph` graph built from the SNN adjacency matrix.
#' @param method Clustering method: `"louvain"` or `"leiden"`.
#' @param resolution Resolution parameter for clustering.
#' @param seed Random seed for clustering.
#'
#' @returns Data frame with gene cluster assignments, numbered from zero.
#' @keywords internal
cluster_genes <- function(
  genes,
  graph,
  method = c("louvain", "leiden"),
  resolution = 1,
  seed = 42
) {
  method <- match.arg(method)
  weights <- igraph::edge_attr(graph, "weight")

  set.seed(seed)
  membership <- switch(
    method,
    louvain = igraph::membership(
      igraph::cluster_louvain(graph, weights = weights, resolution = resolution)
    ),
    leiden = igraph::membership(
      igraph::cluster_leiden(
        graph,
        weights = weights,
        resolution = resolution,
        objective_function = "modularity",
        n_iterations = 2L
      )
    )
  )

  data.frame(
    gene = genes,
    # Number from zero, so the downstream `+ 1` keeps the historical labels.
    cluster = as.integer(membership) - 1L,
    stringsAsFactors = FALSE
  )
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
  
  if (nrow(empty_clusters) > 0) {
    empty_cluster_ids <- unique(empty_clusters[["cluster"]])
    
    to_rename <- lapply(seq_len(nrow(empty_clusters)), function(i) {
      current <- empty_clusters[i, ]

      probabilities <- cons_clustering$.Data[,] |>
        as.data.frame() |>
        tibble::as_tibble(rownames = "gene") |>
        dplyr::filter(!!rlang::sym("gene") == current[["gene"]]) |>
        dplyr::select(-dplyr::any_of(c("gene"))) |>
        tidyr::gather(!!rlang::sym("cluster"), !!rlang::sym("probability")) |>
        dplyr::arrange(dplyr::desc(!!rlang::sym("probability")), !!rlang::sym("cluster")) |>
        dplyr::mutate(
          cluster = as.numeric(sub("V", "", !!rlang::sym("cluster")))
        ) |>
        dplyr::filter(!(!!rlang::sym("cluster") %in% empty_cluster_ids))

      if (nrow(probabilities) > 0) {
        data.frame(
          gene = current[["gene"]],
          new_cluster = probabilities[[1, "cluster"]]
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
      dplyr::left_join(to_rename, by = "gene") |>
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
      dplyr::left_join(mapping_table, by = "new_cluster") |>
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
    tidyr::gather(!!rlang::sym("cons_cluster"), !!rlang::sym("membership"), -1) |>
    dplyr::mutate(
      cons_cluster = as.character(gsub("V", "", !!rlang::sym("cons_cluster")))
    ) |>
    dplyr::left_join(
      mapping_table |>
        dplyr::mutate(
          new_cluster = as.character(!!rlang::sym("new_cluster"))
        ) |>
        dplyr::rename(
          cons_cluster = !!rlang::sym("new_cluster"),
          cluster = !!rlang::sym("renumbered_cluster")
        ),
      by = "cons_cluster"
    ) |>
    dplyr::filter(!is.na(!!rlang::sym("cluster"))) |>
    dplyr::left_join(
      final_clustering |> dplyr::select(
        !!rlang::sym("gene"), 
        hard_cluster = !!rlang::sym("cluster")
      ),
      by = "gene"
    ) |>
    dplyr::filter(
      !!rlang::sym("membership") > 0 | (!!rlang::sym("cluster") == !!rlang::sym("hard_cluster"))
    ) |>
    dplyr::select(-dplyr::any_of(c("hard_cluster", "cons_cluster")))

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
#' @param min_k Lower bound on the median number of communities across seeds
#'   (default is 30). Resolutions that produce fewer are rejected.
#' @param max_k Upper bound on the median number of communities across seeds
#'   (default is 110).
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#'
#' @details
#' The resolution is validated against `min_k` and `max_k`: if the median number
#' of communities across seeds falls outside that window, the function stops and
#' asks for a different `resolution`. These bounds were previously hard-coded.
#'
#' Re-running this function on an object that already carries a clustering
#' replaces the previous `cluster` and `cluster_colors` columns in `obs`.
#'
#' @returns Consensus clustering results stored within the AnnDatR object.
#'
#' @export
#' @examples
#' # Perform consensus clustering on an AnnDatR object after SNN graph computation
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 8, n_seeds = 20)
#' head(adata_res$uns$consensus_clustering)
#' head(adata_res$obs)
hc_cluster_consensus <- function(
  AnnDatR,
  resolution = 6,
  method = c("louvain", "leiden"),
  n_seeds = 100,
  seed = 42,
  min_k = 30,
  max_k = 110,
  verbose = TRUE
) {
  method <- match.arg(method)
  check_installed("igraph", "for community detection")

  if (is.null(AnnDatR[["uns"]][["neighbors"]])) {
    stop(
      "AnnDatR$uns$neighbors not found. Call `hc_snn()` before `hc_cluster_consensus()`.",
      call. = FALSE
    )
  }

  # Build the graph once and reuse it for every seed.
  graph <- snn_to_igraph(AnnDatR[["uns"]][["neighbors"]][["snn"]])
  # Read the column directly rather than through the `obs_names` active
  # binding: R6 active bindings are lost when an object is serialised, so a
  # binding restored from disk can report stale values.
  genes <- AnnDatR[["obs"]][[AnnDatR[["obs_names_col"]]]]

  if (verbose) {
    message(sprintf(
      "Running %s clustering at resolution %s across %d seeds.",
      method,
      format(resolution),
      n_seeds
    ))
  }

  cluster_data <- lapply(seq_len(n_seeds), function(current_seed) {
    partition <- cluster_genes(
      genes = genes,
      graph = graph,
      method = method,
      resolution = resolution,
      seed = current_seed
    )
    partition[["resolution"]] <- resolution
    partition[["seed"]] <- current_seed
    partition
  }) |>
    dplyr::bind_rows() |>
    dplyr::select(dplyr::any_of(c("resolution", "seed", "gene", "cluster")))

  # Label clusters from 1 and check that the resolution yields a usable number
  # of communities before spending time on the consensus step.
  cluster_data <- cluster_data |>
    dplyr::mutate(cluster = as.character(!!rlang::sym("cluster") + 1))

  k_per_seed <- cluster_data |>
    dplyr::group_by(!!rlang::sym("seed")) |>
    dplyr::summarise(k = dplyr::n_distinct(!!rlang::sym("cluster")), .groups = "drop") |>
    dplyr::pull(!!rlang::sym("k"))

  median_k <- ceiling(stats::median(k_per_seed))
  if (verbose) {
    message(sprintf(
      "Median number of communities across seeds: %d (range %d-%d).",
      median_k,
      min(k_per_seed),
      max(k_per_seed)
    ))
  }

  if (median_k < min_k || median_k > max_k) {
    stop(
      sprintf(
        paste0(
          "Resolution %s produces a median of %d clusters, outside the accepted range [%d, %d]. ",
          "%s the `resolution` argument, or widen `min_k`/`max_k`."
        ),
        format(resolution),
        median_k,
        min_k,
        max_k,
        if (median_k < min_k) "Increase" else "Decrease"
      ),
      call. = FALSE
    )
  }

  cluster_data <- cluster_data |>
    dplyr::select(dplyr::any_of(c("resolution", "seed", "gene", "cluster"))) |>
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

  # Drop any clustering left over from a previous call, so that re-running this
  # function on an already-clustered object replaces rather than duplicates it.
  AnnDatR_out[["obs"]] <- AnnDatR_out[["obs"]] |>
    dplyr::select(-dplyr::any_of(c("cluster", "cluster_colors"))) |>
    dplyr::left_join(
      AnnDatR_out[["uns"]][["consensus_clustering"]],
      by = dplyr::join_by(!!AnnDatR_out[["obs_names_col"]] == "gene")
    )

  cluster_levels <- AnnDatR_out[["obs"]][["cluster"]] |>
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
  # Recycle the palette so that every cluster gets a colour.
  cluster_colors <- tibble::tibble(
    cluster = cluster_levels,
    cluster_colors = rep(base_colors, length.out = length(cluster_levels))
  )

  AnnDatR_out[["obs"]] <- AnnDatR_out[["obs"]] |>
    dplyr::left_join(cluster_colors, by = "cluster")

  AnnDatR_out
}
