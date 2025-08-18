run_consensus_clustering <-
  function(adata_t,
           resolution,
           num_seeds,
           clustering_method,
           inplace = TRUE) {
    seeds <- 1:num_seeds

    cluster_data <-
      tidyr::crossing(
        resolution = c(resolution),
        seed = seeds
      ) %>%
      group_by_all() %>%
      do({
        cluster_genes(
          genes = adata_t$obs_names(),
          neighbors = adata_t$uns$neighbors,
          clustering_method = clustering_method,
          resolution = .$resolution,
          seed = .$seed
        )
      }) %>%
      ungroup()

    cluster_data <- cluster_data %>%
      mutate(cluster = as.character(cluster + 1)) %>%
      group_by(resolution, seed) %>%
      mutate(k = n_distinct(cluster)) %>%
      ungroup() %>%
      group_by(resolution) %>%
      mutate(resolution_k = ceiling(median(k))) %>%
      ungroup() %>%
      filter(
        resolution_k >= 30,
        resolution_k <= 110
      ) %>%
      select(-k, -resolution_k) %>%
      spread(seed, cluster, sep = "_")

    cluster_consensus <-
      find_consensus(
        cluster_data,
        n = num_seeds,
        get_membership = T,
        runs = 1
      )


    if (inplace) {
      adata_t$uns$cluster_data <- cluster_data

      adata_t$uns$mapping_table <-
        cluster_consensus$consensus_clustering %>%
        select(cluster, cons_cluster) %>%
        distinct() %>%
        arrange(as.numeric(cluster))

      adata_t$uns$consensus_clustering <-
        cluster_consensus$consensus_clustering %>% select(-cons_cluster)
      adata_t$uns$membership_matrix <-
        cluster_consensus$membership_matrix %>% select(-cons_cluster)

      adata_t$obs <-
        adata_t$obs %>% left_join(adata_t$uns$consensus_clustering,
          by = join_by(!!adata_t$obs_names_col == "gene")
        )

      names <-
        adata_t$obs$cluster %>%
        unique() %>%
        as.double() %>%
        sort() %>%
        as.character()
      cluster_colors <-
        colorRampPalette(ggthemes::tableau_color_pal(palette = "Classic Cyclic")(13))(length(names))
      cluster_colors <-
        as.data.frame(cluster_colors, (names)) %>% rownames_to_column("cluster")
      adata_t$obs <-
        adata_t$obs %>% left_join(cluster_colors, by = "cluster")
    } else {
      return(c(
        "cluster_data" = cluster_data,
        "cluster_consensus" = cluster_consensus
      ))
    }
  }





find_consensus <- function(all_clusterings, n, get_membership = F, runs = 1) {
  require(clue)
  ensemble <- lapply(c(1:n), function(x) {
    all_clusterings %>%
      select(gene, paste("seed_", x, sep = "")) %>%
      deframe() %>%
      as.cl_partition()
  }) %>%
    as.cl_ensemble()

  num_clusters <- lapply(c(1:n), function(x) {
    all_clusterings %>%
      pull(paste("seed_", x, sep = "")) %>%
      n_distinct()
  })

  k <- as.numeric(num_clusters) %>%
    median() %>%
    floor()

  set.seed(1)

  cons_clustering <- cl_consensus(ensemble,
    method = "SE",
    control = list(
      k = k,
      nruns = runs,
      verbose = FALSE
    )
  )
  final_clustering <-
    cl_class_ids(cons_clustering) %>%
    enframe(name = "gene", value = "cluster")

  empty_clusters <-
    final_clustering %>%
    # mutate(cluster = as.numeric(cluster)) %>%
    group_by(cluster) %>%
    mutate(size = n_distinct(gene)) %>%
    filter(size < 5) %>%
    select(-size) %>%
    mutate(cluster = as.numeric(cluster))

  if (dim(empty_clusters)[1] != 0) {
    to_rename <-
      empty_clusters %>%
      group_by(gene) %>%
      do({
        current <- .

        probabilities <-
          cons_clustering$.Data[, ] %>%
          as_tibble(rownames = "gene") %>%
          filter(gene == current$gene) %>%
          select(-gene) %>%
          gather(cluster, probability) %>%
          arrange(-probability) %>%
          mutate(cluster = as.numeric(sub("V", "", cluster))) %>%
          filter(cluster != current$cluster)
        data.frame(
          gene = current$gene,
          new_cluster = probabilities[[1, 1]]
        )
      })

    final_clustering_corrected <-
      final_clustering %>%
      mutate(cluster = as.numeric(cluster)) %>%
      left_join(to_rename) %>%
      mutate(new_cluster = if_else(is.na(new_cluster), cluster, new_cluster))


    mapping_table <-
      final_clustering_corrected %>%
      select(new_cluster) %>%
      distinct() %>%
      arrange(new_cluster) %>%
      rownames_to_column("renumbered_cluster")

    final_clustering <-
      final_clustering_corrected %>%
      left_join(mapping_table) %>%
      select(-new_cluster) %>%
      rename(
        cons_cluster = cluster,
        cluster = renumbered_cluster
      ) %>%
      mutate(
        cluster = as.character(cluster),
        cons_cluster = as.character(cons_cluster)
      )
  } else {
    # final_clustering <- final_clustering %>%
    #   mutate(cluster = as.character(cluster),
    #          cons_cluster = as.character(cluster))
    #
    # mapping_table <-
    #   data.frame(renumbered_cluster = final_clustering$cluster,
    #              new_cluster = final_clustering$cluster) %>%
    #   distinct()
    #
    mapping_table <-
      final_clustering %>%
      select(new_cluster = cluster) %>%
      mutate(new_cluster = as.character(new_cluster)) %>%
      distinct() %>%
      arrange(as.numeric(new_cluster)) %>%
      rownames_to_column("renumbered_cluster")

    final_clustering <-
      final_clustering %>%
      mutate(cons_cluster = as.character(cluster)) %>%
      select(-cluster) %>%
      left_join(mapping_table, by = c("cons_cluster" = "new_cluster")) %>%
      rename(cluster = renumbered_cluster)
  }

  if (get_membership) {
    # Extract cluster membership matrix
    cons_matrix <-
      cons_clustering$.Data[, ] %>%
      as_tibble(rownames = "gene") %>%
      gather(cluster, membership, -1) %>%
      filter(membership > 0) %>%
      mutate(cluster = as.character(gsub("V", "", cluster))) %>%
      rename(cons_cluster = cluster) %>%
      left_join(mapping_table %>%
        mutate(new_cluster = as.character(new_cluster)) %>%
        rename(
          cons_cluster = new_cluster,
          cluster = renumbered_cluster
        )) %>%
      filter(!is.na(cluster)) # %>%
    #   select(-cluster_cons)

    return(list(
      consensus_clustering = final_clustering,
      membership_matrix = cons_matrix
    ))
  } else {
    return(final_clustering)
  }
}
