hpa_gene_clustering <- function(adata_t, resolution, n_neighbors = 20, pruning = 1/15, clustering_method = 'louvain', num_seeds = 100){
  set_PCA(adata_t, scale_by = 'sample', transform = 'none', nPcs = 40) ## note, they are not log_transforming
  # adata_t$uns$pca_stats <- summary(adata_t$uns$pca) %>%
  #   t() %>%
  #   as_tibble(rownames = "PC") %>%
  #   mutate(PC = as.numeric(gsub("PC", "", PC))) %>%
  #   dplyr::rename(R2cum = 3)

  get_distance(adata_t, method = 'spearman', n_comp = kaisers_PCA_rule(adata_t$uns$pca), id = 'zscore_spearman')

  adata_t$uns$neighbors <- adata_t$uns$distance$zscore_spearman %>%
    as.matrix() %>%
    as.dist() %>%
    Seurat::FindNeighbors(
      k.param = n_neighbors,
      compute.SNN = TRUE,
      prune.SNN = pruning,
      nn.method = "annoy", #Distance metric for annoy. Options include: euclidean, cosine, manhattan, and hamming
      annoy.metric = "euclidean",
      nn.eps = 0,
      verbose = TRUE,
      force.recalc = FALSE)

  run_consensus_clustering(adata_t, resolution = resolution, num_seeds = num_seeds, clustering_method = clustering_method)
  run_gene_umap(adata_t)

  adata_t$obs <- adata_t$obs %>% left_join(
    adata_t$obsm$X_umap_scaled  %>%
      as_tibble() %>%  set_names(paste0("UMAP", 1:ncol(.)))  %>% mutate(ensembl_id = rownames(adata_t$obsm$X_umap)),
    by = join_by(ensembl_id)
  )

  adata_t$uns$UMAP_hulls <- generate_cluster_hulls(V1 = adata_t$obs[['UMAP1']],
                                                   V2 = adata_t$obs[['UMAP2']],
                                                   element_id = adata_t$obs [['ensembl_id']],
                                                   cluster_membership =adata_t$obs [['cluster']] ,
                                                   n = 1000,
                                                   poly_concavity = 2,
                                                   poly_smoothing = 1.25,
                                                   relative_bandwidth = 1/150)
  adata_t$obs <- adata_t$obs %>% select(-UMAP1, -UMAP2)


}
