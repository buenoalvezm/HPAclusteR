run_gene_umap <- function(adata_t, seed = 42) {
  set.seed(seed)
  if (length(adata_t$obs_names()) > 10000) {
    n_epochs <- 200
  } else {
    n_epochs <- 500
  }
  adata_t$obsm$X_umap_raw <-
    adata_t$uns$neighbors$snn %>%
    RunUMAP(umap.method = "umap-learn", n.epochs = n_epochs) %>%
    {
      .@cell.embeddings
    }

  adata_t$obsm$X_umap_scaled <-
    adata_t$obsm$X_umap_raw %>%
    as_tibble(rownames = "gene") %>%
    gather(UMAP, UMAP_value, UMAP_1, UMAP_2) %>%
    mutate(UMAP_value = scales::rescale(UMAP_value, c(-1, 1))) %>%
    group_by(UMAP) %>%
    mutate(UMAP_value = UMAP_value - mean(range(UMAP_value))) %>%
    ungroup() %>%
    spread(UMAP, UMAP_value) %>%
    column_to_rownames("gene") %>%
    as.matrix()

  adata_t$obsm$X_umap <- adata_t$obsm$X_umap_scaled
}
