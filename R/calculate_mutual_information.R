calculate_mutual_information <-
  function(clustering,
           GO,
           nrep = 10,
           random_simplify = F) {
    overlap_GO <-
      GO %>%
      filter(ensg_id %in% clustering$gene)
    overlap_clustering <-
      clustering %>%
      filter(gene %in% overlap_GO$ensg_id)

    random_clustering_MI <-
      overlap_GO %>%
      group_by(GO_domain) %>%
      do({
        domain_GOterms <- .
        if (random_simplify) {
          domain_GOterms <-
            domain_GOterms %>%
            group_by(ensg_id) %>%
            summarise(GO_accession = sample(GO_accession, 1)) %>%
            ungroup()
        }
        expand_grid(run = 1:nrep) %>%
          group_by(run) %>%
          do({
            overlap_clustering %>%
              mutate(cluster = sample(cluster)) %>%
              inner_join(domain_GOterms,
                by = c("gene" = "ensg_id")
              ) %>%
              summarise(MI = infotheo::mutinformation(cluster, GO_accession))
          })
      }) %>%
      ungroup()

    random_clustering_MI_specs <-
      random_clustering_MI %>%
      group_by(GO_domain) %>%
      summarise(
        mean_MI = mean(MI),
        sd_MI = sd(MI)
      )
    overlap_clustering_MI <-
      overlap_GO %>%
      group_by(GO_domain) %>%
      do({
        domain_GOterms <- .
        if (random_simplify) {
          domain_GOterms <-
            domain_GOterms %>%
            group_by(ensg_id) %>%
            summarise(GO_accession = sample(GO_accession, 1)) %>%
            ungroup()
        }
        overlap_clustering %>%
          inner_join(domain_GOterms,
            by = c("gene" = "ensg_id")
          ) %>%
          summarise(MI = infotheo::mutinformation(cluster, GO_accession))
      }) %>%
      ungroup() %>%
      left_join(random_clustering_MI_specs,
        by = "GO_domain"
      ) %>%
      mutate(z_score = (MI - mean_MI) / sd_MI)
    return(overlap_clustering_MI)
  }

### Mutual Infomration
#### GO terms
GOterms <-
  vroom::vroom("~/Documents/SciLifeDrive/hs_celltype/other_data/hsapiens_ensembl109_GO-terms.txt") %>%
  select(
    ensg_id = 1,
    gene_name = 2,
    GO_accession = 3,
    GO_domain = 4,
    GO_term_name = 5,
    GO_term_evidence = 6
  ) %>%
  filter(ensg_id %in% adata_t$obs_names())


MI_GOterms <- GOterms %>%
  filter(!is.na(GO_accession), !is.na(GO_domain)) %>%
  distinct() %>%
  group_by(GO_accession) %>%
  mutate(n_genes = n_distinct(ensg_id)) %>%
  ungroup() %>%
  filter(n_genes <= 500, n_genes >= 10) %>%
  select(ensg_id, GO_accession, GO_domain)

MI_scores <- adata_t$uns$consensus_clustering %>%
  select(gene, cluster) %>%
  calculate_MI(., MI_GOterms, nrep = 500)

write_tsv(MI_scores, savefile_MIGO)

gene_distances <- adata_t$obs %>%
  left_join(
    adata_t$obsm$X_umap_scaled %>%
      as_tibble() %>% set_names(paste0("UMAP", 1:ncol(.))) %>% mutate(ensembl_id = rownames(adata_t$obsm$X_umap)),
    by = join_by(ensembl_id)
  ) %>%
  select(gene = ensembl_id, cluster, UMAP_1_scaled = UMAP1, UMAP_2_scaled = UMAP2) %>%
  group_by(cluster) %>%
  do({
    dat <- .
    dat %>%
      select(gene, UMAP_1_scaled, UMAP_2_scaled) %>%
      column_to_rownames("gene") %>%
      dist() %>%
      as.matrix() %>%
      as_tibble(rownames = "gene1") %>%
      gather(gene2, dist, -1) %>%
      filter(gene1 < gene2) %>%
      summarise(
        dist = mean(dist),
        n_genes = n_distinct(dat$gene)
      )
  })
write_tsv(gene_distances, intracluster_dist_file)
