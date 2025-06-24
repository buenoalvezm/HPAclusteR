num_seeds <- 100
resolution <- 6.3
ARI_scores <- combn(1:num_seeds, m = 2) %>%
  t() %>%
  as_tibble() %>%
  set_names(c("seed1", "seed2")) %>%
  group_by_all() %>%
  mutate(ARI = aricode::ARI(
    adata_t$uns$cluster_data %>% select(gene, paste("seed", seed1, sep = "_")) %>%  deframe(),
    adata_t$uns$cluster_data %>% select(gene, paste("seed", seed2, sep = "_")) %>%  deframe()
  ))
write_tsv(ARI_scores, savefile_ARI)
