cluster_membconf_data <-  adata_t$uns$membership_matrix %>% rename(cluster_membership_confidence =  membership,
                                                                   ensg_id = gene)  %>% mutate(dataset_id = 'tissue') %>% select(dataset_id, everything())

heatmap_palette <-
  colorRampPalette(c("white", "orangered"))(100)


# at this step, make sure that all your dataet are in the same data_loc
heatmap_datasets <- c('region_split' = 'tissue_split', 'region_grouped' = 'tissue_grouped')
heatmap_data_scaled <- heatmap_datasets %>% imap(~ {
  print(paste('reading files from:', data_loc))
  adata_temp <- AnnDatR$new(
    prefix_name = file.path(data_loc, .x),
    layer = 'ntpm',
    var_names = 'ensembl_id'
  )
  adata_temp$X %>% gather(sample, value, -1) %>%
    group_by(!!sym(adata_temp$var_names_col)) %>%
    mutate(
      value_scaled = value / max(value),
      value_zscore = (value - mean(value, na.rm = T)) /
        sd(value, na.rm = T)
    ) %>%
    ungroup() %>% left_join(adata_t$uns$consensus_clustering,
                            by = join_by(!!adata_temp$var_names_col == gene)) %>%
    mutate(heatmap_type = .y) %>% rename('ensg_id' = adata_temp$var_names_col)
}) %>%  bind_rows() %>% mutate(dataset_id = 'tissue') %>%  drop_na(cluster)


```

```{r}
ensure_dir(paste(save_loc, 'svg', 'heatmap', sep = '/'))


heatmap_data_scaled  %>% group_by(heatmap_type, cluster) %>%
  do({
    require(ggrastr)
    require(patchwork)
    g_cluster_data <- .
    g_heatmap_type <- unique(g_cluster_data$heatmap_type)
    plot_savename_zscore <- paste(save_loc, "svg", "heatmap", paste0(g_heatmap_type, "_zscore_",  unique(g_cluster_data$cluster), '.svg'), sep = '/')
    plot_savename_scaled <- paste(save_loc, "svg", "heatmap", paste0(g_heatmap_type, "_scaled_",  unique(g_cluster_data$cluster), '.svg'), sep = '/')
    g_cluster_data_wide <- g_cluster_data %>%
      select(ensg_id, sample, value_zscore) %>%
      spread(sample, value_zscore) %>%
      filter(complete.cases(.)) %>%
      column_to_rownames("ensg_id")
    gene_hclust <-
      g_cluster_data_wide %>%
      dist() %>%
      hclust(method = "ward.D2")
    sample_hclust <-
      g_cluster_data_wide %>%
      t() %>%
      dist() %>%
      hclust(method = "ward.D2")

    plot_data <-
      g_cluster_data %>%
      mutate(ensg_id = factor(ensg_id,
                              rev(gene_hclust$labels[gene_hclust$order])),
             sample = factor(sample,
                             sample_hclust$labels[sample_hclust$order]))

    plot_confidence <-
      cluster_membconf_data %>%
      filter(ensg_id %in% gene_hclust$labels) %>%
      filter(cluster == unique(g_cluster_data$cluster)) %>% #-> addition to filter memberships for other clusters
      mutate(ensg_id = factor(ensg_id,
                              rev(gene_hclust$labels[gene_hclust$order]))) %>%
      ggplot(aes("Confidence", ensg_id, fill = cluster_membership_confidence)) +
      geom_tile_rast() +
      theme_void() +
      scale_fill_gradient(low = "#D1EEEA", high = "#2A5674",
                          limits = c(0.2, 1), # before limit to 0.5
                          breaks = c(0.5, 0.75, 1),
                          name = "Confidence") +
      theme(axis.text.x = element_text(angle = -90, face = "bold"),
            legend.position = "left")
    plot_heat_zscore <-
      plot_data %>%
      ggplot(aes(sample, ensg_id, fill = value_zscore)) +
      geom_tile_rast() +
      scale_fill_gradient2(low = "#4575B4", high = "#D73027", mid = "#FFFFFF",
                           name = "Z-score") +
      theme_void() +
      theme(axis.text.x = element_text(angle = -90,
                                       hjust = 0,
                                       size = 6),
            legend.position = "top")

    plot_heat_scaled <-
      plot_data %>%
      ggplot(aes(sample, ensg_id, fill = value_scaled)) +
      geom_tile_rast() +
      scale_fill_gradient2(low = "#4575B4", high = "#D73027", mid = "#FFFFFF",
                           name = "Relative Expression",
                           breaks = c(0, 0.5, 1),
                           limits = c(0, 1)) +
      theme_void() +
      theme(axis.text.x = element_text(angle = -90,
                                       hjust = 0,
                                       size = 6),
            legend.position = "top")

    column_width <- (1/53)
    plot_n <- length(sample_hclust$labels)
    plot_scaling_factor <- column_width * plot_n

    plot_widths <-
      c(3, plot_n) %>%
      {. * column_width} %>%
      {. / max(.)}

    plot_width <-
      1 + 4 * plot_scaling_factor

    p1 <-
      plot_confidence + plot_heat_zscore +
      plot_layout(widths = plot_widths, nrow = 1)
    p2 <-
      plot_confidence + plot_heat_scaled +
      plot_layout(widths = plot_widths, nrow = 1)

    ggsave(plot_savename_zscore, plot = p1,
           width = plot_width,
           height = 5,
           limitsize = F)
    ggsave(plot_savename_scaled, plot = p2,
           width = plot_width,
           height = 5,
           limitsize = F)

    tibble(.)
  }) %>%
  ungroup() %>%
  collect()
