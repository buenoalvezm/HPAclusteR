bubble_heatmap_settings <-
  all_res %>%
  filter(!test_type %in% c("KEGG",
                           "GO_BP_original",
                           "GO_CC_original",
                           "GO_MF_original",
                           "GO_BP_simplified",
                           "GO_CC_simplified",
                           "GO_MF_simplified",
                           "reactome")) %>%
  group_by(#dataset_id,
    test_type) %>%
  summarise(bubbles = length(ID),
            clusters = n_distinct(cluster),
            terms = n_distinct(Description)) %>%
  ungroup() %>%
  mutate(i = row_number(),
         bubble_folder = paste(save_loc, "svg", "bubble", sep = '/'),
         bubble_savefile = paste(save_loc, "svg", "bubble", paste0('bubble_', test_type, '.svg'), sep = '/'))


plot_data_all <-
  all_res %>%
  group_by_all() %>%
  mutate(capped_odds_ratio = min(c(odds_ratio,
                                   50))) %>%
  ungroup()

cluster_levels <-
  clustering_data %>%
  select(#dataset_id,
    cluster) %>%
  distinct() %>%
  arrange(#dataset_id,
    as.numeric(cluster))


plot_size_range <- c(1, 4)

plot_legend_settings <-
  tibble(odds_ratio_str = c(1, 10, 20, 30, 40, ">50"),
         y = 1:6,
         odds_ratio = c(1, 10, 20, 30, 40, 50))

plot_legend <-
  plot_legend_settings %>%
  ggplot(aes(1, y,
             size = odds_ratio,
             fill = odds_ratio,
             label = odds_ratio_str)) +
  geom_point(shape = 21,
             show.legend = F) +
  geom_text(size = 4,
            hjust = 0,
            nudge_x = 0.1) +
  annotate("text", x = 1, y = 7, label = "Odds ratio",
           fontface = "bold", hjust = 0) +
  scale_x_continuous(limits = c(0, 2)) +
  scale_y_continuous(limits = c(0, 7)) +
  scale_fill_gradient(name = "Odds ratio",
                      low = "white", high = "orangered",
                      breaks = plot_legend_settings$odds_ratio,
                      limits = range(plot_legend_settings$odds_ratio)) +
  scale_size_continuous(name = "Odds ratio",
                        range = plot_size_range,
                        breaks = plot_legend_settings$odds_ratio,
                        limits = range(plot_legend_settings$odds_ratio)) +
  theme_void()

for(folder_ in unique(bubble_heatmap_settings$bubble_folder)) {
  ggsave(paste0(folder_, "/legend.svg"),
         plot = plot_legend,
         width = 2,
         height = 2)
}

plots <-
  pbapply::pblapply(1:nrow(bubble_heatmap_settings),
                    function(i) {
                      i <<- i
                      #dataset_id_ <- bubble_heatmap_settings$dataset_id[i]
                      test_type_ <- bubble_heatmap_settings$test_type[i]
                      filename <- bubble_heatmap_settings$bubble_savefile[i]

                      plot_data <-
                        plot_data_all %>%
                        filter(#dataset_id == dataset_id_,
                          test_type == test_type_) %>%
                        mutate(cluster = paste0("cluster_", cluster, "_"))

                      plot_levels <-
                        cluster_levels %>%
                        mutate(cluster = paste0("cluster_", cluster, "_")) %>%
                        # filter(dataset_id == dataset_id_) %>%
                        pull(cluster)

                      plot_clustering <-
                        plot_data %>%
                        select(cluster, ID, odds_ratio) %>%
                        spread(ID, odds_ratio, fill = 0) %>%
                        column_to_rownames("cluster") %>%
                        list(tissue = t(.),
                             cluster = .) %>%
                        map(. %>%
                              dist(method = "binary") %>%
                              hclust(method = "ward.D2"))

                      plot_order <-
                        list(cluster = unique(c(with(plot_clustering$cluster,
                                                     labels[order]),
                                                plot_levels)),
                             ID = with(plot_clustering$tissue,
                                       labels[order]))

                      plot <-
                        plot_data %>%
                        ungroup() %>%
                        mutate(cluster = factor(cluster,
                                                plot_order$cluster),
                               ID = factor(ID,
                                           plot_order$ID)) %>%
                        ggplot(aes(ID, cluster,
                                   size = capped_odds_ratio,
                                   fill = capped_odds_ratio)) +
                        geom_point(data = expand_grid(cluster = plot_order$cluster,
                                                      ID = plot_order$ID) %>%
                                     mutate(cluster = factor(cluster, unique(cluster)),
                                            ID = factor(ID, unique(ID))),
                                   fill = NA,
                                   size = NA) +
                        geom_point(shape = 21,
                                   show.legend = F) +
                        theme_bw() +
                        theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
                              legend.position = "right") +
                        # coord_fixed() +
                        scale_fill_gradient(name = "Odds ratio",
                                            low = "white", high = "orangered",
                                            breaks = plot_legend_settings$odds_ratio,
                                            limits = range(plot_legend_settings$odds_ratio)) +
                        scale_size_continuous(name = "Odds ratio",
                                              range = plot_size_range,
                                              breaks = plot_legend_settings$odds_ratio,
                                              limits = range(plot_legend_settings$odds_ratio)) +
                        theme(axis.title = element_blank(),
                              plot.margin = unit(c(0,0,0,200), "pt")) +
                        coord_fixed()

                      column_width <- (1/30)#(1/53)
                      plot_n <- length(plot_order$ID)
                      plot_scaling_factor_row <- column_width * plot_n

                      plot_n_col <- length(plot_order$cluster)
                      plot_scaling_factor_col <- column_width * plot_n_col

                      #
                      # plot_widths <-
                      #   c(3, length(plot_n)) %>%
                      #   {. * column_width} %>%
                      #   {. / max(.)}

                      plot_width <-
                        4 + 4 * plot_scaling_factor_row
                      plot_height <-
                        2 + 4 * plot_scaling_factor_col

                      ggsave(filename,
                             plot = plot,
                             width = plot_width, height = plot_height)

                    })

