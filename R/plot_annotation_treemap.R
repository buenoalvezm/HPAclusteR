spread_colors <-
  function(color, n, colorrange = 0.5) {
    if (n == 1) {
      return(color)
    } else {
      colorRamp(c("white", color, "black"))(seq(0.5 - colorrange / 2,
        0.5 + colorrange / 2,
        length.out = n
      )) %>%
        as_tibble() %$%
        rgb(V1, V2, V3, maxColorValue = 255)
    }
  }
ontology_pal <-
  c(
    "MF" = "#FF6F00",
    "CC" = "#C71B00",
    "BP" = "#018EA0"
  )



# Generate plot legend:


plot_legend <-
  enframe(
    ontology_pal,
    "Ontology",
    "Color"
  ) %>%
  ggplot(aes(1, 1, fill = Ontology)) +
  geom_tile() +
  scale_fill_manual(values = ontology_pal)
plot_legend <-
  ggpubr::get_legend(plot_legend) %>%
  ggpubr::as_ggplot()

file_structure <-
  create_folder_structure(
    dataset_id = dataset_metadata$dataset_id,
    dataset_run_id = dataset_metadata$version,
    main_folder = "results/Clustering_results",
    folders = c(
      "data" = "data",
      "svg" = "svg",
      "bubbleheatmap" = "svg/bubble",
      "treemap" = "svg/treemap",
      "heatmap" = "svg/heatmap",
      "enrichment" = "enrichment",
      "clustering" = "clustering"
    )
  )

ggsave(paste(save_loc, "svg", "treemaps", "treemap_legend.svg", sep = "/"),
  plot = plot_legend,
  width = 1, height = 1
)


# Generate plot data for treemaps

# plot_data %>%
#   group_by(dataset_id, cluster) %>%
#   do({
#     g_data <- .
#
#     treemap_savefile <-
#       paste0(file_structure[[unique(g_data$dataset_id)]]$treemap,
#              "/",
#              "treemap_",
#              unique(g_data$cluster),
#              ".svg")
#
#
#
#   })

plot_settings <-
  simplified_GO_res %>%
  select(cluster) %>%
  mutate(cluster = as.numeric(cluster)) %>%
  arrange(cluster) %>%
  distinct() %>%
  ungroup() %>%
  mutate(
    i = row_number(),
    treemap_folder = paste(save_loc, "svg", "treemaps", sep = "/"),
    treemap_savefile = paste0(
      treemap_folder,
      "/",
      "treemap_",
      cluster,
      ".svg"
    )
  )
library(magrittr)
library(treemapify)
plots <-
  pbapply::pblapply(
    plot_settings$i,
    function(i) {
      simplified_GO_res %>%
        filter( # dataset_id == plot_settings$dataset_id[i],
          cluster == plot_settings$cluster[i]
        ) %>%
        left_join(
          enframe(
            ontology_pal,
            "ontology",
            "color"
          ),
          by = "ontology"
        ) %>%
        group_by(parent) %>%
        mutate(
          nn = n_distinct(go),
          color2 = spread_colors(unique(color), n_distinct(go), colorrange = 0.2)
        ) %>%
        ggplot(aes(area = score, subgroup = parentTerm)) +
        geom_treemap(aes(fill = color2),
          color = "black",
          show.legend = T
        ) +
        geom_treemap_subgroup_border(color = "black") +
        geom_treemap_text(aes(label = term),
          colour = "black",
          place = "centre",
          alpha = 0.4,
          grow = TRUE
        ) +
        geom_treemap_subgroup_text(
          place = "centre",
          grow = T,
          reflow = T,
          alpha = 1,
          colour = "white",
          fontface = "bold",
          min.size = 0
        ) +
        scale_fill_identity() +
        theme_void()
    }
  )


for (i in plot_settings$i) {
  setting <-
    plot_settings[i, ]
  ggsave(setting$treemap_savefile,
    plots[[i]],
    width = 6, height = 4
  )
}
