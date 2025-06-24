library(multidplyr)
library(clusterProfiler)
worker_cluster <- new_cluster(n = 5)
cluster_library(worker_cluster, c("dplyr",
                                  "tidyverse"))

cluster_copy(worker_cluster, c("compareCluster"))

clustering_data_entrez <-
  clustering_data %>%
  left_join(entrez_dict,
            by = c("gene" = "ensg_id"))

KEGG_enrichment_res <-
  clustering_data_entrez %>%
  partition(worker_cluster) %>%
  do({
    dataset_data <- .
    compareCluster(entrez ~ cluster,
                   fun = "enrichKEGG",
                   data = dataset_data,
                   universe = unique(dataset_data$entrez)) %>%
      as_tibble()
  })  %>% collect()

rm(worker_cluster)


worker_cluster <- new_cluster(n = 5)
cluster_library(worker_cluster, c("dplyr",
                                  "tidyverse"))
cluster_copy(worker_cluster, c("enricher"))

database_enrichment_res <-
  database_terms %>%
  group_by(id) %>%
  do({
    # For each database
    database <- .

    term2gene <-
      database %>%
      select(term_id, ensg_id)

    cluster_copy(worker_cluster, c("term2gene"))

    # For
    clustering_data %>%
      do({
        g_clustering_data <- .

        g_clustering_data %>%
          group_by(cluster) %>%
          partition(worker_cluster) %>%
          do({
            g_data <- .

            pull(g_data, gene) %>%
              enricher(maxGSSize = Inf,
                       universe = unique(g_clustering_data$gene),
                       TERM2GENE = term2gene) %>%
              as_tibble()
          }) %>%
          ungroup() %>%
          collect()

      })
  }) %>%
  ungroup()

rm(worker_cluster)
