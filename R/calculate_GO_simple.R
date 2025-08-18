library(rrvgo)
# Generate similarity matrix for all significant GO terms:
GO_sim <-
  pbapply::pblapply(
    c(
      "BP" = "GO_BP_original",
      "CC" = "GO_CC_original",
      "MF" = "GO_MF_original"
    ),
    function(ont_) {
      GO_terms <-
        database_enrichment_res %>%
        filter(id == ont_) %>%
        pull(ID) %>%
        unique()

      calculateSimMatrix(GO_terms,
        orgdb = "org.Hs.eg.db",
        keytype = "ENSEMBL",
        ont = str_extract(ont_, "BP|CC|MF"),
        method = "Wang"
      )
    }
  )


GO_sim_terms <-
  GO_sim %>%
  map(. %>%
    colnames() %>%
    enframe()) %>%
  bind_rows() %>%
  pull(value)

worker_cluster <- new_cluster(n = 6)
cluster_library(worker_cluster, c(
  "dplyr",
  "tidyverse"
))
cluster_copy(worker_cluster, c(
  "GO_sim",
  "reduceSimMatrix"
))


simplified_GO_res <-
  database_enrichment_res %>%
  filter(id %in% c(
    "GO_BP_original",
    "GO_CC_original",
    "GO_MF_original"
  )) %>%
  filter(ID %in% GO_sim_terms) %>%
  mutate(
    GeneFrac = as.numeric(gsub("\\/.*", "", GeneRatio)) /
      as.numeric(gsub(".*\\/", "", GeneRatio)),
    ontology = str_extract(id, "BP|CC|MF")
  ) %>%
  group_by(
    ontology, id,
    # dataset_id,
    cluster
  ) %>%
  filter(n_distinct(ID) >= 2) %>%
  group_by(
    ontology, id,
    # dataset_id,
    cluster
  ) %>%
  partition(worker_cluster) %>%
  do({
    g_data <<- .

    g_terms <-
      g_data$ID[which(g_data$ID %in% colnames(GO_sim[[unique(g_data$ontology)]]))]

    sim_mat <-
      GO_sim[[unique(g_data$ontology)]][
        g_terms,
        g_terms
      ]

    sim_mat %>%
      reduceSimMatrix(setNames(-log10(g_data$p.adjust), g_data$ID),
        threshold = 0.7,
        orgdb = "org.Hs.eg.db"
      ) %>%
      as_tibble() %>%
      rename(term_cluster = cluster)
  }) %>%
  ungroup() %>%
  collect()

rm(worker_cluster)


database_enrichment_res_simplified_GO <-
  simplified_GO_res %>%
  select(any_of(names(database_enrichment_res)),
    ID = parent
  ) %>%
  distinct() %>%
  left_join(database_enrichment_res) %>%
  bind_rows(database_enrichment_res %>%
    filter(id %in% c(
      "GO_BP_original",
      "GO_CC_original",
      "GO_MF_original"
    )) %>%
    group_by(
      id,
      # dataset_id,
      cluster
    ) %>%
    filter(!(n_distinct(ID) >= 2))) %>%
  mutate(id = gsub("original", "simplified", id))
