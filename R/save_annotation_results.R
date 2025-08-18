all_res <-
  list(
    KEGG_enrichment_res %>%
      mutate(test_type = "KEGG") %>%
      select(-Cluster),
    database_enrichment_res %>%
      rename(test_type = id) %>%
      left_join(database_terms %>%
        select(
          test_type = id,
          ID = term_id,
          Description2 = term
        ) %>%
        distinct()) %>%
      mutate(Description = ifelse(test_type != "KEGG",
        Description2,
        Description
      )) %>%
      select(-Description2),
    database_enrichment_res_simplified_GO %>%
      rename(test_type = id) %>%
      left_join(database_terms %>%
        select(
          ID = term_id,
          Description2 = term
        ) %>%
        distinct()) %>%
      mutate(Description = ifelse(test_type != "KEGG",
        Description2,
        Description
      )) %>%
      select(-Description2)
  ) %>%
  bind_rows() %>%
  mutate(
    GeneFrac = as.numeric(gsub("\\/.*", "", GeneRatio)) /
      as.numeric(gsub(".*\\/", "", GeneRatio)),
    BgFrac = as.numeric(gsub("\\/.*", "", BgRatio)) /
      as.numeric(gsub(".*\\/", "", BgRatio)),
    odds_ratio = GeneFrac / BgFrac
  ) %>%
  left_join(enrichment_settings_versions %>%
    select(
      test_type = id,
      test_description = description,
      database_version = version
    )) %>%
  select( # dataset_id,
    test_type, test_description, database_version,
    cluster, ID, Description,
    GeneRatio, GeneFrac, BgRatio, BgFrac, odds_ratio, everything()
  )


write_tsv(all_res, enrichment_savefile)
