clustering_data <- adata_t$uns$consensus_clustering
gene_info <- read_tsv("~/Documents/SciLifeDrive/hs_celltype/other_data/dbs/geneinfo_109.tsv")
entrez_dict <-
  read_tsv("~/Documents/SciLifeDrive/hs_celltype/other_data/ensembl109_entrez.txt") %>%
  select(
    ensg_id = 1,
    entrez = 2
  ) %>%
  filter(!is.na(entrez))

GO_formatting_function <-
  . %>%
  select(
    ensg_id = `Gene stable ID`,
    term = `GO term name`,
    term_id = `GO term accession`
  ) %>%
  filter(!is.na(term)) %>%
  distinct() %>%
  filter(ensg_id %in% clustering_data$gene) %>%
  group_by(term_id, term) %>%
  filter(n_distinct(ensg_id) >= 8) %>%
  ungroup()
HPAspecificity_formatting_function <-
  . %>%
  select(ensg_id,
    term = enhanced_tissues
  ) %>%
  mutate(term = gsub(", ", ";", term)) %>%
  separate_rows(term, sep = ",") %>%
  mutate(term = gsub(";", ", ", term)) %>%
  filter(!is.na(term)) %>%
  filter(term != "NULL") %>%
  mutate(
    term = trimws(term),
    term_id = term
  )


database_formatting_functions <-
  list(
    secretome_location = . %>%
      select(
        ensg_id = Ensembl,
        term = `Secretome location`
      ) %>%
      filter(!is.na(term)) %>%
      mutate(term_id = term),
    specificity_blood = HPAspecificity_formatting_function,
    specificity_brain = HPAspecificity_formatting_function,
    specificity_tissue = HPAspecificity_formatting_function,
    specificity_celline = HPAspecificity_formatting_function,
    specificity_singlecell = HPAspecificity_formatting_function,
    subcellular_location = . %>%
      filter(Gene %in% clustering_data$gene) %>%
      select(
        ensg_id = Gene,
        Enhanced,
        Supported,
        Approved
      ) %>%
      gather(reliability, location, -1) %>%
      filter(!is.na(location)) %>%
      separate_rows(location, sep = ";") %>%
      select(ensg_id, term = location) %>%
      mutate(term_id = term),
    panglao_cellmarkers = . %>%
      filter(
        grepl("Hs", species),
        `gene type` %in% c(
          "protein-coding gene",
          "protein coding gene"
        )
      ) %>%
      select(
        gene_name = 2,
        term = 3
      ) %>%
      distinct() %>%
      mutate(term_id = term),
    protein_class = . %>% select(
      ensg_id = Ensembl,
      protein_class = `Protein class`
    ) %>%
      separate_rows(protein_class, sep = ", ") %>%
      rename(term = protein_class) %>%
      mutate(term_id = term),
    reactome = . %>%
      filter(species == "Homo sapiens") %>%
      select(ensg_id, term = description, term_id = id) %>%
      distinct() %>%
      filter(ensg_id %in% clustering_data$gene),
    trrust = . %>%
      select(gene_name, term) %>%
      mutate(term_id = term),
    GO_BP_original = . %>%
      filter(`GO domain` == "biological_process") %>%
      GO_formatting_function(),
    GO_MF_original = . %>%
      filter(`GO domain` == "molecular_function") %>%
      GO_formatting_function(),
    GO_CC_original = . %>%
      filter(`GO domain` == "cellular_component") %>%
      GO_formatting_function()
  )



gene_formatting_functions <-
  list(`gene names` = . %>%
    left_join(
      gene_info %>%
        select(ensg_id, gene_name),
      by = c("gene_name")
    ) %>%
    select(ensg_id, term, term_id))





get_db <-
  function(enrichment_settings, db_id, db_loc) {
    db_file <-
      enrichment_settings %>%
      filter(id == db_id) %>%
      pull(file)

    raw_db <-
      read_tsv(paste(db_loc, db_file, sep = "/"))

    db_format_function <-
      database_formatting_functions[[db_id]]

    raw_db %>%
      db_format_function()
  }


db_loc <- "../other_data/dbs"
dataset_id_name <- "tissue"
list.files(db_loc)
enrichment_settings_file <- paste(db_loc, "20250414_enrichment_settings.csv", sep = "/")
enrichment_settings <-
  read_csv(enrichment_settings_file)


enrichment_settings_versions <-
  enrichment_settings %>%
  mutate(version = case_when(
    is.na(version) &
      id == "KEGG" ~ gsub("-", "", Sys.Date()),
    !is.na(version) ~ version
  ))



database_terms <-
  enrichment_settings %>%
  filter(!is.na(file)) %>%
  group_by(id) %>%
  do({
    g_data <<- .

    out_db <-
      get_db(enrichment_settings = enrichment_settings, db_id = g_data$id, db_loc = db_loc)

    if (g_data$gene_format != "ensembl") {
      gene_format_function <- gene_formatting_functions[[g_data$gene_format]]

      out_db <-
        out_db %>%
        gene_format_function()
    }

    stopifnot(all(names(out_db) == c("ensg_id", "term", "term_id")))

    out_db
  }) %>%
  ungroup()


database_terms %>%
  filter(id %in% c(
    "GO_BP_original",
    "GO_CC_original",
    "GO_MF_original"
  )) %>%
  group_by(id, term_id, term) %>%
  summarise(n = n_distinct(ensg_id)) %>%
  # filter(grepl("GO", id)) %>%
  mutate(type = case_when(
    n < 8 ~ as.character(n),
    T ~ "8+"
  )) %>%
  group_by(id, type) %>%
  count() %>%
  ggplot(aes(id, n, label = n, fill = type)) +
  geom_col(color = "black") +
  geom_text(position = position_stack(vjust = 0.5)) +
  facet_wrap(~id, nrow = 1, scales = "free") +
  scale_fill_manual(values = colorRampPalette(c("orangered", "gray"))(8)) +
  theme_bw()

write_tsv(database_terms, paste0("gene_clustering/tissue_HPA25/enrichment/", str_remove(str_remove(enrichment_settings_file, ".*\\/"), " .*")))
