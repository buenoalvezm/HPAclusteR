#' Run enrichment analysis for all annotation databases
#'
#' @param annotation_terms Tibble from build_annotation_terms_tibble() (db_id, ensg_id, term, term_id)
#' @param AnnDatR AnnDatR object containing clustering results
#' @param universe Character vector of background genes (default: NULL, all genes in clustering data)
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#'
#' @return Tibble with enrichment results for all clusters and databases
#' @keywords internal
run_database_enrichment <- function(
  annotation_terms,
  AnnDatR,
  universe = NULL,
  verbose = TRUE
) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop(
      "The 'clusterProfiler' package is required for this function. Please install it using BiocManager::install('clusterProfiler')."
    )
  }
  if (is.null(AnnDatR[["uns"]][["consensus_clustering"]])) {
    stop(
      "AnnDatR$uns$consensus_clustering not found. Call `hc_cluster_consensus()` before `hc_annotate()`."
    )
  }
  clustering_data <- AnnDatR[["uns"]][["consensus_clustering"]]
  if (is.null(universe)) {
    universe <- unique(clustering_data[["gene"]])
  }
  db_ids <- unique(annotation_terms[["db_id"]])
  res_list <- lapply(db_ids, function(db_id) {
    db_terms <- annotation_terms[annotation_terms[["db_id"]] == db_id, ]
    term2gene <- db_terms[, c("term_id", "ensg_id")]
    clusters <- unique(clustering_data[["cluster"]])
    cluster_res <- lapply(clusters, function(clust) {
      genes_in_cluster <- clustering_data[["gene"]][
        clustering_data[["cluster"]] == clust
      ]
      enr <- tryCatch(
        suppressMessages(
          clusterProfiler::enricher(
            gene = genes_in_cluster,
            universe = universe,
            TERM2GENE = term2gene
          )
        ),
        error = function(e) NULL
      )
      if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
        if (verbose) {
          message(sprintf(
            "No enrichment result for db: %s, cluster: %s",
            db_id,
            clust
          ))
        }
        return(NULL)
      }
      tibble::as_tibble(enr) |>
        dplyr::mutate(
          cluster = !!rlang::sym("clust"),
          db_id = !!rlang::sym("db_id")
        )
    })
    dplyr::bind_rows(cluster_res)
  })

  enrich_res <- dplyr::bind_rows(res_list) |>
    dplyr::mutate(
      Database = dplyr::case_when(
        db_id == "specificity_tissue" ~ "Specificity classification Tissue",
        db_id == "specificity_brain" ~ "Specificity classification Brain",
        db_id == "specificity_celline" ~ "Specificity classification Cell line",
        db_id ==
          "specificity_singlecell" ~ "Specificity classification Single cell type",
        db_id == "specificity_blood" ~ "Specificity classification Blood",
        db_id == "secretome_location" ~ "Secretome location",
        db_id == "subcellular_location" ~ "Subcellular location",
        db_id == "trrust" ~ "TRRUST TFs",
        db_id == "protein_class" ~ "HPA protein class",
        db_id ==
          "panglao_cellmarkers" ~ "Cell type markers from the Panglao database",
        TRUE ~ db_id
      )
    ) |>
    dplyr::select(
      `Cluster ID` = !!rlang::sym("cluster"),
      !!rlang::sym("Database"),
      `Term ID` = !!rlang::sym("ID"),
      Term = !!rlang::sym("Description"),
      !!rlang::sym("GeneRatio"),
      !!rlang::sym("BgRatio"),
      `P-value` = !!rlang::sym("pvalue"),
      `Adjusted P-value` = !!rlang::sym("p.adjust"),
      `Gene IDs` = !!rlang::sym("geneID")
    )

  return(enrich_res)
}

#' Map Entrez IDs in a column back to Ensembl IDs
#'
#' @param df Data frame with a column named 'Gene IDs' (slash-separated Entrez IDs)
#' @param geneid_col Name of the column with Entrez IDs (default: "Gene IDs")
#' @return Data frame with an added column of Ensembl IDs (slash-separated)
map_entrez_to_ensembl <- function(
  df,
  geneid_col = "Gene IDs"
) {
  # Get all unique Entrez IDs
  all_entrez <- unique(unlist(strsplit(df[[geneid_col]], "/")))
  # Map Entrez to Ensembl
  mapping <- clusterProfiler::bitr(
    all_entrez,
    fromType = "ENTREZID",
    toType = "ENSEMBL",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  # Helper to map a slash-separated string
  map_ids <- function(entrez_str) {
    ids <- unlist(strsplit(entrez_str, "/"))
    mapped <- mapping[["ENSEMBL"]][match(ids, mapping[["ENTREZID"]])]
    mapped <- mapped[!is.na(mapped)]
    paste(mapped, collapse = "/")
  }
  df[[geneid_col]] <- vapply(df[[geneid_col]], map_ids, character(1))

  return(df)
}

#' Map Ensembl IDs in a column to gene symbols and add as a new column
#'
#' @param df Data frame with a column of Ensembl IDs (slash-separated, default: "Gene IDs")
#' @param geneid_col Name of the column with Ensembl IDs (default: "Gene IDs")
#' @param new_col Name for the new gene symbol column (default: "Gene names")
#' @return Data frame with an added column of gene symbols (slash-separated)
map_ensembl_to_symbol <- function(
  df,
  geneid_col = "Gene IDs",
  new_col = "Gene names"
) {
  # Get all unique Ensembl IDs
  all_ensg <- unique(unlist(strsplit(df[[geneid_col]], "/")))
  # Map Ensembl to gene symbol
  mapping <- clusterProfiler::bitr(
    all_ensg,
    fromType = "ENSEMBL",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  # Helper to map a slash-separated string
  map_ids <- function(ensg_str) {
    ids <- unlist(strsplit(ensg_str, "/"))
    mapped <- mapping[["SYMBOL"]][match(ids, mapping[["ENSEMBL"]])]
    mapped <- mapped[!is.na(mapped)]
    paste(mapped, collapse = "/")
  }
  df[[new_col]] <- vapply(df[[geneid_col]], map_ids, character(1))
  return(df)
}

#' Run KEGG enrichment analysis (with Ensembl-to-Entrez conversion)
#'
#' @param AnnDatR AnnDatR object containing clustering results
#' @param universe Character vector of background genes (default: NULL, all genes in clustering data)
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#'
#' @return Tibble with enrichment results for all clusters
#' @keywords internal
run_kegg_enrichment <- function(
  AnnDatR,
  universe = NULL,
  verbose = TRUE
) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop(
      "The 'clusterProfiler' package is required for this function. Please install it using BiocManager::install('clusterProfiler')."
    )
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop(
      "The 'org.Hs.eg.db' package is required for KEGG enrichment. Please install it using BiocManager::install('org.Hs.eg.db')."
    )
  }
  if (is.null(AnnDatR[["uns"]][["consensus_clustering"]])) {
    stop(
      "AnnDatR$uns$consensus_clustering not found. Call `hc_cluster_consensus()` before `hc_annotate()`."
    )
  }
  clustering_data <- AnnDatR[["uns"]][["consensus_clustering"]]
  gene_map <- clusterProfiler::bitr(
    unique(clustering_data[["gene"]]),
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  clustering_data_entrez <- clustering_data |>
    dplyr::left_join(gene_map, by = c("gene" = "ENSEMBL")) |>
    dplyr::filter(!is.na(!!rlang::sym("ENTREZID")))
  if (is.null(universe)) {
    universe <- clustering_data_entrez[["ENTREZID"]]
  }
  clusters <- unique(clustering_data_entrez[["cluster"]])
  res_list <- lapply(clusters, function(clust) {
    genes_in_cluster <- clustering_data_entrez[["ENTREZID"]][
      clustering_data_entrez[["cluster"]] == clust
    ]
    enr <- tryCatch(
      suppressMessages(
        clusterProfiler::enrichKEGG(
          pvalueCutoff = 0.05,
          gene = genes_in_cluster,
          universe = universe,
          organism = "hsa"
        )
      ),
      error = function(e) NULL
    )
    if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
      if (verbose) {
        message(sprintf(
          "No KEGG enrichment result for cluster: %s",
          clust
        ))
      }
      return(NULL)
    }
    tibble::as_tibble(enr) |>
      dplyr::mutate(
        cluster = clust,
        db_id = "KEGG pathways"
      )
  })

  enrich_res <- dplyr::bind_rows(res_list) |>
    dplyr::select(
      `Cluster ID` = !!rlang::sym("cluster"),
      Database = !!rlang::sym("db_id"),
      `Term ID` = !!rlang::sym("ID"),
      Term = !!rlang::sym("Description"),
      !!rlang::sym("GeneRatio"),
      !!rlang::sym("BgRatio"),
      `P-value` = !!rlang::sym("pvalue"),
      `Adjusted P-value` = !!rlang::sym("p.adjust"),
      `Gene IDs` = !!rlang::sym("geneID")
    )

  enrich_res <- map_entrez_to_ensembl(enrich_res)

  return(enrich_res)
}

#' Run GO enrichment analysis (BP, MF, CC) for all clusters using clusterProfiler::enrichGO
#'
#' @param AnnDatR AnnDatR object containing clustering results
#' @param universe Character vector of background genes (default: NULL, all genes in clustering data)
#' @param verbose Logical indicating whether to print progress messages (default is TRUE).
#'
#' @return Tibble with enrichment results for all clusters and all GO ontologies, with Ensembl gene IDs
#' @keywords internal
run_go_enrichment <- function(
  AnnDatR,
  universe = NULL,
  verbose = TRUE
) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop(
      "The 'clusterProfiler' package is required for this function. Please install it using BiocManager::install('clusterProfiler')."
    )
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop(
      "The 'org.Hs.eg.db' package is required for GO enrichment. Please install it using BiocManager::install('org.Hs.eg.db')."
    )
  }
  if (is.null(AnnDatR[["uns"]][["consensus_clustering"]])) {
    stop(
      "AnnDatR$uns$consensus_clustering not found. Call `hc_cluster_consensus()` before `hc_annotate()`."
    )
  }
  clustering_data <- AnnDatR[["uns"]][["consensus_clustering"]]
  gene_map <- clusterProfiler::bitr(
    unique(clustering_data[["gene"]]),
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  clustering_data_entrez <- clustering_data |>
    dplyr::left_join(gene_map, by = c("gene" = "ENSEMBL")) |>
    dplyr::filter(!is.na(!!rlang::sym("ENTREZID")))
  if (is.null(universe)) {
    universe <- clustering_data_entrez[["ENTREZID"]]
  }
  clusters <- unique(clustering_data_entrez[["cluster"]])
  ontologies <- c("BP", "MF", "CC")
  ont_full <- c(
    BP = "Biological Process",
    MF = "Molecular Function",
    CC = "Cellular Component"
  )
  res_list <- lapply(clusters, function(clust) {
    genes_in_cluster <- clustering_data_entrez[["ENTREZID"]][
      clustering_data_entrez[["cluster"]] == clust
    ]
    ont_res <- lapply(ontologies, function(ont) {
      enr <- tryCatch(
        suppressMessages(
          clusterProfiler::enrichGO(
            gene = genes_in_cluster,
            universe = universe,
            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
            keyType = "ENTREZID",
            ont = ont,
            pvalueCutoff = 0.05
          )
        ),
        error = function(e) NULL
      )
      if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
        if (verbose) {
          message(sprintf(
            "No GO enrichment result for cluster: %s, ontology: %s",
            clust,
            ont
          ))
        }
        return(NULL)
      }
      # Normal results
      normal <- tibble::as_tibble(enr) |>
        dplyr::mutate(
          cluster = clust,
          db_id = paste0("GO analysis ", ont_full[[ont]])
        )
      # Simplified results
      enr_simple <- tryCatch(
        suppressMessages(
          clusterProfiler::simplify(
            enr,
            cutoff = 0.7,
            by = "p.adjust",
            select_fun = min
          )
        ),
        error = function(e) NULL
      )
      if (!is.null(enr_simple) && nrow(as.data.frame(enr_simple)) > 0) {
        simplified <- tibble::as_tibble(enr_simple) |>
          dplyr::mutate(
            cluster = clust,
            db_id = paste0(
              "GO analysis ",
              ont_full[[ont]],
              " (Simplified terms)"
            )
          )
        dplyr::bind_rows(normal, simplified)
      } else {
        normal
      }
    })
    dplyr::bind_rows(ont_res)
  })

  enrich_res <- dplyr::bind_rows(res_list) |>
    dplyr::select(
      `Cluster ID` = !!rlang::sym("cluster"),
      Database = !!rlang::sym("db_id"),
      `Term ID` = !!rlang::sym("ID"),
      Term = !!rlang::sym("Description"),
      !!rlang::sym("GeneRatio"),
      !!rlang::sym("BgRatio"),
      `P-value` = !!rlang::sym("pvalue"),
      `Adjusted P-value` = !!rlang::sym("p.adjust"),
      `Gene IDs` = !!rlang::sym("geneID")
    )
  enrich_res <- map_entrez_to_ensembl(enrich_res)

  return(enrich_res)
}
