#' Master annotation pipeline: download, load, and run all enrichments
#'
#' @param AnnDatR AnnDatR object containing clustering results
#' @param dbs Character vector of databases to use for enrichment ("GO", "KEGG", "Others") (default: all)
#' @param universe Character vector of background genes (default: NULL, all genes in clustering data)
#' @param db_loc Directory to store annotation databases (default: "databases")
#' @param hpa_version Version of the Human Protein Atlas to download (default: 24)
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @details In case some of the databases did not download correctly, please rerun the function.
#' If that does not help, please manually download the databases and place them in the specified db_loc directory.
#' Links:
#' - Human Protein Atlas version 24: https://www.proteinatlas.org/download/rna
#' - Trrust:
#' - Panglao:
#'
#' Or just run the annotation only using KEGG and GO by
#' @return A named list with elements: database_enrichment, kegg_enrichment, go_enrichment
#' @export
#'
#' @examples
#' # Run clustering pipeline
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 6.3)
#'
#' # Run annotation pipeline
#' enrichment_results <- hc_annotate(adata_res, dbs = "KEGG")
#' head(enrichment_results)
hc_annotate <- function(
  AnnDatR,
  dbs = c("GO", "KEGG", "Others"),
  universe = NULL,
  db_loc = "databases",
  hpa_version = 24,
  verbose = TRUE
) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "The 'readr' package is required for this function. Please install it using install.packages('readr')."
    )
  }
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
  if ("Others" %in% dbs) {
    # Download annotation databases
    db_files <- get_annot_dbs(db_loc = db_loc, hpa_version = hpa_version)
    if (verbose) {
      message("Annotation databases ready.")
    }

    # Format all annotation databases
    annotation_terms <- build_annotation_terms_tibble(db_files, AnnDatR)
    if (verbose) {
      message("Annotation databases loaded and formatted.")
    }

    # Run enrichment for custom databases (other than KEGG and GO)
    database_enrichment <- run_database_enrichment(
      annotation_terms = annotation_terms,
      AnnDatR = AnnDatR,
      universe = universe,
      verbose = verbose
    )
    if (verbose) {
      message("Custom database enrichment done.")
    }
  } else {
    database_enrichment <- dplyr::tibble()
  }

  if ("KEGG" %in% dbs) {
    # Run KEGG enrichment
    kegg_enrichment <- run_kegg_enrichment(
      AnnDatR = AnnDatR,
      universe = universe,
      verbose = verbose
    )
    if (verbose) {
      message("KEGG enrichment done.")
    }
  } else {
    kegg_enrichment <- dplyr::tibble()
  }

  if ("GO" %in% dbs) {
    # Run GO enrichment with simplification
    go_enrichment <- run_go_enrichment(
      AnnDatR = AnnDatR,
      universe = universe,
      verbose = verbose
    )
    if (verbose) {
      message("GO enrichment (with simplification) done.")
    }
  } else {
    go_enrichment <- dplyr::tibble()
  }

  all_enrichment <- dplyr::bind_rows(
    database_enrichment,
    kegg_enrichment,
    go_enrichment
  ) |>
    dplyr::filter(!!rlang::sym("Adjusted P-value") < 0.05) |>
    dplyr::arrange(!!rlang::sym("Cluster ID"), !!rlang::sym("Adjusted P-value"))

  # Map Ensembl IDs to Gene Symbols for better interpretability
  all_enrichment <- map_ensembl_to_symbol(all_enrichment)

  return(all_enrichment)
}
