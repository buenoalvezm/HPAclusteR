#' Master annotation pipeline: download, load, and run all enrichments
#'
#' @param AnnDatR AnnDatR object containing clustering results
#' @param universe Character vector of background genes (default: NULL, all genes in clustering data)
#' @param db_loc Directory to store annotation databases (default: "databases")
#' @param hpa_version Version of the Human Protein Atlas to download (default: 24)
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @return A named list with elements: database_enrichment, kegg_enrichment, go_enrichment
#' @export
#'
#' @examples
#' \dontrun{
#' # Run clustering pipeline
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 6.3)
#'
#' # Run annotation pipeline
#' enrichment_results <- hc_annotate(adata_res)
#' head(enrichment_results)
#' }
hc_annotate <- function(
  AnnDatR,
  universe = NULL,
  db_loc = "databases",
  hpa_version = 24,
  verbose = TRUE
) {
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

  kegg_enrichment <- run_kegg_enrichment(
    AnnDatR = AnnDatR,
    universe = universe,
    verbose = verbose
  )
  if (verbose) {
    message("KEGG enrichment done.")
  }

  go_enrichment <- run_go_enrichment(
    AnnDatR = AnnDatR,
    universe = universe,
    verbose = verbose
  )
  if (verbose) {
    message("GO enrichment (with simplification) done.")
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
