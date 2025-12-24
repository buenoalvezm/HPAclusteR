#' Master annotation pipeline: download, load, and run all enrichments
#'
#' `hc_annotate()` is a master function that performs functional enrichment analysis. It downloads necessary annotation databases,
#' formats them, and runs enrichment analyses for Gene Ontology (GO), KEGG pathways, and custom databases (Human Protein Atlas, Reactome, TRRUST, PanglaoDB)
#' based on the clustering results stored in the AnnDatR object. It also generates visualizations such as treemaps for GO terms and bubblemaps for all enrichment results.
#'
#' @param AnnDatR AnnDatR object containing clustering results
#' @param dbs Character vector of databases to use for enrichment ("GO", "KEGG", "Others") (default: all)
#' @param universe Character vector of background genes (default: NULL, all genes in clustering data)
#' @param db_loc Directory to store annotation databases (default: "databases")
#' @param hpa_version Version of the Human Protein Atlas to download (default: 24)
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @details In case some of the databases did not download correctly, please rerun the function.
#' If that does not help, please manually download the databases, place them in the specified db_loc directory and unzip the zipped files.
#' Links:
#' - Human Protein Atlas version 24: https://v%d.proteinatlas.org/download/proteinatlas.tsv.zip
#' - Reactome: https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt
#' - Trrust: https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv
#' - Panglao: https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz
#'
#' Or just run the annotation only using KEGG and GO by setting dbs = c("KEGG", "GO").
#'
#' This function will require time depending on the size of the dataset and the number of clusters, especially if GO is included.
#' Set verbose = TRUE to monitor the progress.
#'
#' For the bubblemaps, the odds ratio is capped at 50 for better visualization. Also, only the top 150 terms (by mean odds ratio across clusters) are shown.
#'
#' @return A named list with elements:
#' - enrichment: Data frame with combined enrichment results from all selected databases
#' - treemaps: List of ggplot objects with GO treemaps per cluster (if GO selected)
#' - bubblemap_go: ggplot object with GO enrichment bubblemap (if GO selected)
#' - bubblemap_kegg: ggplot object with KEGG enrichment bubblemap (if KEGG selected)
#' - bubblemap_others: ggplot object with custom database enrichment bubblemap (if Others selected)
#'
#' @export
#'
#' @examples
#' # Run clustering pipeline
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 7)
#'
#' # Run annotation pipeline
#' enrichment_results <- hc_annotate(adata_res, dbs = "KEGG")
#' head(enrichment_results$enrichment)
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
  if (!requireNamespace("rrvgo", quietly = TRUE) && "GO" %in% dbs) {
    stop(
      "The 'rrvgo' package is required for KEGG enrichment. Please install it using BiocManager::install('rrvgo')."
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
    bubblemap_others <- plot_enrichment_bubblemap(database_enrichment)
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
    bubblemap_kegg <- plot_enrichment_bubblemap(kegg_enrichment)
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
    res_go <- reduce_go_terms(go_enrichment)
    if (verbose) {
      message("Start GO enrichment simplification...")
    }
    go_enrichment <- res_go[["combined"]]
    treemaps <- plot_enrichment_treemap(res_go[["reducedTerms"]])
    rm(res_go)
    bubblemap_go <- plot_enrichment_bubblemap(go_enrichment)
    if (verbose) {
      message("GO enrichment (with simplification) done.")
    }
  } else {
    go_enrichment <- dplyr::tibble()
  }

  if (
    !is.null(database_enrichment) ||
      !is.null(kegg_enrichment) ||
      !is.null(go_enrichment)
  ) {
    all_enrichment <- dplyr::bind_rows(
      database_enrichment,
      kegg_enrichment,
      go_enrichment
    ) |>
      dplyr::filter(!!rlang::sym("Adjusted P-value") < 0.05) |>
      dplyr::arrange(
        !!rlang::sym("Cluster ID"),
        !!rlang::sym("Adjusted P-value")
      )
  } else {
    warning("No enrichment results to combine.")
    return(NULL)
  }
  # Map Ensembl IDs to Gene Symbols for better interpretability
  all_enrichment <- map_ensembl_to_symbol(all_enrichment)

  result <- list()
  result$enrichment <- all_enrichment
  if ("GO" %in% dbs && go_enrichment |> nrow() != 0) {
    result$treemaps <- treemaps
    result$bubblemap_go <- bubblemap_go
  }
  if ("KEGG" %in% dbs && kegg_enrichment |> nrow() != 0) {
    result$bubblemap_kegg <- bubblemap_kegg
  }
  if ("Others" %in% dbs && database_enrichment |> nrow() != 0) {
    result$bubblemap_others <- bubblemap_others
  }
  return(result)
}
