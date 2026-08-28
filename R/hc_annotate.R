#' Master annotation pipeline: download, load, and run all enrichments
#'
#' `hc_annotate()` is a master function that performs functional enrichment analysis. It downloads necessary annotation databases,
#' formats them, and runs enrichment analyses for Gene Ontology (GO), KEGG pathways, and custom databases (Human Protein Atlas, Reactome, TRRUST, PanglaoDB)
#' based on the clustering results stored in the AnnDatR object. It also generates visualizations such as treemaps for GO terms and bubblemaps for all enrichment results.
#'
#' @param AnnDatR AnnDatR object containing clustering results
#' @param dbs Character vector of databases to use for enrichment ("GO", "KEGG", "Others") (default: all)
#' @param universe Character vector of background genes (default: NULL, all genes in clustering data)
#' @param simplify_go Logical, group semantically similar GO terms and plot Tree map
#' @param db_loc Directory to store annotation databases (default: "databases")
#' @param hpa_version Human Protein Atlas release to download. `NULL` (default)
#'   uses the current release; a number such as `24` pins an archived release.
#' @param verbose Logical, print progress messages (default: TRUE)
#'
#' @details In case some of the databases did not download correctly, please rerun the function.
#' If that does not help, please manually download the databases, place them in the specified db_loc directory and unzip the zipped files.
#' Links:
#' - Human Protein Atlas (current release): <https://www.proteinatlas.org/download/proteinatlas.tsv.zip>
#' - Human Protein Atlas (archived release `n`): `https://vn.proteinatlas.org/download/proteinatlas.tsv.zip`
#' - Reactome: <https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt>
#' - TRRUST: <https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv>
#' - PanglaoDB: <https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz>
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
#' @examplesIf requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)
#' \donttest{
#' # Run clustering pipeline
#' adata_res <- hc_pca(example_adata, components = 40)
#' adata_res <- hc_distance(adata_res, components = 20)
#' adata_res <- hc_snn(adata_res, neighbors = 15)
#' adata_res <- hc_cluster_consensus(adata_res, resolution = 8, n_seeds = 20)
#'
#' # Enrichment against KEGG. This queries the KEGG web service, so it needs
#' # network access and takes a while.
#' enrichment_results <- hc_annotate(adata_res, dbs = "KEGG")
#' head(enrichment_results$enrichment)
#' }
hc_annotate <- function(
  AnnDatR,
  dbs = c("GO", "KEGG", "Others"),
  universe = NULL,
  simplify_go = TRUE,
  db_loc = "databases",
  hpa_version = NULL,
  verbose = TRUE
) {
  dbs <- match.arg(dbs, c("GO", "KEGG", "Others"), several.ok = TRUE)

  check_installed("clusterProfiler", "to run enrichment analysis", bioc = TRUE)
  check_installed("org.Hs.eg.db", "to map gene identifiers", bioc = TRUE)
  if ("Others" %in% dbs) {
    check_installed("readr", "to read the annotation databases")
  }
  if ("GO" %in% dbs && isTRUE(simplify_go)) {
    check_installed("rrvgo", "to simplify GO terms", bioc = TRUE)
  }

  if (is.null(AnnDatR[["uns"]][["consensus_clustering"]])) {
    stop(
      paste0(
        "AnnDatR$uns$consensus_clustering not found. ",
        "Call `hc_cluster_consensus()` before `hc_annotate()`."
      ),
      call. = FALSE
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
    if (nrow(database_enrichment) > 0L) {
      bubblemap_others <- plot_enrichment_bubblemap(database_enrichment)
    }
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
    if (nrow(kegg_enrichment) > 0L) {
      bubblemap_kegg <- plot_enrichment_bubblemap(kegg_enrichment)
    }
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
    if (isTRUE(simplify_go) && nrow(go_enrichment) > 0L) {
      if (verbose) {
        message("Start GO enrichment simplification...")
      }
      res_go <- reduce_go_terms(go_enrichment, verbose = verbose)
      go_enrichment <- res_go[["combined"]]
      if (nrow(res_go[["reducedTerms"]]) > 0L) {
        treemaps <- plot_enrichment_treemap(res_go[["reducedTerms"]])
      }
      rm(res_go)
    }
    if (nrow(go_enrichment) > 0L) {
      bubblemap_go <- plot_enrichment_bubblemap(go_enrichment)
    }
    if (verbose) {
      message("GO enrichment done.")
    }
  } else {
    go_enrichment <- dplyr::tibble()
  }

  all_enrichment <- dplyr::bind_rows(
    database_enrichment,
    kegg_enrichment,
    go_enrichment
  )

  if (nrow(all_enrichment) == 0L) {
    warning("No significant enrichment found for any cluster.", call. = FALSE)
    return(list(enrichment = empty_enrichment_result()))
  }

  all_enrichment <- all_enrichment |>
    dplyr::filter(!!rlang::sym("Adjusted P-value") < 0.05) |>
    dplyr::arrange(
      !!rlang::sym("Cluster ID"),
      !!rlang::sym("Adjusted P-value")
    )
  # Map Ensembl IDs to Gene Symbols for better interpretability
  all_enrichment <- map_ensembl_to_symbol(all_enrichment)

  result <- list(enrichment = all_enrichment)
  if ("GO" %in% dbs && nrow(go_enrichment) > 0L) {
    if (isTRUE(simplify_go) && exists("treemaps", inherits = FALSE)) {
      result[["treemaps"]] <- treemaps
    }
    result[["bubblemap_go"]] <- bubblemap_go
  }
  if ("KEGG" %in% dbs && nrow(kegg_enrichment) > 0L) {
    result[["bubblemap_kegg"]] <- bubblemap_kegg
  }
  if ("Others" %in% dbs && nrow(database_enrichment) > 0L) {
    result[["bubblemap_others"]] <- bubblemap_others
  }
  result
}
