enrichment_fixture <- function() {
  tibble::tibble(
    `Cluster ID` = c("1", "1", "2", "2", "3"),
    Database = c(
      "GO analysis Biological Process",
      "GO analysis Molecular Function",
      "GO analysis Biological Process",
      "GO analysis Cellular Component",
      "GO analysis Biological Process"
    ),
    `Term ID` = c("GO:0001", "GO:0002", "GO:0003", "GO:0004", "GO:0005"),
    Term = c("apoptosis", "binding", "cell cycle", "membrane", "transport"),
    GeneRatio = c("5/50", "10/50", "3/40", "8/40", "2/30"),
    BgRatio = c("50/1000", "200/1000", "60/1000", "100/1000", "40/1000"),
    `P-value` = c(1e-5, 1e-3, 1e-4, 1e-2, 1e-6),
    `Adjusted P-value` = c(1e-4, 1e-2, 1e-3, 4e-2, 1e-5),
    `Gene IDs` = c("A/B/C", "D/E", "F", "G/H", "I")
  )
}

test_that("plot_enrichment_bubblemap computes odds ratios from the ratio strings", {
  plot <- plot_enrichment_bubblemap(enrichment_fixture())

  expect_s3_class(plot, "ggplot")
  built <- ggplot2::ggplot_build(plot)

  data <- plot[["data"]]
  # 5/50 over 50/1000 is an odds ratio of 2.
  expect_equal(data[["odds_ratio"]][1], (5 / 50) / (50 / 1000))
  expect_true(all(data[["capped_odds_ratio"]] <= 50))
  expect_no_error(built)
})

test_that("plot_enrichment_bubblemap caps the odds ratio", {
  extreme <- enrichment_fixture()
  extreme[["GeneRatio"]][1] <- "50/50"
  extreme[["BgRatio"]][1] <- "1/1000"

  plot <- plot_enrichment_bubblemap(extreme, odds_ratio_cap = 10)
  expect_lte(max(plot[["data"]][["capped_odds_ratio"]]), 10)
  expect_gt(max(plot[["data"]][["odds_ratio"]]), 10)
})

test_that("plot_enrichment_treemap returns one plot per cluster", {
  skip_if_not_installed("treemapify")

  reduced <- tibble::tibble(
    go = c("GO:0001", "GO:0002", "GO:0003"),
    term = c("apoptosis", "binding", "cell cycle"),
    parentTerm = c("death", "death", "cycle"),
    score = c(4, 2, 3),
    cluster = c("1", "1", "2"),
    ontology = c("BP", "MF", "BP")
  )

  plots <- plot_enrichment_treemap(reduced)

  expect_setequal(names(plots), c("1", "2"))
  expect_s3_class(plots[["1"]], "ggplot")
  expect_no_error(ggplot2::ggplot_build(plots[["1"]]))
})

test_that("reduce_go_terms rejects input without the required columns", {
  skip_if_not_installed("rrvgo")
  skip_if_not_installed("org.Hs.eg.db")

  expect_error(reduce_go_terms(tibble::tibble(a = 1)))
})

test_that("empty_enrichment_result has the standard columns", {
  empty <- empty_enrichment_result()

  expect_equal(nrow(empty), 0L)
  expect_equal(
    names(empty),
    c(
      "Cluster ID", "Database", "Term ID", "Term", "GeneRatio",
      "BgRatio", "P-value", "Adjusted P-value", "Gene IDs"
    )
  )
  # It must survive the operations hc_annotate() performs on results.
  expect_no_error(dplyr::filter(empty, !!rlang::sym("Adjusted P-value") < 0.05))
  expect_no_error(dplyr::bind_rows(empty, enrichment_fixture()))
})

test_that("run_database_enrichment returns an empty result when nothing enriches", {
  skip_if_not_installed("clusterProfiler")

  # Terms that no cluster can be enriched for: one gene, one unique term each.
  clustering <- adata_clustered()[["uns"]][["consensus_clustering"]]
  annotation_terms <- tibble::tibble(
    db_id = "protein_class",
    ensg_id = clustering[["gene"]],
    term = paste0("unique_", seq_len(nrow(clustering))),
    term_id = paste0("unique_", seq_len(nrow(clustering)))
  )

  result <- suppressWarnings(suppressMessages(
    run_database_enrichment(annotation_terms, adata_clustered(), verbose = FALSE)
  ))

  # Regression: this used to fail with "object 'db_id' not found".
  expect_equal(nrow(result), 0L)
  expect_equal(names(result), names(empty_enrichment_result()))
})

test_that("the enrichment entry points require a consensus clustering", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  expect_error(
    run_go_enrichment(adata_snn()),
    "Call `hc_cluster_consensus\\(\\)`"
  )
  expect_error(
    run_kegg_enrichment(adata_snn()),
    "Call `hc_cluster_consensus\\(\\)`"
  )
  expect_error(
    run_database_enrichment(tibble::tibble(db_id = character()), adata_snn()),
    "Call `hc_cluster_consensus\\(\\)`"
  )
})

test_that("hc_annotate rejects unknown databases and unclustered objects", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  expect_error(hc_annotate(adata_clustered(), dbs = "Nonsense"))
  expect_error(
    hc_annotate(adata_snn(), dbs = "KEGG"),
    "Call `hc_cluster_consensus\\(\\)`"
  )
})

test_that("run_database_enrichment finds a planted enrichment", {
  skip_if_not_installed("clusterProfiler")

  # One cluster's genes all carry the same term. Every other gene gets a
  # filler term so that the annotated background is the whole gene set --
  # clusterProfiler restricts the universe to annotated genes.
  clustering <- adata_clustered()[["uns"]][["consensus_clustering"]]
  target_cluster <- clustering[["cluster"]][1]
  target <- clustering[["gene"]][clustering[["cluster"]] == target_cluster]
  others <- setdiff(clustering[["gene"]], target)

  annotation_terms <- dplyr::bind_rows(
    tibble::tibble(
      db_id = "protein_class",
      ensg_id = target,
      term = "planted term",
      term_id = "planted term"
    ),
    tibble::tibble(
      db_id = "protein_class",
      ensg_id = others,
      term = "filler term",
      term_id = "filler term"
    )
  )

  result <- suppressWarnings(run_database_enrichment(
    annotation_terms,
    adata_clustered(),
    verbose = FALSE
  ))

  expect_true("planted term" %in% result[["Term"]])
  expect_true(all(
    c("Cluster ID", "Database", "Term ID", "Term", "P-value", "Adjusted P-value") %in%
      names(result)
  ))
  planted <- result[result[["Term"]] == "planted term", ]
  expect_lt(min(planted[["Adjusted P-value"]]), 0.05)
  # It is enriched in the cluster it was planted in.
  expect_equal(planted[["Cluster ID"]][which.min(planted[["Adjusted P-value"]])],
               target_cluster)
})

test_that("map_ensembl_to_symbol adds a gene symbol column", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")

  input <- tibble::tibble(`Gene IDs` = "ENSG00000141510/ENSG00000087088")
  result <- suppressMessages(suppressWarnings(map_ensembl_to_symbol(input)))

  expect_true("Gene names" %in% names(result))
  expect_type(result[["Gene names"]], "character")
  # TP53 and BAX are the corresponding symbols.
  expect_match(result[["Gene names"]], "TP53")
})
