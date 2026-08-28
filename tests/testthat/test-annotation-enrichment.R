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

# --- reduce_go_terms: rrvgo drops terms it cannot resolve -------------------

go_fixture <- function(term_ids, cluster = "1",
                       database = "GO analysis Biological Process",
                       adjusted_p = 1e-4) {
  n <- length(term_ids)
  tibble::tibble(
    `Cluster ID` = rep(cluster, n),
    Database = rep(database, n),
    `Term ID` = term_ids,
    Term = paste0("term", seq_len(n)),
    GeneRatio = rep("5/50", n),
    BgRatio = rep("50/1000", n),
    `P-value` = rep(1e-5, n),
    `Adjusted P-value` = rep(adjusted_p, n),
    `Gene IDs` = rep("A/B", n)
  )
}

test_that("reduce_go_terms survives when rrvgo resolves fewer than two terms", {
  skip_if_not_installed("rrvgo")
  skip_if_not_installed("org.Hs.eg.db")

  # Regression: calculateSimMatrix() silently drops unresolvable GO ids and
  # returns a bare numeric for a single survivor, so reduceSimMatrix() reached
  # hclust() with one object and failed with
  # "must have n >= 2 objects to cluster". Guarding on the number of *input*
  # terms was not enough.
  cases <- list(
    "one of two resolves" = go_fixture(c("GO:0006915", "GO:9999999")),
    "neither resolves" = go_fixture(c("GO:9999998", "GO:9999999")),
    "ontology mismatch" = go_fixture(
      c("GO:0006915", "GO:0007049"),
      database = "GO analysis Cellular Component"
    )
  )

  for (label in names(cases)) {
    result <- expect_no_error(suppressMessages(reduce_go_terms(cases[[label]])))
    expect_setequal(names(result), c("combined", "reducedTerms"))
    expect_equal(nrow(result[["reducedTerms"]]), 0L, info = label)
    # The enrichment itself must survive untouched.
    expect_equal(nrow(result[["combined"]]), nrow(cases[[label]]), info = label)
  }
})

test_that("reduce_go_terms returns a typed empty table when nothing reduces", {
  skip_if_not_installed("rrvgo")
  skip_if_not_installed("org.Hs.eg.db")

  result <- suppressMessages(reduce_go_terms(go_fixture(c("GO:9999998", "GO:9999999"))))

  expect_equal(
    names(result[["reducedTerms"]]),
    c("go", "parent", "parentTerm", "term", "score", "cluster", "ontology")
  )
  # Regression: an empty simplified table used to be all-logical, so binding it
  # to the character columns of the enrichment failed.
  expect_type(result[["combined"]][["Cluster ID"]], "character")
  expect_no_error(plot_enrichment_treemap(result[["reducedTerms"]]))
})

test_that("reduce_go_terms reduces what it can and skips the rest", {
  skip_if_not_installed("rrvgo")
  skip_if_not_installed("org.Hs.eg.db")

  mixed <- dplyr::bind_rows(
    go_fixture(c("GO:0006915", "GO:0007049", "GO:0006355"), cluster = "1"),
    go_fixture(c("GO:9999998", "GO:9999999"), cluster = "2")
  )

  result <- suppressMessages(reduce_go_terms(mixed))

  expect_gt(nrow(result[["reducedTerms"]]), 0L)
  expect_equal(unique(result[["reducedTerms"]][["cluster"]]), "1")
  expect_gte(nrow(result[["combined"]]), nrow(mixed))
  expect_false(anyNA(result[["reducedTerms"]][["ontology"]]))
})

test_that("reduce_go_terms reports the groups it skipped", {
  skip_if_not_installed("rrvgo")
  skip_if_not_installed("org.Hs.eg.db")

  expect_message(
    reduce_go_terms(go_fixture(c("GO:9999998", "GO:9999999"))),
    "Skipped GO simplification"
  )
  expect_no_message(
    reduce_go_terms(go_fixture(c("GO:9999998", "GO:9999999")), verbose = FALSE)
  )
})

test_that("reduce_go_terms tolerates an adjusted p-value of zero", {
  skip_if_not_installed("rrvgo")
  skip_if_not_installed("org.Hs.eg.db")

  # -log10(0) is infinite and would always win parent selection.
  underflowed <- go_fixture(c("GO:0006915", "GO:0007049", "GO:0006355"))
  underflowed[["Adjusted P-value"]][1] <- 0

  result <- suppressMessages(reduce_go_terms(underflowed))
  expect_gt(nrow(result[["reducedTerms"]]), 0L)
  expect_true(all(is.finite(result[["reducedTerms"]][["score"]])))
})

test_that("reduce_go_terms names the columns it is missing", {
  skip_if_not_installed("rrvgo")
  skip_if_not_installed("org.Hs.eg.db")

  # Regression: this used to be a bare stop() with no message at all.
  expect_error(
    reduce_go_terms(tibble::tibble(a = 1)),
    "missing the required column"
  )
})
