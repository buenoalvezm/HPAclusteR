fake_hpa <- function() {
  tibble::tibble(
    Gene = c("TP53", "BAX", "EGFR"),
    Ensembl = c("ENSG001", "ENSG002", "ENSG003"),
    `Secretome location` = c("Secreted to blood", NA, "Intracellular"),
    `RNA tissue specificity` = c(
      "Tissue enhanced (liver, kidney)",
      "Low tissue specificity",
      NA
    ),
    `RNA blood cell specificity` = c("Immune cell enhanced (T-cell)", NA, NA),
    `RNA single nuclei brain specificity` = c("Region enriched (cortex)", NA, NA),
    `RNA cell line specificity` = c("Cell line enhanced (HeLa)", NA, NA),
    `RNA single cell type specificity` = c("Cell type enhanced (hepatocytes)", NA, NA),
    `Subcellular location` = c("Nucleoplasm;Cytosol", "Mitochondria", NA),
    `Protein class` = c(
      "Transcription factors, Disease related genes",
      "Predicted intracellular proteins",
      NA
    )
  )
}

test_that("annotation_db_urls points at the current release by default", {
  urls <- annotation_db_urls()

  expect_setequal(
    names(urls),
    c("Human_Protein_Atlas", "panglao_cellmarkers", "reactome", "trrust")
  )
  expect_equal(
    unname(urls[["Human_Protein_Atlas"]]),
    "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
  )
  expect_true(all(grepl("^https://", urls)))
})

test_that("annotation_db_urls can pin an archived Human Protein Atlas release", {
  expect_equal(
    unname(annotation_db_urls(24)[["Human_Protein_Atlas"]]),
    "https://v24.proteinatlas.org/download/proteinatlas.tsv.zip"
  )
})

test_that("the Human Protein Atlas formatters emit the standard columns", {
  hpa <- fake_hpa()
  genes <- c("ENSG001", "ENSG002", "ENSG003")

  hpa_databases <- c(
    "secretome_location", "specificity_blood", "specificity_brain",
    "specificity_tissue", "specificity_celline", "specificity_singlecell",
    "subcellular_location", "protein_class"
  )

  for (db_id in hpa_databases) {
    result <- format_annotation_database(hpa, db_id, hpa, genes)
    expect_equal(names(result), c("ensg_id", "term", "term_id"), info = db_id)
    expect_false(anyNA(result[["term"]]), info = db_id)
    expect_true(all(result[["ensg_id"]] %in% genes), info = db_id)
    expect_equal(result[["term"]], result[["term_id"]], info = db_id)
  }
})

test_that("multi-valued Human Protein Atlas columns are split into one row per term", {
  result <- format_annotation_database(
    fake_hpa(),
    "subcellular_location",
    fake_hpa(),
    c("ENSG001", "ENSG002")
  )

  expect_setequal(
    result[["term"]][result[["ensg_id"]] == "ENSG001"],
    c("Nucleoplasm", "Cytosol")
  )
})

test_that("protein_class splits on comma-space", {
  result <- format_annotation_database(
    fake_hpa(),
    "protein_class",
    fake_hpa(),
    "ENSG001"
  )
  expect_setequal(
    result[["term"]],
    c("Transcription factors", "Disease related genes")
  )
})

test_that("formatters restrict output to the requested genes", {
  result <- format_annotation_database(
    fake_hpa(),
    "secretome_location",
    fake_hpa(),
    "ENSG003"
  )
  expect_equal(unique(result[["ensg_id"]]), "ENSG003")
})

test_that("the TRRUST formatter annotates targets with their transcription factors", {
  trrust <- tibble::tibble(
    tf = c("AATF", "AATF", "TP53"),
    target = c("BAX", "EGFR", "BAX"),
    interaction = c("Repression", "Unknown", "Activation"),
    pubmed_id = c("22909821", "17157788", "1")
  )

  result <- format_annotation_database(
    trrust,
    "trrust",
    fake_hpa(),
    c("ENSG002", "ENSG003")
  )

  expect_equal(names(result), c("ensg_id", "term", "term_id"))
  # BAX (ENSG002) is regulated by AATF and TP53.
  expect_setequal(result[["term"]][result[["ensg_id"]] == "ENSG002"], c("AATF", "TP53"))
  expect_equal(result[["term"]][result[["ensg_id"]] == "ENSG003"], "AATF")
})

test_that("read_annotation_source names the headerless TRRUST columns", {
  skip_if_not_installed("readr")

  path <- withr::local_tempfile(fileext = ".tsv")
  # TRRUST ships without a header; the first line is a real record.
  writeLines(
    c("AATF\tBAX\tRepression\t22909821", "AATF\tCDKN1A\tUnknown\t17157788"),
    path
  )

  result <- read_annotation_source("trrust", path)

  expect_equal(names(result), c("tf", "target", "interaction", "pubmed_id"))
  # Regression: reading with a header consumed this first record.
  expect_equal(nrow(result), 2L)
  expect_equal(result[["target"]][1], "BAX")
})

test_that("read_annotation_source names the headerless Reactome columns", {
  skip_if_not_installed("readr")

  path <- withr::local_tempfile(fileext = ".txt")
  writeLines(
    paste(
      c("ENSG001", "R-HSA-1", "https://reactome.org/x", "Signal Transduction", "IEA", "Homo sapiens"),
      collapse = "\t"
    ),
    path
  )

  result <- read_annotation_source("reactome", path)
  expect_equal(
    names(result),
    c("ensg_id", "id", "URL", "description", "evidence_code", "species")
  )
  expect_equal(nrow(result), 1L)
})

test_that("the Reactome formatter keeps only human entries", {
  reactome <- tibble::tibble(
    ensg_id = c("ENSG001", "ENSG002", "ENSG001"),
    id = c("R-HSA-1", "R-CEL-2", "R-HSA-3"),
    URL = rep("https://reactome.org", 3),
    description = c("Signal Transduction", "Worm pathway", "Metabolism"),
    evidence_code = rep("IEA", 3),
    species = c("Homo sapiens", "Caenorhabditis elegans", "Homo sapiens")
  )

  result <- format_annotation_database(
    reactome,
    "reactome",
    NULL,
    c("ENSG001", "ENSG002")
  )

  expect_equal(nrow(result), 2L)
  expect_equal(unique(result[["ensg_id"]]), "ENSG001")
  expect_setequal(result[["term_id"]], c("R-HSA-1", "R-HSA-3"))
})

test_that("the PanglaoDB formatter keeps human protein-coding markers", {
  panglao <- tibble::tibble(
    species = c("Hs", "Mm", "Hs Mm"),
    `official gene symbol` = c("TP53", "BAX", "EGFR"),
    `cell type` = c("Hepatocytes", "Neurons", "Fibroblasts"),
    `gene type` = c("protein-coding gene", "protein-coding gene", "protein-coding gene")
  )

  result <- format_annotation_database(
    panglao,
    "panglao_cellmarkers",
    fake_hpa(),
    c("ENSG001", "ENSG002", "ENSG003")
  )

  expect_equal(names(result), c("ensg_id", "term", "term_id"))
  # The mouse-only marker is dropped; both human ones are kept.
  expect_setequal(result[["term"]], c("Hepatocytes", "Fibroblasts"))
})

test_that("format_annotation_database rejects unknown database IDs", {
  expect_error(
    format_annotation_database(fake_hpa(), "not_a_database", NULL, "ENSG001"),
    "No formatting function found"
  )
})

test_that("get_annot_dbs warns when a database could not be obtained", {
  directory <- withr::local_tempdir()
  testthat::local_mocked_bindings(
    annotation_db_urls = function(hpa_version = NULL) {
      c(broken = "https://example.invalid/nope.tsv")
    }
  )
  expect_warning(
    suppressMessages(get_annot_dbs(db_loc = directory)),
    "Could not obtain"
  )
})

test_that("get_annot_dbs maps discovered files to database IDs", {
  directory <- withr::local_tempdir()
  # Pretend every download already happened, so nothing hits the network.
  for (file in c(
    "proteinatlas.tsv",
    "trrust_rawdata.human.tsv",
    "Ensembl2Reactome_All_Levels.txt",
    "PanglaoDB_markers_27_Mar_2020.tsv"
  )) {
    writeLines("x", file.path(directory, file))
  }

  db_files <- suppressMessages(get_annot_dbs(db_loc = directory))

  # All eight Human Protein Atlas views resolve to the same file.
  hpa_views <- c(
    "secretome_location", "specificity_blood", "specificity_brain",
    "specificity_tissue", "specificity_celline", "specificity_singlecell",
    "subcellular_location", "protein_class"
  )
  expect_true(all(hpa_views %in% names(db_files)))
  expect_equal(
    length(unique(unlist(db_files[hpa_views]))),
    1L
  )
  expect_setequal(
    names(db_files),
    c(hpa_views, "panglao_cellmarkers", "reactome", "trrust")
  )
})

test_that("build_annotation_terms_tibble requires a consensus clustering", {
  expect_error(
    build_annotation_terms_tibble(list(), adata_snn()),
    "Call `hc_cluster_consensus\\(\\)`"
  )
})

test_that("build_annotation_terms_tibble reads the shared file once", {
  skip_if_not_installed("readr")

  directory <- withr::local_tempdir()
  hpa_path <- file.path(directory, "proteinatlas.tsv")
  readr::write_tsv(fake_hpa(), hpa_path)

  adata <- adata_clustered()$clone(deep = TRUE)
  adata[["uns"]][["consensus_clustering"]] <- tibble::tibble(
    gene = c("ENSG001", "ENSG002", "ENSG003"),
    cluster = c("1", "1", "2")
  )

  db_files <- list(
    secretome_location = hpa_path,
    subcellular_location = hpa_path,
    protein_class = hpa_path
  )

  reads <- 0L
  testthat::local_mocked_bindings(
    read_annotation_source = function(db_id, path) {
      reads <<- reads + 1L
      fake_hpa()
    }
  )

  result <- build_annotation_terms_tibble(db_files, adata)

  expect_equal(names(result), c("db_id", "ensg_id", "term", "term_id"))
  expect_setequal(
    unique(result[["db_id"]]),
    c("secretome_location", "subcellular_location", "protein_class")
  )
  # Regression: the table used to be re-read once per database.
  expect_equal(reads, 1L)
})
