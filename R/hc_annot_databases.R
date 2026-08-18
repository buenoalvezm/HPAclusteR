#' Source URLs for the downloadable annotation databases
#'
#' @param hpa_version Human Protein Atlas version. Use `NULL` (the default in
#'   [get_annot_dbs()]) for the current release; a number selects an archived
#'   release host such as `v24.proteinatlas.org`.
#'
#' @returns A named character vector of URLs.
#' @keywords internal
#' @noRd
annotation_db_urls <- function(hpa_version = NULL) {
  hpa_host <- if (is.null(hpa_version)) {
    "https://www.proteinatlas.org"
  } else {
    sprintf("https://v%s.proteinatlas.org", hpa_version)
  }

  c(
    Human_Protein_Atlas = paste0(hpa_host, "/download/proteinatlas.tsv.zip"),
    panglao_cellmarkers = "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz",
    reactome = "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt",
    trrust = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
  )
}

#' Decompress a gzip file
#'
#' @param path Path to the `.gz` file.
#' @param destination Path to write the decompressed file to.
#'
#' @returns `NULL`, invisibly.
#' @keywords internal
#' @noRd
decompress_gz <- function(path, destination) {
  input <- gzfile(path, "rb")
  on.exit(close(input), add = TRUE)
  output <- file(destination, "wb")
  on.exit(close(output), add = TRUE)

  repeat {
    chunk <- readBin(input, what = "raw", n = 65536L)
    if (length(chunk) == 0L) {
      break
    }
    writeBin(chunk, output)
  }

  invisible(NULL)
}

#' Download annotation databases
#'
#' @param db_loc Destination directory to save the downloaded databases
#'   (default: `"databases"`).
#' @param hpa_version Human Protein Atlas release to download. `NULL` (default)
#'   uses the current release; a number such as `24` pins an archived release.
#' @param timeout Download timeout in seconds (default: 1800). These files are
#'   large, and R's 60 second default is not enough.
#'
#' @returns A named list mapping database IDs to file paths.
#' @keywords internal
get_annot_dbs <- function(
  db_loc = "databases",
  hpa_version = NULL,
  timeout = 1800
) {
  dir.create(db_loc, showWarnings = FALSE, recursive = TRUE)

  previous_timeout <- options(timeout = timeout)
  on.exit(options(previous_timeout), add = TRUE)

  urls <- annotation_db_urls(hpa_version)
  failures <- character(0)

  for (db in names(urls)) {
    url <- urls[[db]]
    destination <- file.path(db_loc, basename(url))
    uncompressed <- sub("\\.(gz|zip)$", "", destination, ignore.case = TRUE)

    if (file.exists(destination) || file.exists(uncompressed)) {
      message(sprintf("%s already exists, skipping.", db))
      next
    }

    message(sprintf("Downloading %s...", db))
    downloaded <- tryCatch(
      {
        utils::download.file(url, destfile = destination, mode = "wb", quiet = TRUE)
        TRUE
      },
      error = function(e) {
        message(sprintf("Failed to download %s: %s", db, conditionMessage(e)))
        FALSE
      },
      warning = function(w) {
        message(sprintf("Failed to download %s: %s", db, conditionMessage(w)))
        FALSE
      }
    )

    if (!downloaded) {
      failures <- c(failures, sprintf("%s (%s)", db, url))
      if (file.exists(destination)) {
        file.remove(destination)
      }
      next
    }

    extracted <- tryCatch(
      {
        if (grepl("\\.zip$", destination, ignore.case = TRUE)) {
          message(sprintf("Unzipping %s...", basename(destination)))
          utils::unzip(destination, exdir = db_loc)
          file.remove(destination)
        } else if (grepl("\\.gz$", destination, ignore.case = TRUE)) {
          message(sprintf("Decompressing %s...", basename(destination)))
          decompress_gz(destination, uncompressed)
          file.remove(destination)
        }
        TRUE
      },
      error = function(e) {
        message(sprintf("Failed to extract %s: %s", db, conditionMessage(e)))
        FALSE
      }
    )

    if (!extracted) {
      failures <- c(failures, sprintf("%s (%s)", db, url))
    }
  }

  if (length(failures) > 0L) {
    warning(
      "Could not obtain the following annotation databases:\n  ",
      paste(failures, collapse = "\n  "),
      "\nEnrichment will run on the databases that were downloaded. ",
      "Download the missing files manually into '",
      db_loc,
      "' and re-run, or restrict `dbs` to \"GO\" and \"KEGG\".",
      call. = FALSE
    )
  }

  db_files <- list()

  hpa_file <- file.path(db_loc, "proteinatlas.tsv")
  if (file.exists(hpa_file)) {
    hpa_databases <- c(
      "secretome_location",
      "specificity_blood",
      "specificity_brain",
      "specificity_tissue",
      "specificity_celline",
      "specificity_singlecell",
      "subcellular_location",
      "protein_class"
    )
    for (db in hpa_databases) {
      db_files[[db]] <- hpa_file
    }
  }

  other_files <- c(
    panglao_cellmarkers = "PanglaoDB_markers_27_Mar_2020.tsv",
    reactome = "Ensembl2Reactome_All_Levels.txt",
    trrust = "trrust_rawdata.human.tsv"
  )
  for (db in names(other_files)) {
    path <- file.path(db_loc, other_files[[db]])
    if (file.exists(path)) {
      db_files[[db]] <- path
    }
  }

  db_files
}

#' Read a raw annotation database file
#'
#' Applies the per-source reading options. Reactome and TRRUST ship without a
#' header row, so their column names must be supplied explicitly; reading
#' TRRUST with `col_names = TRUE` silently consumed its first record and named
#' the columns after whichever transcription factor happened to sort first.
#'
#' @param db_id The database ID.
#' @param path Path to the file.
#'
#' @returns A tibble of the raw file contents.
#' @keywords internal
#' @noRd
read_annotation_source <- function(db_id, path) {
  check_installed("readr", "to read the annotation databases")

  if (db_id == "reactome") {
    return(readr::read_tsv(
      path,
      col_names = c(
        "ensg_id",
        "id",
        "URL",
        "description",
        "evidence_code",
        "species"
      ),
      col_types = readr::cols(.default = readr::col_character())
    ))
  }

  if (db_id == "trrust") {
    return(readr::read_tsv(
      path,
      col_names = c("tf", "target", "interaction", "pubmed_id"),
      col_types = readr::cols(.default = readr::col_character())
    ))
  }

  readr::read_tsv(path, show_col_types = FALSE)
}

#' Map Human Protein Atlas gene symbols to Ensembl gene IDs
#'
#' @param hpa The Human Protein Atlas table.
#'
#' @returns A tibble with `gene_name` and `ensg_id`.
#' @keywords internal
#' @noRd
hpa_symbol_map <- function(hpa) {
  hpa |>
    dplyr::select(
      gene_name = !!rlang::sym("Gene"),
      ensg_id = !!rlang::sym("Ensembl")
    ) |>
    dplyr::distinct()
}

#' Reshape a Human Protein Atlas classification column into long form
#'
#' @param db The Human Protein Atlas table.
#' @param column Name of the column holding the classification.
#' @param cluster_genes Character vector of Ensembl IDs to keep.
#' @param split Logical, whether to split multi-valued entries on `;`.
#'
#' @returns A tibble with `ensg_id`, `term` and `term_id`.
#' @keywords internal
#' @noRd
format_hpa_column <- function(db, column, cluster_genes, split = TRUE) {
  result <- db |>
    dplyr::select(
      ensg_id = !!rlang::sym("Ensembl"),
      term = !!rlang::sym(column)
    ) |>
    dplyr::filter(!is.na(!!rlang::sym("term")))

  if (split) {
    result <- result |>
      dplyr::mutate(term = gsub(", ", ";", !!rlang::sym("term"))) |>
      tidyr::separate_longer_delim(!!rlang::sym("term"), delim = ";")
  }

  result |>
    dplyr::filter(!is.na(!!rlang::sym("term"))) |>
    dplyr::mutate(term_id = !!rlang::sym("term")) |>
    dplyr::filter(!!rlang::sym("ensg_id") %in% cluster_genes)
}

#' Formatting functions for annotation databases (excluding GO)
#'
#' Each function takes `(db, hpa, cluster_genes)` and returns a tibble with the
#' standardised columns `ensg_id`, `term` and `term_id`. `hpa` is the already
#' loaded Human Protein Atlas table, shared by every formatter that needs it.
#' @keywords internal
database_formatting_functions <- list(
  secretome_location = function(db, hpa, cluster_genes) {
    format_hpa_column(db, "Secretome location", cluster_genes, split = FALSE)
  },
  specificity_blood = function(db, hpa, cluster_genes) {
    format_hpa_column(db, "RNA blood cell specificity", cluster_genes)
  },
  specificity_brain = function(db, hpa, cluster_genes) {
    format_hpa_column(db, "RNA single nuclei brain specificity", cluster_genes)
  },
  specificity_tissue = function(db, hpa, cluster_genes) {
    format_hpa_column(db, "RNA tissue specificity", cluster_genes)
  },
  specificity_celline = function(db, hpa, cluster_genes) {
    format_hpa_column(db, "RNA cell line specificity", cluster_genes)
  },
  specificity_singlecell = function(db, hpa, cluster_genes) {
    format_hpa_column(db, "RNA single cell type specificity", cluster_genes)
  },
  subcellular_location = function(db, hpa, cluster_genes) {
    format_hpa_column(db, "Subcellular location", cluster_genes)
  },
  protein_class = function(db, hpa, cluster_genes) {
    db |>
      dplyr::select(
        ensg_id = !!rlang::sym("Ensembl"),
        term = !!rlang::sym("Protein class")
      ) |>
      tidyr::separate_longer_delim(!!rlang::sym("term"), delim = ", ") |>
      dplyr::filter(!is.na(!!rlang::sym("term"))) |>
      dplyr::mutate(term_id = !!rlang::sym("term")) |>
      dplyr::filter(!!rlang::sym("ensg_id") %in% cluster_genes)
  },
  panglao_cellmarkers = function(db, hpa, cluster_genes) {
    db |>
      dplyr::filter(
        grepl("Hs", !!rlang::sym("species")),
        !!rlang::sym("gene type") %in%
          c("protein-coding gene", "protein coding gene")
      ) |>
      dplyr::select(gene_name = 2, term = 3) |>
      dplyr::distinct() |>
      dplyr::inner_join(hpa_symbol_map(hpa), by = "gene_name") |>
      dplyr::filter(!is.na(!!rlang::sym("ensg_id"))) |>
      dplyr::mutate(term_id = !!rlang::sym("term")) |>
      dplyr::select(dplyr::all_of(c("ensg_id", "term", "term_id"))) |>
      dplyr::filter(!!rlang::sym("ensg_id") %in% cluster_genes)
  },
  reactome = function(db, hpa, cluster_genes) {
    db |>
      dplyr::filter(!!rlang::sym("species") == "Homo sapiens") |>
      dplyr::select(
        ensg_id = !!rlang::sym("ensg_id"),
        term = !!rlang::sym("description"),
        term_id = !!rlang::sym("id")
      ) |>
      dplyr::distinct() |>
      dplyr::filter(!!rlang::sym("ensg_id") %in% cluster_genes)
  },
  trrust = function(db, hpa, cluster_genes) {
    # Genes are annotated with the transcription factors that regulate them.
    db |>
      dplyr::select(
        gene_name = !!rlang::sym("target"),
        term = !!rlang::sym("tf")
      ) |>
      dplyr::distinct() |>
      dplyr::inner_join(hpa_symbol_map(hpa), by = "gene_name") |>
      dplyr::filter(!is.na(!!rlang::sym("ensg_id"))) |>
      dplyr::mutate(term_id = !!rlang::sym("term")) |>
      dplyr::select(dplyr::all_of(c("ensg_id", "term", "term_id"))) |>
      dplyr::filter(!!rlang::sym("ensg_id") %in% cluster_genes)
  }
)

#' Format an annotation database (excluding GO)
#'
#' @param raw_db A tibble containing the raw database data.
#' @param db_id The ID of the database to format (e.g. `"secretome_location"`).
#' @param hpa The loaded Human Protein Atlas table, or `NULL` when not needed.
#' @param cluster_genes Character vector of Ensembl IDs to keep.
#'
#' @return A tibble with standardized columns: `ensg_id`, `term`, `term_id`.
#' @keywords internal
format_annotation_database <- function(raw_db, db_id, hpa, cluster_genes) {
  db_format_function <- database_formatting_functions[[db_id]]
  if (is.null(db_format_function)) {
    stop(
      sprintf("No formatting function found for database ID: %s", db_id),
      call. = FALSE
    )
  }

  formatted_db <- db_format_function(raw_db, hpa, cluster_genes)
  stopifnot(identical(names(formatted_db), c("ensg_id", "term", "term_id")))

  formatted_db
}

#' Load and format a single annotation database
#'
#' @param db_id The ID of the database to load (e.g. `"reactome"`).
#' @param db_files Named list mapping database IDs to file paths.
#' @param AnnDatR AnnDatR object containing clustering results.
#'
#' @return A tibble with standardized columns: `ensg_id`, `term`, `term_id`.
#' @keywords internal
load_annotation_database <- function(db_id, db_files, AnnDatR) {
  cluster_genes <- AnnDatR[["uns"]][["consensus_clustering"]][["gene"]]
  raw_db <- read_annotation_source(db_id, db_files[[db_id]])

  hpa <- if (db_id %in% c("panglao_cellmarkers", "trrust")) {
    read_annotation_source(
      "secretome_location",
      db_files[["secretome_location"]]
    )
  } else {
    NULL
  }

  format_annotation_database(raw_db, db_id, hpa, cluster_genes)
}

#' Build annotation terms tibble from all databases in db_files
#'
#' @param db_files Named list mapping database IDs to file paths (as returned
#'   by [get_annot_dbs()]).
#' @param AnnDatR AnnDatR object containing clustering results.
#'
#' @details
#' Eight of the supported databases are different columns of the same Human
#' Protein Atlas table, and two more need it to map gene symbols to Ensembl
#' IDs. That file is several hundred megabytes, so it is read once here and
#' shared across all formatters rather than re-read for each database.
#'
#' @return A tibble with columns: `db_id`, `ensg_id`, `term`, `term_id`.
#' @keywords internal
build_annotation_terms_tibble <- function(db_files, AnnDatR) {
  check_installed("readr", "to read the annotation databases")

  if (is.null(AnnDatR[["uns"]][["consensus_clustering"]])) {
    stop(
      paste0(
        "AnnDatR$uns$consensus_clustering not found. ",
        "Call `hc_cluster_consensus()` before `hc_annotate()`."
      ),
      call. = FALSE
    )
  }

  cluster_genes <- AnnDatR[["uns"]][["consensus_clustering"]][["gene"]]
  db_ids <- names(db_files)

  # Read the Human Protein Atlas table at most once.
  hpa_path <- db_files[["secretome_location"]]
  hpa <- if (!is.null(hpa_path)) {
    read_annotation_source("secretome_location", hpa_path)
  } else {
    NULL
  }

  hpa_backed <- names(database_formatting_functions)[
    vapply(
      names(database_formatting_functions),
      function(db_id) identical(db_files[[db_id]], hpa_path),
      logical(1)
    )
  ]

  results <- lapply(db_ids, function(db_id) {
    raw_db <- if (db_id %in% hpa_backed && !is.null(hpa)) {
      hpa
    } else {
      read_annotation_source(db_id, db_files[[db_id]])
    }

    formatted_db <- format_annotation_database(
      raw_db,
      db_id,
      hpa,
      cluster_genes
    )
    formatted_db[["db_id"]] <- db_id
    formatted_db
  })

  dplyr::bind_rows(results) |>
    dplyr::select(dplyr::all_of(c("db_id", "ensg_id", "term", "term_id")))
}
