#' Download annotation databases
#'
#' @param db_loc Destination directory to save the downloaded databases (default: "databases").
#' @param hpa_version Version of the Human Protein Atlas to download (default: 24).
#'
#' @returns Invisible TRUE if all downloads are successful.
#' @keywords internal
get_annot_dbs <- function(
  db_loc = "databases",
  hpa_version = 24
) {
  if (!requireNamespace("utils", quietly = TRUE)) {
    stop(
      "The 'utils' package is required for this function. Please install it using install.packages('utils')."
    )
  }

  dir.create(db_loc, showWarnings = FALSE, recursive = TRUE)

  urls <- list(
    HumanProteinAtlas = sprintf(
      "https://v%d.proteinatlas.org/download/proteinatlas.tsv.zip",
      hpa_version
    ),
    panglao_cellmarkers = "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz",
    reactome = "https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt",
    trrust = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"
  )

  for (db in names(urls)) {
    url <- urls[[db]]
    dest_file <- file.path(db_loc, basename(url))
    uncompressed_file <- sub("\\.gz$", "", dest_file, ignore.case = TRUE)
    uncompressed_file <- sub(
      "\\.zip$",
      "",
      uncompressed_file,
      ignore.case = TRUE
    )
    # Skip if either file exists
    if (file.exists(dest_file) || file.exists(uncompressed_file)) {
      message(sprintf(
        "%s already exists, skipping.",
        db
      ))
      next
    }

    message(sprintf("Downloading %s...", db))
    download_success <- FALSE
    tryCatch(
      {
        utils::download.file(
          url,
          destfile = dest_file,
          mode = "wb",
          quiet = TRUE
        )
        download_success <- TRUE
      },
      error = function(e) {
        message(sprintf("Failed to download %s: %s", db, e$message))
      }
    )
    if (download_success) {
      # Unzip .zip files and delete after extraction
      if (
        grepl("\\.zip$", dest_file, ignore.case = TRUE) &&
          file.exists(dest_file)
      ) {
        message(sprintf("Unzipping %s...", dest_file))
        unzip_success <- FALSE
        tryCatch(
          {
            utils::unzip(dest_file, exdir = db_loc)
            unzip_success <- TRUE
          },
          error = function(e) {
            message(sprintf("Failed to unzip %s: %s", dest_file, e$message))
          }
        )
        if (unzip_success) {
          file.remove(dest_file)
          message(sprintf("Deleted zip file: %s", dest_file))
        }
      }
      # Decompress .gz files and delete after extraction
      if (
        grepl("\\.gz$", dest_file, ignore.case = TRUE) &&
          file.exists(dest_file)
      ) {
        message(sprintf("Decompressing %s...", dest_file))
        decompress_success <- FALSE
        tryCatch(
          {
            out_file <- sub("\\.gz$", "", dest_file, ignore.case = TRUE)
            con_in <- gzfile(dest_file, "rb")
            con_out <- file(out_file, "wb")
            while (
              length(buf <- readBin(con_in, what = "raw", n = 65536)) > 0
            ) {
              writeBin(buf, con_out)
            }
            close(con_in)
            close(con_out)
            decompress_success <- TRUE
          },
          error = function(e) {
            message(sprintf(
              "Failed to decompress %s: %s",
              dest_file,
              e$message
            ))
          }
        )
        if (decompress_success) {
          file.remove(dest_file)
          message(sprintf("Deleted gz file: %s", dest_file))
        }
      }
    }
  }

  db_files <- list(
    secretome_location = file.path(db_loc, "proteinatlas.tsv"),
    specificity_blood = file.path(db_loc, "proteinatlas.tsv"),
    specificity_brain = file.path(db_loc, "proteinatlas.tsv"),
    specificity_tissue = file.path(db_loc, "proteinatlas.tsv"),
    specificity_celline = file.path(db_loc, "proteinatlas.tsv"),
    specificity_singlecell = file.path(db_loc, "proteinatlas.tsv"),
    subcellular_location = file.path(db_loc, "proteinatlas.tsv"),
    protein_class = file.path(db_loc, "proteinatlas.tsv"),
    panglao_cellmarkers = file.path(
      db_loc,
      "PanglaoDB_markers_27_Mar_2020.tsv"
    ),
    reactome = file.path(db_loc, "Ensembl2Reactome_All_Levels.txt"),
    trrust = file.path(db_loc, "trrust_rawdata.human.tsv")
  )

  return(db_files)
}

#' Load and format an annotation database
#'
#' @param db_id The ID of the database to load (e.g., "reactome").
#' @param db_files Named list mapping database IDs to file paths.
#' @param AnnDatR AnnDatR object containing clustering results.
#'
#' @return A tibble with standardized columns: `ensg_id`, `term`, `term_id`.
#' @keywords internal
load_annotation_database <- function(db_id, db_files, AnnDatR) {
  cluster_genes <- AnnDatR$uns$consensus_clustering$gene
  if (db_id == "reactome") {
    raw_db <- readr::read_tsv(
      db_files[[db_id]],
      col_names = c(
        "ensg_id",
        "id",
        "URL",
        "description",
        "evidence_code",
        "species"
      ),
      show_col_types = FALSE
    )
  } else {
    raw_db <- readr::read_tsv(db_files[[db_id]], show_col_types = FALSE)
  }
  formatted_db <- format_annotation_database(
    raw_db,
    db_id,
    db_files,
    cluster_genes
  )
  stopifnot(all(names(formatted_db) == c("ensg_id", "term", "term_id")))

  return(formatted_db)
}

#' Format an annotation database (excluding GO)
#'
#' @param raw_db A tibble containing the raw database data.
#' @param db_id The ID of the database to format (e.g., "secretome_location").
#' @param db_files Named list mapping database IDs to file paths.
#' @param cluster_genes Character vector of Ensembl IDs to keep.
#'
#' @return A tibble with standardized columns: `ensg_id`, `term`, `term_id`.
#' @keywords internal
format_annotation_database <- function(raw_db, db_id, db_files, cluster_genes) {
  db_format_function <- database_formatting_functions[[db_id]]
  if (is.null(db_format_function)) {
    stop(sprintf("No formatting function found for database ID: %s", db_id))
  }
  formatted_db <- db_format_function(raw_db, db_files, cluster_genes)
  stopifnot(all(names(formatted_db) == c("ensg_id", "term", "term_id")))

  return(formatted_db)
}

#' Formatting functions for annotation databases (excluding GO)
#'
#' Each function takes (db, db_files, cluster_genes) as arguments.
#' @keywords internal
database_formatting_functions <- list(
  secretome_location = function(db, db_files, cluster_genes) {
    db |>
      dplyr::select(ensg_id = Ensembl, term = `Secretome location`) |>
      dplyr::filter(!is.na(term)) |>
      dplyr::mutate(term_id = term) |>
      dplyr::filter(ensg_id %in% cluster_genes)
  },
  specificity_blood = function(db, db_files, cluster_genes) {
    db |>
      dplyr::select(ensg_id = Ensembl, term = `RNA blood cell specificity`) |>
      dplyr::mutate(term = gsub(", ", ";", term)) |>
      tidyr::separate_rows(term, sep = ";") |>
      dplyr::filter(!is.na(term)) |>
      dplyr::mutate(term_id = term) |>
      dplyr::filter(ensg_id %in% cluster_genes)
  },
  specificity_brain = function(db, db_files, cluster_genes) {
    db |>
      dplyr::select(
        ensg_id = Ensembl,
        term = `RNA single nuclei brain specificity`
      ) |>
      dplyr::mutate(term = gsub(", ", ";", term)) |>
      tidyr::separate_rows(term, sep = ";") |>
      dplyr::filter(!is.na(term)) |>
      dplyr::mutate(term_id = term) |>
      dplyr::filter(ensg_id %in% cluster_genes)
  },
  specificity_tissue = function(db, db_files, cluster_genes) {
    db |>
      dplyr::select(ensg_id = Ensembl, term = `RNA tissue specificity`) |>
      dplyr::mutate(term = gsub(", ", ";", term)) |>
      tidyr::separate_rows(term, sep = ";") |>
      dplyr::filter(!is.na(term)) |>
      dplyr::mutate(term_id = term) |>
      dplyr::filter(ensg_id %in% cluster_genes)
  },
  specificity_celline = function(db, db_files, cluster_genes) {
    db |>
      dplyr::select(ensg_id = Ensembl, term = `RNA cell line specificity`) |>
      dplyr::mutate(term = gsub(", ", ";", term)) |>
      tidyr::separate_rows(term, sep = ";") |>
      dplyr::filter(!is.na(term)) |>
      dplyr::mutate(term_id = term) |>
      dplyr::filter(ensg_id %in% cluster_genes)
  },
  specificity_singlecell = function(db, db_files, cluster_genes) {
    db |>
      dplyr::select(
        ensg_id = Ensembl,
        term = `RNA single cell type specificity`
      ) |>
      dplyr::mutate(term = gsub(", ", ";", term)) |>
      tidyr::separate_rows(term, sep = ";") |>
      dplyr::filter(!is.na(term)) |>
      dplyr::mutate(term_id = term) |>
      dplyr::filter(ensg_id %in% cluster_genes)
  },
  subcellular_location = function(db, db_files, cluster_genes) {
    db |>
      dplyr::select(ensg_id = Ensembl, term = `Subcellular location`) |>
      dplyr::mutate(term = gsub(", ", ";", term)) |>
      tidyr::separate_rows(term, sep = ";") |>
      dplyr::filter(!is.na(term)) |>
      dplyr::mutate(term_id = term) |>
      dplyr::filter(ensg_id %in% cluster_genes)
  },
  protein_class = function(db, db_files, cluster_genes) {
    db |>
      dplyr::select(ensg_id = Ensembl, protein_class = `Protein class`) |>
      tidyr::separate_rows(protein_class, sep = ", ") |>
      dplyr::rename(term = protein_class) |>
      dplyr::filter(!is.na(term)) |>
      dplyr::mutate(term_id = term) |>
      dplyr::filter(ensg_id %in% cluster_genes)
  },
  panglao_cellmarkers = function(db, db_files, cluster_genes) {
    hpa_file <- db_files[["secretome_location"]]
    hpa <- readr::read_tsv(hpa_file, show_col_types = FALSE)
    hpa_map <- hpa |>
      dplyr::select(gene_name = Gene, ensg_id = Ensembl) |>
      dplyr::distinct()
    db |>
      dplyr::filter(
        grepl("Hs", species),
        `gene type` %in% c("protein-coding gene", "protein coding gene")
      ) |>
      dplyr::select(gene_name = 2, term = 3) |>
      dplyr::distinct() |>
      dplyr::inner_join(hpa_map, by = "gene_name") |>
      dplyr::filter(!is.na(ensg_id)) |>
      dplyr::mutate(term_id = term) |>
      dplyr::select(ensg_id, term, term_id) |>
      dplyr::filter(ensg_id %in% cluster_genes)
  },
  reactome = function(db, db_files, cluster_genes) {
    db |>
      dplyr::filter(species == "Homo sapiens") |>
      dplyr::select(ensg_id, term = description, term_id = id) |>
      dplyr::distinct() |>
      dplyr::filter(ensg_id %in% cluster_genes)
  },
  trrust = function(db, db_files, cluster_genes) {
    hpa_file <- db_files[["secretome_location"]]
    hpa <- readr::read_tsv(hpa_file, show_col_types = FALSE)
    hpa_map <- hpa |>
      dplyr::select(gene_name = Gene, ensg_id = Ensembl) |>
      dplyr::distinct()
    db |>
      dplyr::select(gene_name = BAX, term = AATF) |>
      dplyr::inner_join(hpa_map, by = "gene_name") |>
      dplyr::filter(!is.na(ensg_id)) |>
      dplyr::mutate(term_id = term) |>
      dplyr::select(ensg_id, term, term_id) |>
      dplyr::filter(ensg_id %in% cluster_genes)
  }
)

#' Build annotation terms tibble from all databases in db_files
#'
#' @param db_files Named list mapping database IDs to file paths (as returned by get_annot_dbs()).
#' @param AnnDatR AnnDatR object containing clustering results.
#'
#' @return A tibble with columns: db_id, ensg_id, term, term_id.
#' @keywords internal
build_annotation_terms_tibble <- function(db_files, AnnDatR) {
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop(
      "The 'readr' package is required for this function. Please install it using install.packages('readr')."
    )
  }
  if (is.null(AnnDatR[["uns"]][["consensus_clustering"]])) {
    stop(
      "AnnDatR$uns$consensus_clustering not found. Call `hc_cluster_consensus()` before `hc_annotate()`."
    )
  }
  cluster_genes <- AnnDatR[["uns"]][["consensus_clustering"]][["gene"]]
  db_ids <- names(db_files)
  res_list <- lapply(db_ids, function(db_id) {
    if (db_id == "reactome") {
      raw_db <- readr::read_tsv(
        db_files[[db_id]],
        col_names = c(
          "ensg_id",
          "id",
          "URL",
          "description",
          "evidence_code",
          "species"
        ),
        show_col_types = FALSE
      )
    } else {
      raw_db <- readr::read_tsv(db_files[[db_id]], show_col_types = FALSE)
    }
    formatted_db <- format_annotation_database(
      raw_db,
      db_id,
      db_files,
      cluster_genes
    )
    formatted_db[["db_id"]] <- db_id
    formatted_db
  })

  res_tibble <- dplyr::bind_rows(res_list) |>
    dplyr::select(
      dplyr::any_of(c("db_id", "ensg_id", "term", "term_id")),
      dplyr::everything()
    )

  return(res_tibble)
}
