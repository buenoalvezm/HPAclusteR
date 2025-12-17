#' Download annotation databases for HPAclusteR
#'
#' `hc_annot_dbs()` downloads several annotation databases commonly used in single-cell RNA-seq analysis,
#' including the Human Protein Atlas, PanglaoDB cell markers, Reactome pathways, and TRRUST transcriptional regulatory relationships.
#' The databases are downloaded to the specified destination directory and unzipped or decompressed if necessary.
#'
#' @param dest_dir Destination directory to save the downloaded databases (default: "databases").
#' @param hpa_version Version of the Human Protein Atlas to download (default: 24).
#'
#' @returns Invisible TRUE if all downloads are successful.
#'
#' @export
#' @examples
#' hc_annot_dbs()
hc_annot_dbs <- function(
  dest_dir = "databases",
  hpa_version = 24
) {
  if (!requireNamespace("utils", quietly = TRUE)) {
    stop(
      "The 'utils' package is required for this function. Please install it using install.packages('utils')."
    )
  }

  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

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
    dest_file <- file.path(dest_dir, basename(url))
    if (!file.exists(dest_file)) {
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
              utils::unzip(dest_file, exdir = dest_dir)
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
    } else {
      message(sprintf("%s already exists, skipping.", dest_file))
    }
  }
  invisible(TRUE)
}
