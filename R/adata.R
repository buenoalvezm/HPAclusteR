#' Human Protein Atlas Tissue Data
#'
#' Human Protein Atlas tissue bulk-RNA expression data in AnnDatR format.
#'
#'
#' @format
#' An AnnDatR object with the following structure:
#' \describe{
#'   \item{X}{A data frame with 20,162 rows and 79 columns representing the main data matrix (e.g., gene expression values).}
#'   \item{layers}{A list containing additional data matrices, such as normalized counts (`ntpm`).}
#'   \item{obs}{A data frame with 78 rows and 7 columns containing metadata for observations (e.g., tissues).}
#'   \item{obs_names_col}{A string specifying the column used for observation names (`tissue_id`).}
#'   \item{var}{A data frame with 20,162 rows and 11 columns containing metadata for variables (e.g., genes).}
#'   \item{var_names_col}{A string specifying the column used for variable names (`ensembl_id`).}
#'   \item{uns}{A list for storing unstructured metadata (currently `NULL`).}
#'   \item{obsm}{A list for storing multidimensional annotations for observations (currently `NULL`).}
#'   \item{raw}{A slot for storing raw data (currently `NULL`).}
#' }
#' @source <https://www.science.org/doi/10.1126/science.1260419>
#'
#' @examples
#' adata
"adata"
