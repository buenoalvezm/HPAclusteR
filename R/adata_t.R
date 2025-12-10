#' Human Protein Atlas Tissue Data (transposed)
#'
#' Human Protein Atlas tissue bulk-RNA expression data in AnnDatR format transposed.
#'
#'
#' @format
#' An AnnDatR object with the following structure:
#' \describe{
#'   \item{X}{A data frame with 78 rows and 20,163 columns representing the main data matrix (e.g., gene expression values).}
#'   \item{obs}{A data frame with 78 rows and 20,162 columns containing metadata for observations (e.g., genes).}
#'   \item{obs_names_col}{A string specifying the column used for observation names (`ensembl_id`).}
#'   \item{var}{A data frame with 78 rows and 7 columns containing metadata for variables (e.g., samples).}
#'   \item{var_names_col}{A string specifying the column used for variable names (`sample_id`).}
#'   \item{uns}{A list for storing unstructured metadata (currently `NULL`).}
#'   \item{obsm}{A list for storing multidimensional annotations for observations (currently `NULL`).}
#'   \item{raw}{A slot for storing raw data (currently `NULL`).}
#' }
#' @source <https://www.science.org/doi/10.1126/science.1260419>
#'
#' @examples
#' adata_t
"adata_t"
