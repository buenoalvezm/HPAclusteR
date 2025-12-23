#' Initialize AnnDatR object for gene clustering (genes as samples)
#'
#' `hc_initialize()` takes three data frames: an expression matrix, gene metadata, and sample metadata.
#' The expression matrix must have genes as rows and samples as columns, with the first column as gene IDs.
#' The gene metadata (`obs`) must have a column matching the gene IDs in the expression matrix.
#' The sample metadata (`var`) must have a column matching the sample IDs (column names) in the expression matrix.
#'
#' @param X Expression matrix data frame (genes as rows, samples as columns, first column = gene IDs)
#' @param obs Gene metadata data frame (must contain a column matching gene IDs)
#' @param var Sample metadata data frame (must contain a column matching sample IDs)
#' @param obs_names_col Name of the column in gene metadata to use as gene IDs (default: "ensembl_id")
#' @param var_names_col Name of the column in sample metadata to use as sample IDs (default: "sample_id")
#'
#' @return AnnDatR object with genes as samples, ready for clustering
#' @export
hc_initialize <- function(
  X,
  obs,
  var,
  obs_names_col = "ensembl_id",
  var_names_col = "sample_id"
) {
  X_t <- X |>
    tibble::column_to_rownames(obs_names_col) |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column(var_names_col)

  var <- var[match(X_t[[var_names_col]], var[[var_names_col]]), ]
  stopifnot(all(X_t[[var_names_col]] == var[[var_names_col]]))
  obs <- obs[match(colnames(X_t)[-1], obs[[obs_names_col]]), ]
  stopifnot(all(colnames(X_t)[-1] == obs[[obs_names_col]]))

  adata_t <- AnnDatR$new(
    X = X_t,
    obs = obs,
    var = var,
    var_names = var_names_col
  )

  return(adata_t)
}
