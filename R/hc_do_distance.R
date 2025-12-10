#' Title
#'
#' @param AnnDatR
#' @param method
#' @param n_comp
#' @param id
#' @param pca_id
#'
#' @returns
#'
#' @export
#' @examples
hc_do_distance <- function(
  AnnDatR,
  method,
  n_comp = NULL,
  id = NULL,
  pca_id = 'pca'
) {
  if (is.null(AnnDatR$uns[[pca_id]])) {
    stop('AnnDatR$uns$pca not found. Call set_PCA function before plotting.')
  }
  pca_results <- AnnDatR$uns[[pca_id]]
  if (is.null(n_comp)) {
    n_comp <- which(pca_results@R2cum > 0.8)[1]
  }

  if (is.null(AnnDatR$uns$distance)) {
    AnnDatR$uns$distance = list()
  }

  distance <- pca_results@scores[, 1:n_comp] %>%
    factoextra::get_dist(method = method)

  key <- if (!is.null(id)) id else method
  AnnDatR$uns$distance[[key]] <- distance
}
