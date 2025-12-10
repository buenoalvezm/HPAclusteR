#' Scale data matrix
#'
#' @param data Data to be scaled
#' @param method Scaling method: 'sample' (default), 'gene', 'pareto', 'vector', 'min-max', 'max_sequencial', or 'max'
#'
#' @returns Scaled data matrix
#' @keywords internal
scale_data <- function(data, method = 'sample') {
  if (method == 'gene') {
    scaled_data <- data |>
      t() |>
      pcaMethods::prep(scale = 'uv', eps = .Machine[["double.eps"]])
  } else if (method == 'sample') {
    scaled_data <- data |>
      pcaMethods::prep(scale = 'uv', eps = .Machine[["double.eps"]]) |>
      t()
  } else if (method == 'pareto') {
    scaled_data <- data |> t() |> pcaMethods::prep(scale = 'pareto')
  } else if (method == 'vector') {
    scaled_data <- data |> t() |> pcaMethods::prep(scale = 'vector')
  } else if (method == 'min-max') {
    scaled_data <- apply(t(data), 1, function(x) {
      (x - min(x)) / (max(x) - min(x))
    }) |>
      t()
  } else if (method == 'max_sequencial') {
    scaled_data <- apply(t(data), 1, function(x) x / max(x)) |> t()
  } else if (method == 'max') {
    data <- t(data)
    scaled_data <- data / matrixStats::rowMaxs(data)
  } else {
    stop(
      'Scaling method not defined, only "sample", "gene", "pareto", "vector", "min-max", max_sequencial, and "max" are allowed'
    )
  }

  return(scaled_data)
}
