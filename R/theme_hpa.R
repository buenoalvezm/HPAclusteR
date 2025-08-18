#' A simple theme for visualizations
#'
#' @return A ggplot2 theme object that can be added to plots with `+`.
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot2::ggplot(mpg, aes(displ, hwy)) +
#'   ggplot2::geom_point() +
#'   theme_cluster()
theme_cluster <- function() {
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = NA, colour = NA),
    plot.background = ggplot2::element_rect(fill = NA, colour = NA),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black", size = 0.5)
  )
}
