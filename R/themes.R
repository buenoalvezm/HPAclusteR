#' A simple theme for visualizations
#'
#' `theme_hc()` provides a minimalist ggplot2 theme with no background or grid lines. Usefull for clean, publication-ready point, bar, line and boxplots among others.
#'
#' @return A ggplot2 theme object that can be added to plots with `+`.
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot2::ggplot(mpg, aes(displ, hwy)) +
#'   ggplot2::geom_point() +
#'   theme_hc()
theme_hc <- function() {
  ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = NA, colour = NA),
    plot.background = ggplot2::element_rect(fill = NA, colour = NA),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_blank(),
    axis.line = ggplot2::element_line(colour = "black", linewidth = 0.5)
  )
}
