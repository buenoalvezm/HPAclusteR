#' A simple theme for visualizations
#'
#' @return A ggplot2 theme object that can be added to plots with `+`.
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mpg, aes(displ, hwy)) +
#'   geom_point() +
#'   theme_cluster()
#'
theme_cluster <-
  function() {
    theme(
      panel.background = element_rect(fill = NA, colour = NA),
      plot.background = element_rect(fill = NA, color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black", size = 0.5)
    )
  }
