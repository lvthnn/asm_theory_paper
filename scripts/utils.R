library(ggplot2)
library(ggthemes)

set_knoop_theme <- function() {
  theme_set(
    theme_par(base_size = 10, base_family = "Times New Roman") +
      theme(plot.margin = margin(0, 5.5, 0, 5.5)) +
      theme(legend.position = "top") +
      theme(
        plot.title = element_text(
          size = 10,
          face = "bold",
          hjust = 0.5
        ),
        legend.title = element_text(size = 8),
        legend.text = element_text(
          size = 8,
          margin = margin(l = 3, r = 0)
        ),
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(10, "pt"),
        legend.key.spacing.x = unit(3, "pt"),
        panel.grid.major.y = element_line(
          colour = "grey92",
          linewidth = 0.3
        ),
        panel.grid.major.x = element_line(
          colour = "grey92",
          linewidth = 0.3
        ),
        panel.grid.minor = element_blank()
      )
  )

  invisible(NULL)
}

vech <- function(mat) {
  mat[lower.tri(mat, diag = TRUE)]
}

unvech <- function(vec, dim) {
  mat <- matrix(0, dim, dim)
  mat[lower.tri(mat, diag = TRUE)] <- vec
  mat + t(mat) - diag(diag(mat))
}
