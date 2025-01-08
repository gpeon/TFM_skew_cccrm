plot_see <- function(x, size_point = 2, panel = TRUE, ...) {
  model <- see:::.retrieve_data(x)
  x <- x[-which(x$p_Residuals == 0 & x$p_Response == 0), ]
  
  dat <- data.frame(
    x = factor(c(x$Distribution, x$Distribution), levels = rev(sort(unique(x$Distribution)))),
    y = c(x$p_Response, x$p_Residuals),
    group = factor(
      c(rep("Response", length(x$p_Response)), rep("Residuals", length(x$p_Residuals))),
      levels = c("Response", "Residuals")
    ),
    stringsAsFactors = FALSE
  )
  
  # remove all zero-probabilities
  dat$y[dat$y == 0] <- NA
  
  # find max-range, add some tolerance
  max_y <- max(dat$y, na.rm = TRUE) * 1.1
  
  # default legend-position
  lp <- ifelse(isTRUE(panel), "right", "bottom")
  
  p1 <- ggplot(
    dat,
    aes(
      y = .data$x,
      x = .data$y,
      colour = .data$group
    )
  ) +
    geom_linerange(
      aes(xmin = 0, xmax = .data$y),
      position = position_dodge(0.4),
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    geom_point(
      size = size_point,
      position = position_dodge(0.4),
      na.rm = TRUE
    ) +
    labs(
      y = NULL,
      x = NULL,
      fill = NULL,
      colour = NULL,
      title = "Predicted Distribution of Residuals and Response"
    ) +
    scale_x_continuous(
      labels = see:::.percents,
      expand = c(0, 0),
      limits = c(0, max_y)
    ) +
    see::scale_color_material_d(reverse = TRUE) +
    guides(colour = guide_legend(reverse = TRUE)) +
    see::theme_lucid(legend.position = lp)
  
  dat1 <- as.data.frame(stats::density(stats::residuals(model)))
  dat2 <- data.frame(x = see:::.factor_to_numeric(insight::get_response(model)))
  
  # this ensures that we have integers for response variables with values
  # much greater than 1 and with a large range, so we have proper bar plots
  if (diff(range(dat2$x)) > 10) dat2$x <- round(dat2$x)
  
  p2 <- ggplot(dat1, aes(x = .data$x, y = .data$y)) +
    geom_line(colour = "#2196F3") +
    labs(x = NULL, y = NULL, title = "Density of Residuals") +
    see::theme_lucid()
  
  # usually, we have an outline for the bars in the same color as the
  # theme background, however, for very thin bars, this results in "invisible"
  # bars, because the fill colour is not visible. For larger ranges, use fill
  # colour as color, too.
  
  if (abs(diff(range(dat2$x))) > 100) {
    bar_color <- "#f44336"
  } else {
    bar_color <- see::theme_lucid()$panel.background$fill
  }
  
  p3 <- ggplot(dat2, aes(x = .data$x)) +
    geom_histogram(
      fill = "#f44336", colour = bar_color,
      na.rm = TRUE
    ) +
    labs(x = NULL, y = NULL, title = "Distribution of Response") +
    see::theme_lucid()
  
  if (panel) {
    insight::check_if_installed("patchwork")
    return(p1 / (p2 | p3) + patchwork::plot_layout(nrow = 2))
  } else {
    return(list(p1, p2, p3))
  }
}