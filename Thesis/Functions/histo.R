gghistostats <- function(
    data,
    x,
    grouping.var = NULL, # New parameter for grouping variable
    binwidth = NULL,
    xlab = NULL,
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    type = "parametric",
    test.value = 0,
    bf.prior = 0.707,
    bf.message = TRUE,
    effsize.type = "g",
    conf.level = 0.95,
    tr = 0.2,
    digits = 2L,
    ggtheme = ggstatsplot::theme_ggstatsplot(),
    results.subtitle = TRUE,
    bin.args = list(color = "black", fill = "grey50", alpha = 0.7),
    centrality.plotting = TRUE,
    centrality.type = type,
    centrality.line.args = list(color = "blue", linewidth = 1, linetype = "dashed"),
    ggplot.component = NULL,
    facet = TRUE, # New parameter: use facets for grouped plots
    ...
) {
  x <- ensym(x)
  grouping.var <- enquo(grouping.var) # Enclose the grouping variable
  data <- tidyr::drop_na(select(data, {{ x }}, !!grouping.var))
  x_vec <- pull(data, {{ x }})
  type <- stats_type_switch(type)
  
  # Statistical analysis
  if (results.subtitle) {
    .f.args <- list(
      data = data,
      x = {{ x }},
      test.value = test.value,
      effsize.type = effsize.type,
      conf.level = conf.level,
      digits = digits,
      tr = tr,
      bf.prior = bf.prior
    )
    
    subtitle_df <- .eval_f(one_sample_test, !!!.f.args, type = type)
    subtitle <- .extract_expression(subtitle_df)
    
    if (type == "parametric" && bf.message) {
      caption_df <- .eval_f(one_sample_test, !!!.f.args, type = "bayes")
      caption <- .extract_expression(caption_df)
    }
  }
  
  # Plotting
  plot_hist <- ggplot(data, mapping = aes(x = {{ x }}, fill = !!grouping.var)) +
    exec(
      stat_bin,
      mapping  = aes(y = after_stat(count), fill = after_stat(count)),
      binwidth = binwidth %||% .binwidth(x_vec),
      !!!bin.args
    ) +
    scale_y_continuous(
      sec.axis = sec_axis(
        transform = ~ . / nrow(data),
        labels = function(x) insight::format_percent(x, digits = 0L),
        name = "proportion"
      )
    ) +
    labs(
      x        = xlab %||% as_name(x),
      y        = "count",
      title    = title,
      subtitle = subtitle,
      caption  = caption
    ) +
    ggtheme +
    ggplot.component
  
  if (isTRUE(facet) && !is.null(grouping.var)) {
    plot_hist <- plot_hist + facet_wrap(vars(!!grouping.var), scales = "free_y")
  } else if (!is.null(grouping.var)) {
    plot_hist <- plot_hist + guides(fill = "none")
  }
  
  # Centrality line
  if (isTRUE(centrality.plotting)) {
    plot_hist <- .histo_labeller(
      plot_hist,
      x = x_vec,
      type = stats_type_switch(centrality.type),
      tr = tr,
      digits = digits,
      centrality.line.args = centrality.line.args
    )
  }
  
  plot_hist
}

grouped_gghistostats <- function(
    data,
    x,
    grouping.var,
    binwidth = NULL,
    ...
) {
  gghistostats(
    data = data,
    x = {{ x }},
    grouping.var = {{ grouping.var }},
    binwidth = binwidth,
    ...
  )
}
