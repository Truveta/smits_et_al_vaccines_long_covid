# useful function for summarizing cox models
cox_table <- function(cox_fit) {
  cox_fit_tidy <- 
    dplyr::inner_join(
      broom::tidy(cox_fit, conf.int = TRUE),
      broom::tidy(cox_fit, conf.int = TRUE, exponentiate = TRUE), 
      by = 'term', 
      suffix = c('_base', '_exp')
    ) |>
    janitor::clean_names()

  cox_fit_tidy
}

write_hazard_table <- function(df, name_tag, output_dir, caption) {
  temp <- 
    df |>
    dplyr::select(comparison, ratio, conf.low, conf.high) |>
    dplyr::mutate(
      dplyr::across(ratio:conf.high, ~ round(.x, 2)),
      hr = paste0(ratio, ' [', conf.low, ', ', conf.high, ']')
    ) |>
    dplyr::select(comparison, hr) |>
    setNames(c('Comparison', 'Hazard Ratio [95% CI]'))

  csv_filename <- paste0(name_tag, '.csv')
  write.csv(temp, file = file.path(output_dir, csv_filename))

  label <- paste0('tab:', name_tag)
  tex_filename <- paste0(name_tag, '.tex')
  
  xtable::print.xtable(
    xtable::xtable(
      as.data.frame(temp, make.names = FALSE),
      label = label,
      caption = caption
    ),
    type = 'latex', 
    file = file.path(output_dir, tex_filename), 
    include.rownames = FALSE,
    comment = FALSE,
    timestamp = NULL,
    table.placement = '!htbp'
  )
  
  print(temp)
}
