write_as_xtable <- function(table1, filepath, type = 'latex', coerce_to_df = TRUE, ...) {
  if(coerce_to_df) {
    table1 <- as.data.frame(table1)
  } 

  rownames(table1) <- NULL

  xtable::xtable(table1, ...) |>
    xtable::print.xtable(
      type = type, 
      file = filepath,
      include.rownames = FALSE,
      comment = FALSE,
      table.placement = '!htbp',
      format.args = list(big.mark = ',')
    )

  return(filepath)
}