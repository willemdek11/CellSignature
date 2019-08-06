pp_cell_type <- function (sig)
{
  cell_type <- levels(sig$cell)
  cell_type <- sub("-", " ", cell_type, perl = TRUE)
  cell_type <- unlist(strsplit(cell_type, ";"))
  cell_type <- sub(".*_", "", cell_type, perl = TRUE)
  cell_type <- unique(cell_type)
  return(cell_type)
}