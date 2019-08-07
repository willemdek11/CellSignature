#' @title Pre-processing signature genes to return cell types present.
#' @description List of cell types present in the Signature
#' @param exp List of cell types present in Signature table. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @return Signature Cell types list

pp_cell_type <- function (sig)
{
  cell_type <- levels(sig$cell)
  cell_type <- sub("-", " ", cell_type, perl = TRUE)
  cell_type <- unlist(strsplit(cell_type, ";"))
  cell_type <- sub(".*_", "", cell_type, perl = TRUE)
  cell_type <- unique(cell_type)
  return(cell_type)
}