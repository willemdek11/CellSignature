#' @title Pre-processing correlation matrix
#' @description Subsets the correlation dataset based on the genes that are common to Cell Signature and total correlation matrix.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @return Correlation dataframe
#' @export

pp_cor <- function (cor_m, sig, cell_type)
{
  sig <- sig
  g <- Reduce(intersect, list(as.character(row.names(cor_m)), sig[sig$cell %in% cell_type,]$gene))
  cor_m <- cor_m[as.character(g), as.character(g), drop=FALSE]
  return(cor_m)
}