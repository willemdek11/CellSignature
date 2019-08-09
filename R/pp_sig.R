#' @title Pre-processing signature genes file
#' @description Subsets signature genes based on the genes that are common to the users dataset and signature genes
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @return Signature genes dataframe
#' @export

pp_sig <- function (exp, sig)
{
  g <- Reduce(intersect, list(as.character(row.names(exp)), sig$gene))
  sig <- sig[sig$gene %in% as.character(g), ]
  output_sig <- data.frame(gene = character(), cell = character())
  for (gene in g) {
    for (cell in strsplit(as.character(sig[sig$gene == gene,]$cell), ";")) {
      cell_x <- sub("-", " ", cell, perl = TRUE)
      cell_x <- sub(".*_", "", cell_x, perl = TRUE)
      output_sig <- rbind(output_sig, data.frame(gene = gene, cell = cell_x))
    }
  }
  output_sig <- unique(output_sig)

  output_sig <- mutate_if(output_sig, is.character, as.factor)
  return(output_sig)
}