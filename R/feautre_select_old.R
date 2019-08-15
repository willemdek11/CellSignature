#' @title Feature selection of signature genes
#' @description Signature genes are expected to be co-expressed in tissue transcriptomic data. However, depending on the dataset some of the genes may not co-express with the dominant module. In order to remove such deviant genes, a feature selection can be carried out based on correlation. This function removes genes that exhibit a poor correlation (less than the defined r value) with in the experimental data. This step of feature selection is recommended to enrich the prediction of relative abundance of immune cells.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @param sig Dataframe of Cell Signature containing genes and their cell signature in rows.
#' @param r Use a value between 0 and 1. Default is 0.6. This is a user defined correlation cut-off to perform feature selection. To get an idea of what cut-off to use check the results of (\code{\link{gene_stat}}) and choose a cut-off that displays high median correlation and maintains a high proportion of genes after feature selection.
#' @return Returns a table of 'feature selected' genes based on the set r value.
#' @examples
#' feature_select_cs (exp = example_data, sig = signature_table, r = 0.7)
#' @export

feature_select_old <- function(exp,sig, r = 0.7){
  exp <- pp_exp(exp, sig)
  sig <- pp_sig(exp, sig)
  fg_all <- NULL
  
  for (i in levels(sig$cell)){
    sig_subb <- sig[sig$cell %in% i,]
    exp_subb <- exp[as.character(sig_subb$gene),]
    cor_data <- corr_matrix_cs(exp_subb, sig, i, r)
    fg_subb <- apply(cor_data,1,function(x)!all(x==0))
    fg_subb <- tryCatch(
      {
        data.frame(row.names(cor_data[fg_subb,fg_subb]), i)
      },
      error = function(cond) {
        data.frame(sig_subb$gene, i)
      })
    
    try(fg_all <- rbindlist(list(fg_all, fg_subb), use.names = FALSE), silent = TRUE)
  }
  colnames(fg_all) <- c("gene", "cell")
  fg_all <- data.frame(lapply(fg_all, as.character), stringsAsFactors = F)
  return(fg_all)
}