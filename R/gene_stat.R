#' @title General stastitics of ImSig analysis
#' @description [Total genes in ImSig]: The total number of genes in ImSig list. [No. of ImSig genes in user dataset]: The number of ImSig genes found in user's dataset. Like all signatures, ImSig works best when this overlap is high, preferably over 75%. [Median correlation of ImSig genes]: This is the number of remaining ImSig genes after \code{\link{feature_select}}. If this number drops drastically or if the median correlation is low, it may indicate the absence of the particular cell type in users dataset. [Median correlation of feature selected ImSig genes]: These values again can be used gauge the presence or absence of a cell type. As ImSig genes were designed to be co-expressed when the cell type is present, poor correlations may indicate absence of the cell type in the dataset. A network graph (\code{\link{plot_network }}) may be generated to visually confirm the observation.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @param r Use a value between 0 and 1. Default is 0.6. This is a user defined correlation cut-off to perform feature selection (\code{\link{feature_select}}). Feature selection aids to enrich the prediction of relative abundance of immune cells by filtering off poorly correlated ImSig genes. To get an idea of what cut-off to use check the results of (\code{\link{gene_stat}}) and choose a cut-off that displays high median correlation and maintains a high proportion of genes after feature selection.
#' @return Dataframe of general statistics of ImSig analysis.
#' @import stats
#' @importFrom HiClimR fastCor
#' @examples
#' gene_stat (exp = example_data, r = 0.7)
#' @seealso \code{\link{feature_select}}
#' @export

gene_stat <- function(exp, sig = sig, r = 0.6){
  sig_total <- sig
  exp <- pp_exp(exp, sig)
  sig <- pp_sig(exp, sig)
  fg <- feature_select(exp, sig, r)
  e_cor_data <- fastCor(t(exp))
  diag(e_cor_data) = NA
  genestat <- data.frame(matrix(nrow = 5))
  genestat <- genestat[,-1]
  g <- row.names(exp)
  for (i in levels(sig$cell)){
    sig_total_sub <- sig_total[sig_total$cell %in% i,]
    sig_len <- length(sig_total_sub$gene)

    exp_sub <- sig[sig$cell %in% i,]
    user_len <- length(exp_sub$gene)

    feature_sub <- exp_sub[as.character(exp_sub$gene) %in% fg,]
    feature_len <- length(feature_sub$gene)

    sig_cor <- e_cor_data[as.character(exp_sub$gene), as.character(exp_sub$gene)]

    s_cor <- tryCatch(
      {
        round(median(apply(sig_cor, 1, median, na.rm = T)), 2)
      },
      error = function(cond) {
        return(NA)
      })

    feature_cor <- e_cor_data[as.character(feature_sub$gene), as.character(feature_sub$gene)]
    f_cor <- round(median(apply(feature_cor, 1, median, na.rm = T)), 2)

    stat <- data.frame(c(sig_len, user_len, feature_len, s_cor, f_cor))
    genestat <- cbind(genestat,stat)
  }
  row.names(genestat) <- c("Total signature genes", "No. of Signature genes in user dataset", "No. of signature genes after feature selection", "Median correlation of signature genes", "Median correlation of feature selected signature genes")
  names(genestat) <- levels(sig$cell)
  return(genestat)
}