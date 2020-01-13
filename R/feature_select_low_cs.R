#' @title Feature selection of signature genes
#' @description Signature genes are expected to be co-expressed in tissue transcriptomic data. However, depending on the dataset some of the genes may not co-express with the dominant module. In order to remove such deviant genes, a feature selection can be carried out based on correlation. This function removes genes that exhibit a poor correlation (less than the defined r value) with in the experimental data. This step of feature selection is recommended to enrich the prediction of relative abundance of immune cells.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @param sig Dataframe of Cell Signature containing genes and their cell signature in rows.
#' @param r Use a value between 0 and 1. Default is 0.6. This is a user defined correlation cut-off to perform feature selection. To get an idea of what cut-off to use check the results of (\code{\link{gene_stat}}) and choose a cut-off that displays high median correlation and maintains a high proportion of genes after feature selection.
#' @return Returns a table of 'feature selected' genes based on the set r value.
#' @examples
#' feature_select_cs (exp = example_data, sig = signature_table, r = 0.7)
#' @export

feature_select_low_cs <- function(exp,sig, r = 0.55, l = 0.2){
  exp <- pp_exp(exp, sig)
  sig <- pp_sig(exp, sig)
  fg_all <- NULL
  cor_data <- fastCor_cs(t(exp))

  for (i in levels(sig$cell)){
    HIGH <- TRUE
    cor_subb <- pp_cor(cor_data, sig, i)

    npass = c()
    for(j in 1:dim(cor_subb)[1])
    {
      npass[j] = sum(cor_subb[j,]> r,na.rm=T)
    }

    if(sum(npass) == dim(cor_subb)[1]){
      npass = c()
      for(j in 1:dim(cor_subb)[1])
      {
        npass[j] = sum(cor_subb[j,]> l,na.rm=T)
      }
      i <- paste0(i , "*")
      HIGH <- FALSE
    }

    to.keep = (1:dim(cor_subb)[1])[(npass == max(npass))&(npass>0)]
    # only keep the one with the highest total cor:
    if(length(to.keep)>1){
      to.keep = (to.keep)[rowSums(cor_subb[to.keep,,drop=FALSE],na.r=T)==max(rowSums(cor_subb[to.keep,,drop=FALSE],na.rm=T))]
      }

    if(length(to.keep)>0 & HIGH)
    {
      genes = row.names(cor_subb[,to.keep,drop=FALSE])[cor_subb[,to.keep] > r]
    } else
      genes = row.names(cor_subb[,to.keep,drop=FALSE])[cor_subb[,to.keep] > l]

    fg_subb <- tryCatch(
      {
        data.frame(genes, i)
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