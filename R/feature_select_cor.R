feature_select_cor <- function(exp,sig, r = 0.6){
  exp <- pp_exp(exp, sig)
  sig <- pp_sig(exp, sig)
  fg_all <- NULL

  for (i in levels(sig$cell)){
    sig_subb <- sig[sig$cell %in% i,]
    exp_subb <- exp[as.character(sig_subb$gene),]
    cor_data <- corr_matrix(exp_subb, sig, i, r)
    fg_subb <- apply(cor_data,1,function(x)!all(x==0))
    fg_subb <- tryCatch(
      {
        data.frame(row.names(cor_data[fg_subb,fg_subb]), i)
      },
      error = function(cond) {
        data.frame(sig_subb$gene, i)
      })

    fg_all <- rbindlist(list(fg_all, fg_subb), use.names = FALSE)
  }
  colnames(fg_all) <- c("gene", "cell")
  fg_all <- data.frame(lapply(fg_all, as.character), stringsAsFactors = F)
  return(fg_all)
}