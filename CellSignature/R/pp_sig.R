pp_sig <- function (exp)
{
  g <- Reduce(intersect, list(as.character(row.names(exp),
                                           sig$gene)))
  sig <- sig[sig$gene %in% as.character(g), ]
  return(sig)
}
