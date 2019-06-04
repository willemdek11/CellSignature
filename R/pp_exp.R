pp_exp <- function (exp)
{
  sig <- sig
  g <- Reduce(intersect, list(as.character(row.names(exp),
                                           sig$gene)))
  exp <- exp[as.character(g), ]
  return(exp)
}
