cor_CellSignature <- function(xt)
{
  n <- ncol(xt)
  m <- matrix(NA, nrow = n, ncol = n)
  colnames(m) <- colnames(xt)
  rownames(m) <- colnames(xt)

  for (i in 1:ncol(xt)) {
    x1 <- xt[,i]
    for (j in 1:ncol(xt)) {
      x2 <- xt[,j]
      sharedvar <- mean(c(var(x1),var(x2)))
      m[i,j] <- sum((x1-mean(x1))*(x2-mean(x2)))/((length(x1)-1)*sharedvar)
    }
  }
  return(m)
}