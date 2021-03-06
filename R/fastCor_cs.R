#' @title Pre-processing expression matrix
#' @description Subsets the user's dataset based on the genes that are common to the users dataset and Signature genes.
#' @param exp Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.
#' @return Expression dataframe
#' @export

fastCor_cs <- function (xt, nSplit = 1, upperTri = FALSE, optBLAS = FALSE,
          verbose = TRUE)
{
  varnames <- colnames(xt)
  nn <- ncol(xt)
  ii <- which(apply(xt, 2, var) > 0)
  if (verbose && length(ii) < nn) {
    write("---> Checking zero-variance data...", "")
    write(paste("--->\t Total number of variables: ", nn),
          "")
    write(paste("--->\t WARNING:", nn - length(nn[ii]),
                "variables found with zero variance"), "")
  }
  xt <- xt[, ii]
  x <- t(xt) - colMeans(xt)
  m <- nrow(xt)
  n <- ncol(xt)
  r <- matrix(NA, nrow = n, ncol = n)
  nSplitMax <- floor(n/2)
  if (nSplit > nSplitMax) {
    if (verbose) {
      write(paste("---> Maximum number of splits: floor(n/2) =",
                  nSplitMax), "")
      write(paste("---> WARNING: number of splits nSplit >",
                  nSplitMax), "")
      write(paste("---> WARNING: using maximum number of splits: nSplit =",
                  nSplitMax), "")
    }
    nSplit <- nSplitMax
  }
  if (nSplit == 1) {
    if (optBLAS & .Machine$sizeof.pointer == 8) {
      r <- tcrossprod(x/sqrt(rowSums(x^2)))
    }
    else {
      r <- cor_cs(xt)
    }
  }
  else if (nSplit > 1) {
    if (verbose)
      write("---> Computing split sizes...", "")
    lSplit <- floor(n/nSplit)
    iSplit <- vector("list", nSplit)
    for (i in 1:(nSplit - 1)) {
      iSplit[[i]] <- (lSplit * (i - 1) + 1):(lSplit *
                                               i)
    }
    iSplit[[nSplit]] <- (lSplit * (nSplit - 1) + 1):n
    if (verbose)
      write("---> Computing split combinations...", "")
    cSplit <- cbind(combn(nSplit, 2), rbind(c(1:nSplit),
                                            c(1:nSplit)))
    cSplit <- cSplit[, order(cSplit[1, ], cSplit[2, ])]
    if (verbose && n%%nSplit == 0) {
      write(paste("---> Splitting data matrix:", nSplit,
                  "splits of", paste(m, "x", length(iSplit[[1]]),
                                     sep = ""), "size"), "")
    }
    else {
      write(paste("---> Splitting data matrix:", (nSplit -
                                                    1), ifelse((nSplit - 1) == 1, "split of", "splits of"),
                  paste(m, "x", length(iSplit[[1]]), sep = ""),
                  "size"), "")
      write(paste("---> Splitting data matrix:", 1, "split of",
                  paste(m, "x", length(iSplit[[nSplit]]), sep = ""),
                  "size"), "")
    }
    for (nc in 1:max(nSplit, ncol(cSplit))) {
      i1 <- iSplit[[cSplit[1, nc]]]
      i2 <- iSplit[[cSplit[2, nc]]]
      if (verbose)
        write(paste("---> Correlation matrix: split #",
                    sprintf("%3.0f", cSplit[1, nc]), "   x   split #",
                    sprintf("%3.0f", cSplit[2, nc]), sep = ""),
              "")
      if (optBLAS & .Machine$sizeof.pointer == 8) {
        r[i2, i1] <- t(tcrossprod(x[i1, ]/sqrt(rowSums(x[i1,
                                                         ]^2)), x[i2, ]/sqrt(rowSums(x[i2, ]^2))))
      }
      else {
        r[i2, i1] <- t(cor(xt[, i1], xt[, i2]))
      }
      if (!upperTri)
        r[i1, i2] <- t(r[i2, i1])
    }
  }
  else {
    stop(paste("invalid nSplit:", nSplit))
  }
  if (upperTri)
    r[col(r) >= row(r)] <- NA
  rr <- matrix(NA, nrow = nn, ncol = nn)
  rr[ii, ii] <- r
  rownames(rr) <- varnames
  colnames(rr) <- varnames
  return(rr)
}