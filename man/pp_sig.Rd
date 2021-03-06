% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pp_sig.R
\name{pp_sig}
\alias{pp_sig}
\title{Pre-processing signature genes file}
\usage{
pp_sig(exp, sig)
}
\arguments{
\item{exp}{Dataframe of transcriptomic data (natural scale) containing genes as rows and samples as columns. Note: Gene names should be set as row names and duplicates are not allowed. Missing values are not allowed within the expression matrix. Check example- head(example_data): \code{\link{example_data}}.}
}
\value{
Signature genes dataframe
}
\description{
Subsets signature genes based on the genes that are common to the users dataset and signature genes
}
