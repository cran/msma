% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/src.r
\docType{methods}
\name{hcmsma}
\alias{hcmsma}
\title{Hierarchical cluster analysis}
\usage{
hcmsma(
  object,
  nclust = 4,
  graph = FALSE,
  hmethod = "ward.D2",
  axes = c(1, 2),
  block = "block",
  XY = "X"
)
}
\arguments{
\item{object}{an object of class "\code{msma}", usually, a result of a call to \code{\link{msma}}}

\item{nclust}{a numeric scalar, number of clusters.}

\item{graph}{a numeric vector, numbers of columns for Y. The length of vector corresponds to the number of blocks.}

\item{hmethod, }{the agglomeration method to be used in the function "\code{hclust}".}

\item{axes}{a numeric (or vector), specifying the component(s) to analyze.}

\item{block}{a character, indicating which the "block" or "super" is used.}

\item{XY}{a character, indicating "X" or "Y".}
}
\value{
\item{hcout}{An object of class hclust}

\item{clusters}{a vector with group memberships}

\item{object}{an object of class "\code{msma}", usually, a result of a call to \code{\link{msma}}}
}
\description{
This is a function for performing a hierarchical cluster analysis using scores
}
\details{
This function performs a hierarchical cluster analysis using scores.
}
