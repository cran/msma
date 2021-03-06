% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/src.r
\docType{methods}
\name{summary.msma}
\alias{summary.msma}
\alias{print.summary.msma}
\title{Summarizing Fits}
\usage{
\method{summary}{msma}(object, ...)

\method{print}{summary.msma}(x, ...)
}
\arguments{
\item{object, x}{an object of class "\code{msma}", usually, a result of a call to \code{\link{msma}}}

\item{...}{further arguments passed to or from other methods.}
}
\description{
summary method for class "msma".
}
\details{
This function provide the summary of results .
}
\examples{
##### data #####
tmpdata = simdata(n = 50, rho = 0.8, Yps = c(10, 12, 15), Xps = 20, seed=1)
X = tmpdata$X; Y = tmpdata$Y 

##### One Component #####
fit1 = msma(X, Y, comp=1, lambdaX=2, lambdaY=1:3)
summary(fit1)

}
\concept{print}
