% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/src.r
\docType{methods}
\name{regparasearch}
\alias{regparasearch}
\alias{print.regparasearch}
\title{Regularized Parameters Search}
\usage{
regparasearch(
  X,
  Y = NULL,
  Z = NULL,
  eta = 1,
  type = "lasso",
  inX = NULL,
  inY = NULL,
  inXsup = NULL,
  inYsup = NULL,
  muX = 0,
  muY = 0,
  comp = 1,
  nfold = 5,
  maxrep = 3,
  minpct = 0,
  maxpct = 1,
  criterion = c("CV", "BIC")[1],
  whichselect = NULL,
  homo = NULL,
  intseed = 1
)

\method{print}{regparasearch}(x, ...)
}
\arguments{
\item{X}{a matrix or list of matrices indicating the explanatory variable(s). This parameter is required.}

\item{Y}{a matrix or list of matrices indicating objective variable(s). This is optional. If there is no input for Y, then PCA is implemented.}

\item{Z}{a vector, response variable(s) for implementing the supervised version of (multiblock) PCA or PLS. This is optional. The length of Z is the number of subjects. If there is no input for Z, then unsupervised PLS/PCA is implemented.}

\item{eta}{numeric scalar indicating the parameter indexing the penalty family. This version contains only choice 1.}

\item{type}{a character, indicating the penalty family. In this version, only one choice is available: "lasso."}

\item{inX}{a (list of) numeric vector to specify the variables of X which are always in the model.}

\item{inY}{a (list of) numeric vector to specify the variables of X which are always in the model.}

\item{inXsup}{a (list of) numeric vector to specify the blocks of X which are always in the model.}

\item{inYsup}{a (list of) numeric vector to specify the blocks of Y which are always in the model.}

\item{muX}{a numeric scalar for the weight of X for the supervised case. 0 <= muX <= 1.}

\item{muY}{a numeric scalar for the weight of Y for the supervised case. 0 <= muY <= 1.}

\item{comp}{numeric scalar for the maximum number of componets to be considered.}

\item{nfold}{number of folds. Default is 5.}

\item{maxrep}{numeric scalar for the number of iteration.}

\item{minpct}{percent of minimum candidate parameters.}

\item{maxpct}{percent of maximum candidate parameters.}

\item{criterion}{a character, the evaluation criterion, "CV" for cross-validation, based on a matrix element-wise error, and "BIC" for Bayesian information criteria. The "BIC" is the default.}

\item{whichselect}{which blocks selected.}

\item{homo}{same parameters.}

\item{intseed}{seed number for the random number in the parameter estimation algorithm.}

\item{x}{an object of class "\code{regparasearch}", usually, a result of a call to \code{\link{regparasearch}}}

\item{...}{further arguments passed to or from other methods.}
}
\value{
\item{optlambdaX}{Optimal parameters for X}

\item{optlambdaY}{Optimal parameters for Y}

\item{mincriterion}{Minimum of criterion values}

\item{criterions}{Resulting criterion value}

\item{pararange}{Range of candidates parameters}
}
\description{
Regularized parameters search method for "msma".
}
\details{
This is a function for identifying the regularized parameters of sparseness lambdaX and lambdaY for \code{msma}. The initial range of candidates is computed based on fit, with regularized parameter values of 0. A binary search is conducted for dividing the parameter range into two regions. The representative value for the region is a median value, and the optimal region is selected using the minimum criteria obtained from the fit with that median value. The CV error or BIC can be used as criteria. The selected region is also divided into two region and the same process is iterated by maxrep times. Thus, the final median value in the selected region is set to be the optimal regularized parameter. The search is conducted with combinations of parameters for X and Y. The range of candidates for regularized parameters can be restricted, with a percentile of the limit (minimum or maximum) for the range.
}
\examples{
##### data #####
tmpdata = simdata(n = 50, rho = 0.8, Yps = c(10, 12, 15), Xps = c(20, 15), seed=1)
X = tmpdata$X; Y = tmpdata$Y 

##### Regularized parameters search #####
opt1 = regparasearch(X, Y, comp=1, criterion="BIC", maxrep=2, 
whichselect=c("X", "Y", "Xsup", "Ysup"))
opt1
fit4 = msma(X, Y, comp=1, lambdaX=opt1$optlambdaX, lambdaY=opt1$optlambdaY, 
lambdaXsup=opt1$optlambdaXsup, lambdaYsup=opt1$optlambdaYsup)
fit4
summary(fit4)

}
