## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 5)
library(msma)

## ----data---------------------------------------------------------------------
dat <- simdata(n = 35, rho = 0.8, Xps = c(4, 4), Yps = c(3, 3), seed = 4)
X <- dat$X
Y <- dat$Y

## ----ncomp-bic----------------------------------------------------------------
search_comp <- ncompsearch(X, comps = 1:3, criterion = "BIC", intseed = 1)
search_comp

## ----ncomp-plot---------------------------------------------------------------
plot(search_comp)

## ----nested-search, eval=FALSE------------------------------------------------
# search_nested <- ncompsearch(
#   X,
#   comps = list(1:4, 1:3),
#   criterion = "BIC",
#   intseed = 1
# )

## ----regularization-search, eval=FALSE----------------------------------------
# search_lambda <- regparasearch(
#   X = X,
#   comp = 2,
#   criterion = "BIC",
#   maxrep = 5,
#   intseed = 1
# )
# search_lambda

## ----combined-search, eval=FALSE----------------------------------------------
# opt <- optparasearch(
#   X = X,
#   search.method = "ncomp1st",
#   criterion = "BIC",
#   intseed = 1
# )
# 
# fit <- msma(
#   X = X,
#   comp = opt$optncomp,
#   lambdaX = opt$optlambdaX,
#   lambdaXsup = opt$optlambdaXsup,
#   intseed = 1
# )

## ----pls-selection, eval=FALSE------------------------------------------------
# opt_pls <- optparasearch(
#   X = X, Y = Y,
#   search.method = "regparaonly",
#   criterion = "BIC",
#   intseed = 1
# )
# 
# fit_pls <- msma(
#   X = X, Y = Y,
#   comp = opt_pls$optncomp,
#   lambdaX = opt_pls$optlambdaX,
#   lambdaY = opt_pls$optlambdaY,
#   lambdaXsup = opt_pls$optlambdaXsup,
#   lambdaYsup = opt_pls$optlambdaYsup,
#   intseed = 1
# )

## ----cv-example, eval=FALSE---------------------------------------------------
# cv <- cvmsma(
#   X = X, Y = Y,
#   comp = 1,
#   lambdaX = c(0.1, 0.1),
#   lambdaY = c(0.1, 0.1),
#   nfold = 5,
#   seed = 1,
#   intseed = 1
# )
# cv

## ----snmf-selection, eval=FALSE-----------------------------------------------
# search_snmf <- ncompsearch(
#   X,
#   comps = list(1:3, 1:3),
#   criterion = "BIC",
#   sprmethod = "sNMF",
#   nneg = "posneg",
#   intseed = 1
# )

## ----session-info-------------------------------------------------------------
sessionInfo()

