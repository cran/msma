## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#>")
library(msma)

## ----example------------------------------------------------------------------
set.seed(4)
X <- scale(matrix(rnorm(50*12),50,12))
Z <- as.numeric(scale(.7*X[,1]-.4*X[,5]+rnorm(50)))
fit <- msma(X,Z=Z,comp=3,lambdaX=.1,muX=.8,
            structure.method="soft",gammaX=.1,niterS4=30,
            scaling=FALSE,intseed=4)
fit$W
fit$overlap_all
fit$overlap_selected
fit$diagnostics

## ----session------------------------------------------------------------------
sessionInfo()

## ----selection, eval=FALSE----------------------------------------------------
# selection <- s4pca_conformal_select(
#   X, Z,
#   lambdaX = c(0.05, 0.10),
#   gammaX = c(0, 0.10),
#   comp = 2:3,
#   repeats = 5,
#   alpha = 0.10,
#   max_overlap = 0.50,
#   muX = 0.8,
#   intseed = 4
# )
# selection$selected
# selection$summary
# selection$fit

