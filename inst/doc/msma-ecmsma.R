## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 5)
library(msma)

## ----data---------------------------------------------------------------------
dat <- simdata(n = 45, rho = 0.8, Xps = c(5, 5, 5), Yps = 3, seed = 3)
X <- dat$X
names(X) <- paste0("block", seq_along(X))

## ----compare-methods----------------------------------------------------------
fit_pca <- msma(X, comp = c(2, 2), sprmethod = "PCA", intseed = 1)
fit_nmf <- msma(X, comp = c(2, 2), sprmethod = "NMF", intseed = 1)
fit_snmf <- msma(
  X, comp = c(2, 2),
  sprmethod = "sNMF",
  lambdaXsup = 0.05,
  intseed = 1
)

## ----nneg-options, eval=FALSE-------------------------------------------------
# fit_posneg <- msma(X, comp = c(2, 2), sprmethod = "NMF", nneg = "posneg")
# fit_absolute <- msma(X, comp = c(2, 2), sprmethod = "NMF", nneg = "absolute")
# fit_min <- msma(X, comp = c(2, 2), sprmethod = "NMF", nneg = "min")

## ----nonnegative--------------------------------------------------------------
all(unlist(fit_snmf$ssX) >= 0)
all(unlist(fit_snmf$wsX) >= 0)

## ----reproducibility----------------------------------------------------------
fit_snmf_2 <- msma(
  X, comp = c(2, 2),
  sprmethod = "sNMF",
  lambdaXsup = 0.05,
  intseed = 1
)
all.equal(fit_snmf$ssX, fit_snmf_2$ssX)
all.equal(fit_snmf$wsX, fit_snmf_2$wsX)

## ----clustering---------------------------------------------------------------
cluster_matrix <- vapply(
  fit_snmf$ssX,
  function(score) max.col(score, ties.method = "first"),
  integer(nrow(fit_snmf$ssX[[1]]))
)
colnames(cluster_matrix) <- names(fit_snmf$ssX)
head(cluster_matrix)
apply(cluster_matrix, 2, table)

## ----multivariate-z-----------------------------------------------------------
set.seed(3)
z1 <- rnorm(nrow(X[[1]]))
z2 <- 0.5 * z1 + rnorm(nrow(X[[1]]), sd = 0.5)
Z <- cbind(clinical = z1, biomarker = z2)

fit_multi_z <- msma(
  X = X, Z = Z,
  con4spv = c(0.7, 0.3),
  comp = 2,
  muX = 0.20,
  intseed = 1
)
fit_multi_z$predictiv

## ----one-column-z-------------------------------------------------------------
fit_z_vector <- msma(X, Z = z1, comp = 1, muX = 0.20, intseed = 1)
fit_z_matrix <- msma(X, Z = cbind(z1), comp = 1, muX = 0.20, intseed = 1)
fit_z_vector$call <- NULL
fit_z_matrix$call <- NULL
isTRUE(all.equal(fit_z_vector, fit_z_matrix, tolerance = 1e-8))

## ----session-info-------------------------------------------------------------
sessionInfo()

