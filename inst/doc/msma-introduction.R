## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 5)
library(msma)

## ----data---------------------------------------------------------------------
dat <- simdata(n = 40, rho = 0.8, Xps = 5, Yps = 4, seed = 1)
X <- dat$X[[1]]
Y <- dat$Y[[1]]
set.seed(1)
Z <- rbinom(nrow(X), 1, 0.5)

dim(X)
dim(Y)

## ----pca----------------------------------------------------------------------
fit_pca <- msma(X, comp = 2)
fit_pca
summary(fit_pca)

## ----pca-results--------------------------------------------------------------
fit_pca$wbX
head(fit_pca$sbX[[1]])
fit_pca$cpevX

## ----pca-plot, fig.show='hold'------------------------------------------------
plot(fit_pca, axes = 1, plottype = "bar", las = 2)
plot(fit_pca, v = "score", axes = 1:2, plottype = "scatter")

## ----sparse-pca---------------------------------------------------------------
fit_spca <- msma(X, comp = 2, lambdaX = 0.10)
fit_spca$nzwbX
fit_spca$selectXnames

## ----supervised-pca-----------------------------------------------------------
fit_sup_pca <- msma(
  X = X, Z = Z, comp = 2,
  lambdaX = 0.05, muX = 0.20,
  intseed = 1
)
fit_sup_pca$predictiv

## ----pls----------------------------------------------------------------------
fit_pls <- msma(X = X, Y = Y, comp = 2)
fit_pls

## ----pls-plot, fig.show='hold'------------------------------------------------
plot(fit_pls, axes = 1, XY = "XY")
plot(fit_pls, axes = 2, XY = "XY")

## ----sparse-supervised-pls----------------------------------------------------
fit_spls <- msma(
  X = X, Y = Y, Z = Z, comp = 2,
  lambdaX = 0.10, lambdaY = 0.10,
  muX = 0.10, muY = 0.10,
  intseed = 1
)
fit_spls$nzwbX
fit_spls$nzwbY

## ----prediction---------------------------------------------------------------
pred <- predict(fit_pls, newX = X, newY = Y)
names(pred)

## ----session-info-------------------------------------------------------------
sessionInfo()

