## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7, fig.height = 5)
library(msma)

## ----data---------------------------------------------------------------------
dat <- simdata(
  n = 40, rho = 0.8,
  Xps = c(4, 5), Yps = c(3, 4),
  seed = 2
)
X <- dat$X
Y <- dat$Y
names(X) <- c("X_block_1", "X_block_2")
names(Y) <- c("Y_block_1", "Y_block_2")
lapply(X, dim)
lapply(Y, dim)

## ----mb-pca-------------------------------------------------------------------
fit_mb_pca <- msma(X = X, comp = 2, intseed = 1)
fit_mb_pca

## ----mb-components------------------------------------------------------------
lapply(fit_mb_pca$wbX, dim)
lapply(fit_mb_pca$sbX, dim)
lapply(fit_mb_pca$wsX, dim)
lapply(fit_mb_pca$ssX, dim)

## ----mb-pca-plot, fig.show='hold'---------------------------------------------
plot(fit_mb_pca, axes = 1, plottype = "bar", block = "block", las = 2)
plot(fit_mb_pca, axes = 1, plottype = "bar", block = "super")

## ----sparse-mb-pca------------------------------------------------------------
fit_sparse <- msma(
  X = X, comp = 2,
  lambdaX = c(0.10, 0.15),
  lambdaXsup = 0.05,
  intseed = 1
)
fit_sparse$nzwbX
fit_sparse$nzwsX

## ----nested-------------------------------------------------------------------
fit_nested <- msma(X = X, comp = c(2, 3), intseed = 1)
lapply(fit_nested$wsX, dim)
lapply(fit_nested$ssX, dim)

## ----nested-plot, fig.show='hold'---------------------------------------------
plot(fit_nested, axes = 1, axes2 = 1, plottype = "bar", block = "super")
plot(fit_nested, axes = 1, axes2 = 2, plottype = "bar", block = "super")

## ----supervised-mb-pca--------------------------------------------------------
set.seed(2)
Z <- rnorm(nrow(X[[1]]))
fit_supervised <- msma(
  X = X, Z = Z, comp = 2,
  lambdaX = c(0.10, 0.10),
  muX = 0.20,
  intseed = 1
)
fit_supervised$predictiv

## ----mb-pls-------------------------------------------------------------------
fit_mb_pls <- msma(
  X = X, Y = Y, comp = 2,
  lambdaX = c(0.10, 0.10),
  lambdaY = c(0.10, 0.10),
  intseed = 1
)
fit_mb_pls

## ----nested-mb-pls------------------------------------------------------------
fit_nested_pls <- msma(
  X = X, Y = Y, comp = c(2, 2),
  lambdaX = c(0.10, 0.10),
  lambdaY = c(0.10, 0.10),
  lambdaXsup = 0.05,
  lambdaYsup = 0.05,
  intseed = 1
)
lapply(fit_nested_pls$ssX, dim)
lapply(fit_nested_pls$ssY, dim)

## ----mb-prediction------------------------------------------------------------
Xnew <- lapply(X, function(block) block[1:5, , drop = FALSE])

prediction_mb <- predict(
  fit_mb_pca,
  newX = Xnew,
  preprocessing = "fitted",
  reconstruction = "fitted",
  output_scale = "original"
)

names(prediction_mb$X)
lapply(prediction_mb$X, colnames)
lapply(prediction_mb$sbX, dim)

## ----session-info-------------------------------------------------------------
sessionInfo()

