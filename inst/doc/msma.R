## ----eval=FALSE---------------------------------------------------------------
#  if(!require("msma")) install.packages("msma")

## ----eval=TRUE----------------------------------------------------------------
library(msma)

## ----echo=FALSE, error=FALSE, eval=FALSE--------------------------------------
#  #source("../package/Ver3_nest/src.r")
#  #predir = unlist(strsplit(getwd(), "ADNI"))[1]
#  #source(paste(predir, "ADNI/Multimodal/program/package/Ver2/src.r", sep="")); library(mvtnorm)

## -----------------------------------------------------------------------------
dataset0 = simdata(n = 50, rho = 0.8, Yps = c(3, 4, 5), Xps = 3, seed=1)
X0 = dataset0$X; Y0 = dataset0$Y 

## -----------------------------------------------------------------------------
fit01 = msma(X0, Y0, comp=1, lambdaX=0.05, lambdaY=1:3)
fit01

## -----------------------------------------------------------------------------
plot(fit01)

## -----------------------------------------------------------------------------
fit02 = msma(X0, Y0, comp=2, lambdaX=0.03, lambdaY=0.01*(1:3))
fit02

## -----------------------------------------------------------------------------
dataset1 = simdata(n = 50, rho = 0.8, Yps = 5, Xps = 5, seed=1)
X1 = dataset1$X[[1]]; Y1 = dataset1$Y 

## -----------------------------------------------------------------------------
(fit111 = msma(X1, comp=5))

## -----------------------------------------------------------------------------
fit111$wbX

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit111, axes = 1, plottype="bar", cex.names=0.7, las=2)
plot(fit111, axes = 2, plottype="bar", cex.names=0.7, las=2)

## -----------------------------------------------------------------------------
lapply(fit111$sbX, head)

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit111, v="score", axes = 1:2, plottype="scatter")
plot(fit111, v="score", axes = 2:3, plottype="scatter")

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit111, v="cpev", ylim=c(0.7, 1))

## -----------------------------------------------------------------------------
(fit1112 = prcomp(X1, scale=TRUE))
summary(fit1112)

## -----------------------------------------------------------------------------
fit1113 = msma(X1, comp=5, ceps=0.0000001)
fit1113$wbX

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
biplot(fit1112)
plot(-fit1113$sbX[[1]][,1:2],xlab="Component 1",ylab="Component 2")

## -----------------------------------------------------------------------------
(fit112 = msma(X1, comp=5, lambdaX=0.1))
par(mfrow=c(1,2))
plot(fit112, axes = 1, plottype="bar", las=2)
plot(fit112, axes = 2, plottype="bar", las=2)

## -----------------------------------------------------------------------------
set.seed(1); Z = rbinom(50, 1, 0.5)

## -----------------------------------------------------------------------------
(fit113 = msma(X1, Z=Z, comp=5, lambdaX=0.02))

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit113, axes = 1, plottype="bar", las=2)
plot(fit113, axes = 2, plottype="bar", las=2)

## -----------------------------------------------------------------------------
(fit121 = msma(X1, Y1, comp=2))

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit121, axes = 1, XY="XY")
plot(fit121, axes = 2, XY="XY")

## -----------------------------------------------------------------------------
(fit122 = msma(X1, Y1, comp=2, lambdaX=0.5, lambdaY=0.5))
par(mfrow=c(1,2))
plot(fit122, axes = 1, XY="XY")
plot(fit122, axes = 2, XY="XY")

## -----------------------------------------------------------------------------
(fit123 = msma(X1, Y1, Z, comp=2, lambdaX=0.5, lambdaY=0.5))
par(mfrow=c(1,2))
plot(fit123, axes = 1, XY="XY")
plot(fit123, axes = 2, XY="XY")

## -----------------------------------------------------------------------------
dataset2 = simdata(n = 50, rho = 0.8, Yps = c(2, 3), Xps = c(3, 4), seed=1)
X2 = dataset2$X; Y2 = dataset2$Y 

## -----------------------------------------------------------------------------
class(X2)

## -----------------------------------------------------------------------------
length(X2)

## -----------------------------------------------------------------------------
lapply(X2, dim)

## -----------------------------------------------------------------------------
(fit211 = msma(X2, comp=1))

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit211, axes = 1, plottype="bar", block="block", las=2)
plot(fit211, axes = 1, plottype="bar", block="super")

## -----------------------------------------------------------------------------
(fit212 = msma(X2, comp=1, lambdaX=c(0.5, 0.5)))

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit212, axes = 1, plottype="bar", block="block", las=2)
plot(fit212, axes = 1, plottype="bar", block="super")

## -----------------------------------------------------------------------------
(fit213 = msma(X2, Z=Z, comp=1, lambdaX=c(0.5, 0.5)))

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit213, axes = 1, plottype="bar", block="block", las=2)
plot(fit213, axes = 1, plottype="bar", block="super")

## -----------------------------------------------------------------------------
(fit214 = msma(X2, comp=c(2,3)))

## -----------------------------------------------------------------------------
fit214$wbX

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit214, axes = 1, axes2 = 1, plottype="bar", block="block", las=2)
plot(fit214, axes = 2, axes2 = 1, plottype="bar", block="block", las=2)

## -----------------------------------------------------------------------------
fit214$wsX

## -----------------------------------------------------------------------------
par(mfrow=c(2,3))
for(j in 1:2) for(i in 1:3) plot(fit214, axes = j, axes2 = i, plottype="bar", block="super")

## -----------------------------------------------------------------------------
(fit221 = msma(X2, Y2, comp=1))

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit221, axes = 1, plottype="bar", block="block", XY="X", las=2)
plot(fit221, axes = 1, plottype="bar", block="super", XY="X")

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit221, axes = 1, plottype="bar", block="block", XY="Y", las=2)
plot(fit221, axes = 1, plottype="bar", block="super", XY="Y")

## -----------------------------------------------------------------------------
(fit222 = msma(X2, Y2, comp=1, lambdaX=c(0.5, 0.5), lambdaY=c(0.5, 0.5)))

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit222, axes = 1, plottype="bar", block="block", XY="X", las=2)
plot(fit222, axes = 1, plottype="bar", block="super", XY="X")

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit222, axes = 1, plottype="bar", block="block", XY="Y", las=2)
plot(fit222, axes = 1, plottype="bar", block="super", XY="Y")

## -----------------------------------------------------------------------------
(fit223 = msma(X2, Y2, Z, comp=1, lambdaX=c(0.5, 0.5), lambdaY=c(0.5, 0.5)))

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit223, axes = 1, plottype="bar", block="block", XY="X", las=2)
plot(fit223, axes = 1, plottype="bar", block="super", XY="X")

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(fit223, axes = 1, plottype="bar", block="block", XY="Y", las=2)
plot(fit223, axes = 1, plottype="bar", block="super", XY="Y")

## -----------------------------------------------------------------------------
(ncomp11 = ncompsearch(X1, comps = c(1, 5, 10*(1:2)), nfold=5))
plot(ncomp11)

## -----------------------------------------------------------------------------
(ncomp12 = ncompsearch(X1, comps = 20, criterion="BIC"))
plot(ncomp12)

## -----------------------------------------------------------------------------
ncomp21 = ncompsearch(X2, Y2, comps = c(1, 5, 10*(1:2)), nfold=5)
plot(ncomp21)

## -----------------------------------------------------------------------------
dataset3 = simdata(n = 50, rho = 0.8, Yps = rep(4, 5), Xps = rep(4, 5), seed=1)
X3 = dataset3$X; Y3 = dataset3$Y 

## -----------------------------------------------------------------------------
ncomp31 = ncompsearch(X3, comps = 20, criterion="BIC")
plot(ncomp31)

## -----------------------------------------------------------------------------
(ncomp32 = ncompsearch(X3, comps = list(20, 20), criterion="BIC"))

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(ncomp32,1)
plot(ncomp32,2)

## -----------------------------------------------------------------------------
ncomp41 = ncompsearch(X3, Y3, comps = c(1, 5, 10*(1:2)), criterion="BIC")
plot(ncomp41)

## -----------------------------------------------------------------------------
criteria = c("BIC", "CV")
search.methods = c("regparaonly", "regpara1st", "ncomp1st", "simultaneous")

## -----------------------------------------------------------------------------
(opt11 = optparasearch(X1, search.method = "regparaonly", criterion="BIC", comp=ncomp11$optncomp))
(fit311 = msma(X1, comp=opt11$optncomp, lambdaX=opt11$optlambdaX))

## -----------------------------------------------------------------------------
(opt12 = optparasearch(X1, search.method = "ncomp1st", criterion="BIC"))
(fit312 = msma(X1, comp=opt12$optncomp, lambdaX=opt12$optlambdaX))

## -----------------------------------------------------------------------------
(opt13 = optparasearch(X1, search.method = "regpara1st", criterion="BIC"))
(fit313 = msma(X1, comp=opt13$optncomp, lambdaX=opt13$optlambdaX))

## -----------------------------------------------------------------------------
(opt14 = optparasearch(X1, search.method = "simultaneous", criterion="BIC"))
(fit314 = msma(X1, comp=opt14$optncomp, lambdaX=opt14$optlambdaX))

## -----------------------------------------------------------------------------
(opt132 = optparasearch(X1, search.method = "ncomp1st", criterion="BIC", maxpct4ncomp=0.5))
(fit3132 = msma(X1, comp=opt132$optncomp, lambdaX=opt132$optlambdaX))

## -----------------------------------------------------------------------------
(opt133 = optparasearch(X1, search.method = "regpara1st", criterion="BIC", comp=5))
(fit3133 = msma(X1, comp=opt133$optncomp, lambdaX=opt133$optlambdaX))

## -----------------------------------------------------------------------------
(opt21 = optparasearch(X2, Y2, search.method = "regparaonly", criterion="BIC"))
(fit321 = msma(X2, Y2, comp=opt21$optncomp, lambdaX=opt21$optlambdaX, lambdaY=opt21$optlambdaY))

## -----------------------------------------------------------------------------
(opt31 = optparasearch(X3, search.method = "regparaonly", criterion="BIC"))
(fit331 = msma(X3, comp=opt31$optncomp, lambdaX=opt31$optlambdaX, lambdaXsup=opt31$optlambdaXsup))

## -----------------------------------------------------------------------------
(opt32 = optparasearch(X3, search.method = "regparaonly", criterion="BIC", whichselect="X"))
(fit332 = msma(X3, comp=opt32$optncomp, lambdaX=opt32$optlambdaX, lambdaXsup=opt32$optlambdaXsup))

## -----------------------------------------------------------------------------
(opt33 = optparasearch(X3, search.method = "regparaonly", criterion="BIC", whichselect="Xsup"))
(fit333 = msma(X3, comp=opt33$optncomp, lambdaX=opt33$optlambdaX, lambdaXsup=opt33$optlambdaXsup))

## -----------------------------------------------------------------------------
(opt341 = optparasearch(X3, search.method = "ncomp1st", criterion="BIC", comp=c(8, 8)))
(fit341 = msma(X3, comp=opt341$optncomp, lambdaX=opt341$optlambdaX, lambdaXsup=opt341$optlambdaXsup))

## -----------------------------------------------------------------------------
(opt342 = optparasearch(X3, search.method = "regparaonly", criterion="BIC", comp=c(4, 5)))
(fit342 = msma(X3, comp=opt342$optncomp, lambdaX=opt342$optlambdaX, lambdaXsup=opt342$optlambdaXsup))

## -----------------------------------------------------------------------------
(opt344 = optparasearch(X3, search.method = "regpara1st", criterion="BIC", comp=c(8, 8)))
(fit344 = msma(X3, comp=opt344$optncomp, lambdaX=opt344$optlambdaX, lambdaXsup=opt344$optlambdaXsup))

## ----eval=FALSE---------------------------------------------------------------
#  (opt345 = optparasearch(X3, search.method = "simultaneous", criterion="BIC", comp=c(8, 8)))
#  (fit345 = msma(X3, comp=opt345$optncomp, lambdaX=opt345$optlambdaX))

## ----eval=FALSE---------------------------------------------------------------
#  (opt41 = optparasearch(X3, Y3, search.method = "regparaonly", criterion="BIC"))
#  (fit341 = msma(X3, Y3, comp=opt41$optncomp, lambdaX=opt41$optlambdaX, lambdaY=opt41$optlambdaY, lambdaXsup=opt41$optlambdaXsup, lambdaYsup=opt41$optlambdaYsup))

## -----------------------------------------------------------------------------
(opt42 = optparasearch(X3, Y3, search.method = "regparaonly", criterion="BIC", whichselect=c("Xsup","Ysup")))
(fit342 = msma(X3, Y3, comp=opt42$optncomp, lambdaX=opt42$optlambdaX, lambdaY=opt42$optlambdaY, lambdaXsup=opt42$optlambdaXsup, lambdaYsup=opt42$optlambdaYsup))

## -----------------------------------------------------------------------------
dataset4 = simdata(n = 50, rho = 0.8, Yps = rep(4, 2), Xps = rep(4, 3), seed=1)
X4 = dataset4$X; Y4 = dataset4$Y 

## -----------------------------------------------------------------------------
(opt43 = optparasearch(X4, Y4, search.method = "regparaonly", criterion="BIC"))
(fit343 = msma(X4, Y4, comp=opt43$optncomp, lambdaX=opt43$optlambdaX, lambdaY=opt43$optlambdaY, lambdaXsup=opt43$optlambdaXsup, lambdaYsup=opt43$optlambdaYsup))

