#' Multiblock Sparse Multivariable Analysis Package
#'
#' A Package for implementation of multiblock multivariable data analysis.
#'
#' @name msma-package
#' @aliases msma-package
#' @rdname msma-package
#' @docType package
#' @keywords documentation
#'
#' @author Atsushi Kawaguchi. \email{kawa_a24@@yahoo.co.jp}
#' @seealso \code{\link{msma}}
#' @references 
#' Kawaguchi A, Yamashita F (2017). Supervised Multiblock Sparse Multivariable Analysis with Application to Multimodal Brain Imaging Genetics. Biostatistics, 18(4) 651-665. 
#' @import mvtnorm
#' @importFrom graphics abline matplot plot 
#' @importFrom stats cor predict quantile rbinom rnorm runif
NULL

##########################################################################################################################################
#' Multiblock Sparse Multivariable Analysis
#'
#' This is a function for a matrix decomposition method incorporating sparse and supervised modeling for a multiblock multivariable data analysis 
#' 
#' \code{msma} requires at least one input X (a matrix or list). In this case, (multiblock) PCA is conducted. If Y is also specified, then a PLS is conducted using X as explanatory variables and Y as objective variables. This function scales each data matrix to a mean of 0 and variance of 1 in the default. The block structure can be represented as a list. If Z is also specified, a supervised version is implemented, and the degree is controlled by muX or muY, where 0 <= muX <= 1, 0 <= muY <= 1, and 0 <= muX + muY < 1. If a positive lambdaX or lambdaY is specified, then a sparse estimation based on the L1 penalty is implemented.
#'
#' @name msma
#' @aliases msma
#' @rdname msma
#' @docType methods
#' @export
#'
#' @param X a matrix or list of matrices indicating the explanatory variable(s). This parameter is required.
#' @param Y a matrix or list of matrices indicating objective variable(s). This is optional. If there is no input for Y, then PCA is implemented.
#' @param Z a vector, response variable(s) for implementing the supervised version of (multiblock) PCA or PLS. This is optional. The length of Z is the number of subjects. If there is no input for Z, then unsupervised PLS/PCA is implemented.
#' @param comp numeric scalar for the number of components to be considered.
#' @param lambdaX numeric vector of regularized parameters for X, with a length equal to the number of blocks. If lambdaX is omitted, no regularization is conducted.
#' @param lambdaY numeric vector of regularized parameters for Y, with a length equal to the number of blocks. If lambdaY is omitted, no regularization is conducted.
#' @param eta numeric scalar indicating the parameter indexing the penalty family. This version contains only choice 1.
#' @param type a character, indicating the penalty family. In this version, only one choice is available: "lasso."
#' @param inX a vector or list of numeric vectors specifying the variables in X, always included in the model
#' @param inY a vector or list of numeric vectors specifying the variables in Y, always included in the model
#' @param muX a numeric scalar for the weight of X for the supervised case. 0 <= muX <= 1.
#' @param muY a numeric scalar for the weight of Y for the supervised case. 0 <= muY <= 1. 
#' @param defmethod a character representing the deflation method. This version has only the choice "canonical." 
#' @param scaling a logical, indicating whether or not data scaling is performed. The default is TRUE.
#' @param verbose information
#' @param intseed seed number for the random number in the parameter estimation algorithm.
#' @param x an object of class "\code{msma}." Usually, a result of a call to \code{\link{msma}}
#' @param ... further arguments passed to or from other methods.
#' @return \item{dmode}{Indicates mode "PLS" or "PCA"}
#' @return \item{X}{Scaled X, which has a list form.}
#' @return \item{Y}{Scaled Y, which has a list form.}
#' @return \item{Xscale}{Scaling information for X. The mean and standard deviation values for each block of X are returned.}
#' @return \item{Yscale}{Scaling information for Y. The mean and standard deviation values for each block of Y are returned.}
#' @return \item{comp}{Number of components}
#' @return \item{wbX}{Block loading for X. The list has the same length as that of the input list X (number of blocks) and consists of a matrix. The number of variables is present in the row and the number of components is present in the column.}
#' @return \item{sbX}{Block score for X. The list has the same length as that of the input list X (number of blocks) and consists of a matrix, with the number of subjects in the row and the number of components in the column.}
#' @return \item{wbY}{Block loading for Y. The list has same length as that of the input list Y (number of blocks) and consists of a matrix, with the number of variables in the row and the number of components in the column.}
#' @return \item{sbY}{Block score for Y. The list has same length as that of the input list Y (number of blocks) and consists of a matrix, with the number of subjects in the row and the number of components in the column.}
#' @return \item{ssX}{Super score for X. In the matrix, the number of subjects is in the row and the number of components is in the column.}
#' @return \item{wsX}{Super loading for X. In the matrix, the number of blocks is in the row and the number of components is in the column.}
#' @return \item{ssY}{Super score for Y. In the matrix, the number of subjects is in the row and the number of components is in the column.}
#' @return \item{wsY}{Super loading for Y. In the matrix, the number of blocks is in the row and the number of components is in the column.}
#' @return \item{nzwbX}{Number of nonzeros in block loading for X}
#' @return \item{nzwbY}{Number of nonzeros in block loading for Y}
#' @return \item{selectXnames}{Names of selected variables for X. This returns from the original names of X}
#' @return \item{selectYnames}{Names of selected variables for Y. This returns from the original names of Y}
#'
#' @examples
#' ##### data #####
#' tmpdata = simdata(n = 50, rho = 0.8, Yps = c(10, 12, 15), Xps = 20, seed=1)
#' X = tmpdata$X; Y = tmpdata$Y 
#' 
#' ##### One Component #####
#' fit1 = msma(X, Y, comp=1, lambdaX=2, lambdaY=1:3)
#' fit1
#' 
#' ##### Two Component #####
#' fit2 = msma(X, Y, comp=2, lambdaX=2, lambdaY=1:3)
#' fit2
#' 
#' ##### Matrix data #####
#' sigma = matrix(0.8, 10, 10)
#' diag(sigma) = 1
#' X2 = rmvnorm(50, rep(0, 10), sigma)
#' Y2 = rmvnorm(50, rep(0, 10), sigma)
#' 
#' fit3 = msma(X2, Y2, comp=1, lambdaX=2, lambdaY=2)
#' fit3
#' 
#' ##### Sparse Principal Component Analysis #####
#' fit5 = msma(X2, comp=5, lambdaX=2.5)
#' summary(fit5)
#' 
msma = function(X, ...) UseMethod("msma")

#' @rdname msma
#' @method msma default
#' @export
msma.default = function(X, Y=NULL, Z=NULL, comp=2, lambdaX=NULL, lambdaY=NULL, eta=1, type="lasso", inX=NULL, inY=NULL, muX = 0, muY = 0, defmethod = "canonical", scaling = TRUE, verbose=FALSE, intseed=1, ...)
{
##### Requirement #####
if(missing(X)) stop("X should be specified")
#if(missing(Y)) stop("Y should be specified")
if(is.null(Y)){ 
dmode = "PCA"
Y = X
if(!is.null(lambdaY)) stop("lambdaY should be specified")
}else
{
dmode = "PLS"
}

if(dmode == "PCA"){ lambdaY=lambdaX}

if(class(X) != "list") X = list(X)
if(class(Y) != "list") Y = list(Y)

if(class(X) != "list") stop("X should be a list")
if(class(Y) != "list") stop("Y should be a list")

Ynb = length(Y)
Xnb = length(X)

if(!is.null(lambdaX)){if(length(lambdaX) != Xnb) stop("lambdaX should be same length as the number of blocks of X")}
if(!is.null(lambdaY)){if(length(lambdaY) != Ynb) stop("lambdaY should be same length as the number of blocks of Y")}

if(!is.null(inX)){if(class(inX) != "list") inX = list(inX)}
if(!is.null(inY)){if(class(inY) != "list") inY = list(inY)}

Yps = unlist(lapply(Y, ncol))
Xps = unlist(lapply(X, ncol))

Yns = unlist(lapply(Y, nrow))
Xns = unlist(lapply(X, nrow))

if(length(unique(Yns)) == 1){Yns = unique(Yns)}else{stop("numbers of rows should be same among Y")}
if(length(unique(Xns)) == 1){Xns = unique(Xns)}else{stop("numbers of rows should be same among X")}
if(Yns == Xns){n = Xns}else{stop("numbers of rows should be same between X and Y")}
if(!is.null(Z)){ if(length(Z) != n) stop("length of Z should be same as row of X or Y")}

if(!is.null(inX)){
extest = unlist(lapply(1:Xnb, function(a) inX[[a]] %in% (1:Xps[a]) ))
if(!(all(extest))) stop("inX should be specified the column number of X")
}
if(!is.null(inY)){
extest = unlist(lapply(1:Ynb, function(a) inY[[a]] %in% (1:Yps[a]) ))
if(!(all(extest))) stop("inY should be specified the column number of Y")
}

if(muX < 0) stop("muX should be positive")
if(muY < 0) stop("muY should be positive")

if(eta != 1) stop("Sorry! in this version, eta should be only 1")
if(type != "lasso") stop("Sorry! in this version, type only ''lasso'' ")

########## Start ##########
X = lapply(X, function(x) scale(x, scale=scaling))
Y = lapply(Y, function(x) scale(x, scale=scaling))
if(!is.null(Z)) Z = scale(Z, scale=scaling)

Xscale = lapply(X, attributes)
Yscale = lapply(Y, attributes)

if(scaling){
for(x in 1:length(Xscale)) Xscale[[x]]$`scaled:scale` = ifelse(Xscale[[x]]$`scaled:scale`==0, 1, Xscale[[x]]$`scaled:scale`)
for(x in 1:length(Yscale)) Yscale[[x]]$`scaled:scale` = ifelse(Yscale[[x]]$`scaled:scale`==0, 1, Yscale[[x]]$`scaled:scale`)
}

X = lapply(X, function(x) apply(x, 2, function(x2){x2[is.na(x2)] = 0; x2}))
Y = lapply(Y, function(y) apply(y, 2, function(y2){y2[is.na(y2)] = 0; y2}))

tvarX = lapply(X, function(x) sum(x^2))
tvarY = lapply(Y, function(x) sum(x^2))

########## Iteration for components Start ##########
Xd = X; Yd = Y; out = cpevX = cpevY = list(); bic = bicX = bicY = reproduct = predictiv = rep(0, comp)
for(ncomp in 1:comp){

##### Fit #####
out[[ncomp]] = out1 = msma_OneComp(X=Xd, Yd, Z, lambdaX, lambdaY, eta, type, inX, inY, muX, muY, dmode, verbose, intseed)

##### arrange #####
sbX = lapply(1:Xnb, function(x){ do.call(cbind, lapply(out, function(y) y$sbX[[x]]))})
sbY = lapply(1:Ynb, function(x){ do.call(cbind, lapply(out, function(y) y$sbY[[x]]))})

wbX = lapply(1:Xnb, function(x){ do.call(cbind, lapply(out, function(y) y$wbX[[x]]))})
wbY = lapply(1:Ynb, function(x){ do.call(cbind, lapply(out, function(y) y$wbY[[x]]))})

##### Number of nonzero #####
nzwbX = do.call(rbind, lapply(wbX, function(z) apply(z, 2, function(z2) sum(z2!=0))))
nzwbY = do.call(rbind, lapply(wbY, function(z) apply(z, 2, function(z2) sum(z2!=0))))

##### Prediction #####
predX = lapply(1:Xnb, function(x) sbX[[x]] %*% project(X[[x]], sbX[[x]]))
predY = lapply(1:Ynb, function(x) sbY[[x]] %*% project(Y[[x]], sbY[[x]]))

cpevX[[ncomp]] = sapply(1:Xnb, function(x) sum(predX[[x]]^2)/tvarX[[x]] )
cpevY[[ncomp]] = sapply(1:Ynb, function(x) sum(predY[[x]]^2)/tvarY[[x]] )

##### Deflation #####
if(defmethod == "canonical"){
Xd = lapply(1:Xnb, function(x) X[[x]] - predX[[x]])
Yd = lapply(1:Ynb, function(x) Y[[x]] - predY[[x]])
}else
if(defmethod == "pls"){
Xd = lapply(1:Xnb, function(x) X[[x]] - predX[[x]])
Yd = lapply(1:Ynb, function(x) Y[[x]] - predY[[x]])
}

##### Predictivity and Reproductivity #####
reproduct[ncomp] = cor(unlist(lapply(Y, c)), unlist(lapply(predY, c)))^2
predictiv[ncomp] = cor(out1$ssY, Z)^2

##### Information Criteria #####
if(dmode == "PCA"){ 
tmpse = sum(unlist(lapply(Xd, function(x) sum(x^2))))
#print(tmpse)
bic[ncomp] = log(tmpse/(n*sum(Xps))) + log(n*sum(Xps))/(n*sum(Xps))*(sum(nzwbX))
}else
{
#tmpse = sum(unlist(lapply(Xd, function(x) sum(x^2)))) + sum(unlist(lapply(Yd, function(x) sum(x^2))))
tmpse = sum(unlist(lapply(Yd, function(x) sum(x^2))))
bicY[ncomp] = log(tmpse/(n*sum(Yps))) + log(n*sum(Yps))/(n*sum(Yps))*(sum(nzwbY))
#bic[ncomp] = log(tmpse/(n*sum(c(Xps,Yps)))) + log(n*sum(c(Xps,Yps)))/(n*sum(c(Xps,Yps)))*(sum(c(nzwbX,nzwbY)))
tmpse = sum(unlist(lapply(Xd, function(x) sum(x^2))))
bicX[ncomp] = log(tmpse/(n*sum(Xps))) + log(n*sum(Xps))/(n*sum(Xps))*(sum(nzwbX))
}

}

##### arrange #####
ssX = do.call(cbind, lapply(out, function(y) y$ssX))
wsX = do.call(cbind, lapply(out, function(y) y$wsX))

ssY = do.call(cbind, lapply(out, function(y) y$ssY))
wsY = do.call(cbind, lapply(out, function(y) y$wsY))

cpevX = do.call(cbind, cpevX)
cpevY = do.call(cbind, cpevY)

if(ncomp==1){
avX = cpevX; avY = cpevY
}else
if(ncomp==2){
avX = cbind(cpevX[,1], cpevX[,2]-cpevX[,1])
avY = cbind(cpevY[,1], cpevY[,2]-cpevY[,1])
}else{
avX = cbind(cpevX[,1], t(apply(cpevX, 1, diff)))
avY = cbind(cpevY[,1], t(apply(cpevY, 1, diff)))
}
colnames(cpevY) = colnames(cpevX) = colnames(avY) = colnames(avX) = paste("comp", 1:comp, sep="")
rownames(cpevX) = rownames(avX) = paste("block", 1:nrow(nzwbX), sep="")
rownames(cpevY) = rownames(avY) = paste("block", 1:nrow(nzwbY), sep="")

colnames(nzwbX) = colnames(nzwbY) = paste("comp", 1:comp, sep="")
rownames(nzwbX) = paste("block", 1:nrow(nzwbX), sep="")
rownames(nzwbY) = paste("block", 1:nrow(nzwbY), sep="")

##### Names of variables #####
selectXnames = lapply(1:length(wbX), function(z) apply(wbX[[z]], 2, function(z2){ 
if(is.null(colnames(X[[z]]))) colnames(X[[z]]) = paste("X", z, 1:ncol(X[[z]]), sep=".")
colnames(X[[z]])[z2!=0]
}))
selectYnames = lapply(1:length(wbY), function(z) apply(wbY[[z]], 2, function(z2){ 
if(is.null(colnames(Y[[z]]))) colnames(Y[[z]]) = paste("Y", z, 1:ncol(Y[[z]]), sep=".")
colnames(Y[[z]])[z2!=0]
}))

##### Output #####
if(dmode == "PCA"){ 
out = list(call = match.call(), 
dmode=dmode,
X=X, Z=Z, Xscale=Xscale, comp = comp,
muX=muX, 
wbX=wbX, sbX=sbX, ssX=ssX, wsX=wsX,
nzwbX = nzwbX,
selectXnames = selectXnames,
cpevX = cpevX,
avX = avX,
bic = bic,
reproduct = reproduct, predictiv=predictiv
)
}else
{
out = list(call = match.call(), 
dmode=dmode,
X=X, Y=Y, Z=Z,
Xscale=Xscale, Yscale=Yscale, comp = comp,
muX=muX, muY=muY,
wbX=wbX, sbX=sbX, wbY=wbY, sbY=sbY, ssX=ssX, wsX=wsX, ssY=ssY, wsY=wsY, 
nzwbX = nzwbX, nzwbY = nzwbY,
selectXnames = selectXnames, selectYnames = selectYnames,
cpevX = cpevX, cpevY = cpevY,
avX = avX, avY = avY,
bic = bicY, bicX = bicX, bicY = bicY,
reproduct = reproduct, predictiv=predictiv
)
}
class(out) = "msma"
out
}

#' @rdname msma
#' @method print msma
#' @export
print.msma = function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\n")
cat("Numbers of non-zeros for X: \n")
print(x$nzwbX)
cat("\n")
if(x$dmode == "PLS"){
cat("Numbers of non-zeros for Y: \n")
print(x$nzwbY)
cat("\n")
}
}

#' Summarizing Fits
#'
#' summary method for class "msma". 
#'
#' This function provides a summary of results.
#'
#' @name summary.msma
#' @aliases summary.msma
#' @rdname summary.msma
#' @method summary msma
#' @docType methods
#' @export
#'
#' @param object,x an object of class "\code{msma}." Usually, a result of a call to \code{\link{msma}}
#' @param ... further arguments passed to or from other methods.
#'
#' @examples
#'##### data #####
#'tmpdata = simdata(n = 50, rho = 0.8, Yps = c(10, 12, 15), Xps = 20, seed=1)
#'X = tmpdata$X; Y = tmpdata$Y 
#'
#'##### One Component #####
#'fit1 = msma(X, Y, comp=1, lambdaX=2, lambdaY=1:3)
#'summary(fit1)
#'
summary.msma = function(object, ...)
{
if(object$dmode == "PLS"){
test = predict(object, newX=object$X, newY=object$Y)
err = (matserr(object$X, test$X) + matserr(object$Y, test$Y))/2
}else{
test = predict(object, newX=object$X, newY=NULL)
err = matserr(object$X, test$X)
}

res = list(call=object$call, err=err)
class(res) = "summary.msma"
res
}

#' @rdname msma
#' @method plot msma
#' @export
plot.msma = function(x,...){
matplot(t(x$cpevX), type="b", xlab="Components", ylab="", ylim=c(0,1), ...)
if(x$dmode == "PLS"){
matplot(t(x$cpevY), type="b", lty=2, add = TRUE)
}
}

#' @rdname summary.msma
#' @method print summary.msma
#' @family print
#' @export
print.summary.msma = function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\n")
cat("Error : ")
cat(round(x$err, 3))
cat("\n")
}

#' Prediction
#'
#' predict method for class "msma". 
#'
#' This function produces a prediction from new data based on \code{\link{msma}} fit. It is mainly used in cross-validation
#'
#' @name predict.msma
#' @aliases predict.msma
#' @rdname predict.msma
#' @method predict msma
#' @docType methods
#' @export
#'
#' @param object an object of class "\code{msma}." Usually, a result of a call to \code{msma}
#' @param newX a matrix in which to look for the variables used to predict X. This is required.
#' @param newY a matrix in which to look for the variables used to predict Y.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \item{X}{predicted X}
#' @return \item{sbX}{block score for X}
#' @return \item{Y}{predicted Y}
#' @return \item{sbY}{block score for Y}
#'
#' @examples
#' ##### data #####
#' tmpdata = simdata(n = 50, rho = 0.8, Yps = c(10, 12, 15), Xps = 20, seed=1)
#' X = tmpdata$X; Y = tmpdata$Y 
#' 
#' ##### Two Component #####
#' fit2 = msma(X, Y, comp=2, lambdaX=2, lambdaY=1:3)
#' summary(fit2)
#' 
#' ##### Predict #####
#' test = predict(fit2, newX=X, newY=Y)
#' 
predict.msma = function(object, newX, newY=NULL,...){
##### Requirement #####
if(missing(newX)) stop("newX should be specified")
#if(missing(newY)) stop("newY should be specified")

if(is.null(newY) | (object$dmode == "PCA")) newY = newX

if(class(newX) != "list") newX = list(newX)
if(class(newY) != "list") newY = list(newY)

if(class(newX) != "list") stop("newX should be a list")
if(class(newY) != "list") stop("newY should be a list")

Ynb = length(newY)
Xnb = length(newX)

########## Start ##########
#newX = lapply(1:length(newX), function(x) scale(newX[[x]], center=object$Xscale[[x]]$`scaled:center`, scale=object$Xscale[[x]]$`scaled:scale`))
newX = lapply(1:length(newX), function(x) scale(newX[[x]], center=object$Xscale[[x]]$`scaled:center`))
newX = lapply(newX, function(x) apply(x, 2, function(x2){ x2[is.na(x2)] = 0; x2}))
sbX = lapply(1:Xnb, function(x) newX[[x]] %*% object$wbX[[x]] )
predX = lapply(1:Xnb, function(x) cbind(sbX[[x]]) %*% project(newX[[x]], sbX[[x]]))

if(!(is.null(newY) | object$dmode == "PCA")) {
#newY = lapply(1:length(newY), function(x) scale(newY[[x]], center=object$Yscale[[x]]$`scaled:center`, scale=object$Yscale[[x]]$`scaled:scale`))
newY = lapply(1:length(newY), function(x) scale(newY[[x]], center=object$Yscale[[x]]$`scaled:center`))
newY = lapply(newY, function(y) apply(y, 2, function(y2){ y2[is.na(y2)] = 0; y2}))
sbY = lapply(1:Ynb, function(x) newY[[x]] %*% object$wbY[[x]] )
predY = lapply(1:Ynb, function(x) cbind(sbY[[x]]) %*% project(newY[[x]], sbY[[x]]))
}else
{
predY = sbY = NULL
}

##### Output #####
list(X = predX, sbX = sbX, Y = predY, sbY = sbY)
}

#' Cross-Validation
#'
#' Cross-validated method to evaluate the fit of msma.
#'
#' k-fold cross-validation for \code{msma}. The evaluation is based on the matrix element-wise errors.
#'
#' @name cvmsma
#' @aliases cvmsma
#' @rdname cvmsma
#' @docType methods
#' @export
#'
#' @param X a matrix or list of matrices indicating the explanatory variable(s). This parameter is required.
#' @param Y a matrix or list of matrices indicating objective variable(s). This is optional. If there is no input for Y, then PCA is implemented.
#' @param Z a vector, response variable(s) for implementing the supervised version of (multiblock) PCA or PLS. This is optional. The length of Z is the number of subjects. If there is no input for Z, then unsupervised PLS/PCA is implemented.
#' @param comp numeric scalar for the number of components to be considered.
#' @param lambdaX numeric vector of regularized parameters for X, with a length equal to the number of blocks. If lambdaX is omitted, no regularization is conducted.
#' @param lambdaY numeric vector of regularized parameters for Y, with a length equal to the number of blocks. If lambdaY is omitted, no regularization is conducted.
#' @param eta numeric scalar indicating the parameter indexing the penalty family. This version contains only choice 1.
#' @param type a character, indicating the penalty family. In this version, only one choice is available: "lasso."
#' @param inX a vector or list of numeric vectors specifying the variables in X, always included in the model
#' @param inY a vector or list of numeric vectors specifying the variables in Y, always included in the model
#' @param muX a numeric scalar for the weight of X for the supervised case. 0 <= muX <= 1.
#' @param muY a numeric scalar for the weight of Y for the supervised case. 0 <= muY <= 1. 
#' @param nfold number of folds - default is 5.
#' @param seed seed number for the random number in the cross-validation.
#' @param intseed seed number for the random number in the parameter estimation algorithm.
#' 
#' @return \item{err}{The mean cross-validated errors which has three elements consisting of the mean of predict errors for X and Y, the errors for X and for Y in the PLS and only the errors for X in the PCA.}
#'
#' @examples
#'##### data #####
#'tmpdata = simdata(n = 50, rho = 0.8, Yps = c(10, 12, 15), Xps = 20, seed=1)
#'X = tmpdata$X; Y = tmpdata$Y 
#'
#'##### One Component CV #####
#'cv1 = cvmsma(X, Y, comp = 1, lambdaX=2, lambdaY=1:3, nfold=5, seed=1)
#'cv1
#'
#'##### Two Component CV #####
#'cv2 = cvmsma(X, Y, comp = 2, lambdaX=2, lambdaY=1:3, nfold=5, seed=1)
#'cv2
#'
cvmsma = function(X, Y=NULL, Z=NULL, comp=1, lambdaX, lambdaY=NULL, eta=1, type="lasso", inX=NULL, inY=NULL, muX = 0, muY = 0, nfold=5, seed=1, intseed=1){
##### Requirement #####
if(missing(X)) stop("X should be specified")
if(is.null(Y)){ 
dmode = "PCA"
Y = X
if(!is.null(lambdaY)) stop("lambdaY can not be specified")
}else{
dmode = "PLS"
}

if(class(X) != "list") X = list(X)
if(class(Y) != "list") Y = list(Y)
if(class(X) != "list") stop("X should be a list")
if(class(Y) != "list") stop("Y should be a list")

Ynb = length(Y)
Xnb = length(X)

if(!is.null(lambdaX)){if(length(lambdaX) != Xnb) stop("lambdaX should be same length as the number of blocks of X")}
if(!is.null(lambdaY)){if(length(lambdaY) != Ynb) stop("lambdaY should be same length as the number of blocks of Y")}

if(!is.null(inX)){if(class(inX) != "list") inX = list(inX)}
if(!is.null(inY)){if(class(inY) != "list") inY = list(inY)}

Yps = unlist(lapply(Y, ncol))
Xps = unlist(lapply(X, ncol))

if(!is.null(inX)){
extest = unlist(lapply(1:Xnb, function(a) inX[[a]] %in% (1:Xps[a]) ))
if(!(all(extest))) stop("inX should be specified the column number of X")
}
if(!is.null(inY)){
extest = unlist(lapply(1:Ynb, function(a) inY[[a]] %in% (1:Yps[a]) ))
if(!(all(extest))) stop("inY should be specified the column number of Y")
}

########## Start ##########
##### Set Up #####
n = nrow(Y[[1]])
id = 1:n
divid = cut(id, c(-Inf, floor(quantile(id, seq(1/nfold, 1-1/nfold, length=(nfold-1)))), Inf), labels = FALSE)
set.seed(seed); fold = sample(divid)

##### Cross-Validation #####
cv.err = sapply(1:nfold, function(i){

tmpX = lapply(X, function(tmpset) tmpset[!(fold == i),])
testX = lapply(X, function(testset) testset[(fold == i),])

if(dmode == "PCA"){
tmpY = testY = NULL
}else
{
tmpY = lapply(Y, function(tmpset) tmpset[!(fold == i),])
testY = lapply(Y, function(testset) testset[(fold == i),])
}

if(is.null(Z)){
tmpZ = testZ = NULL
}else{
if(class(Z) != "matrix") dim(Z) = c(length(Z),1)
tmpZ = Z[!(fold == i),]
testZ = Z[(fold == i),]
}

fit = msma(tmpX, tmpY, tmpZ, comp, lambdaX=lambdaX, lambdaY=lambdaY, eta=eta, type=type, inX=inX, inY=inY, muX = muX, muY = muY, intseed=intseed)
test = predict(fit, newX=testX, newY=testY)

#testX = lapply(1:length(testX), function(x) scale(testX[[x]], center=fit$Xscale[[x]]$`scaled:center`, scale=fit$Xscale[[x]]$`scaled:scale`))
testX = lapply(1:length(testX), function(x) scale(testX[[x]], center=fit$Xscale[[x]]$`scaled:center`))
Xerr = matserr(testX, test$X)
out = Xerr

if(dmode == "PLS"){
testY = lapply(1:length(testY), function(x) scale(testY[[x]], center=fit$Yscale[[x]]$`scaled:center`, scale=fit$Yscale[[x]]$`scaled:scale`))
testY = lapply(1:length(testY), function(x) scale(testY[[x]], center=fit$Yscale[[x]]$`scaled:center`))
Yerr = matserr(testY, test$Y)
out = c(mean(c(Xerr, Yerr)), Xerr, Yerr)
}

out
})

if(dmode == "PLS"){ 
err = rowMeans(cv.err)
names(err) = c("Average", "X", "Y")
}else
{
err = mean(cv.err)
}
err
##### Output #####
list(err = err)
}


#' Search for Number of Components
#'
#' Determination of the number of components based on cross-validated method or Bayesian information criterion (BIC)
#'
#' This function searches for the optimal number of components.
#'
#' @name ncompsearch
#' @aliases ncompsearch
#' @rdname ncompsearch
#' @docType methods
#' @export
#'
#' @param X a matrix or list of matrices indicating the explanatory variable(s). This parameter is required.
#' @param Y a matrix or list of matrices indicating objective variable(s). This is optional. If there is no input for Y, then PCA is implemented.
#' @param Z a vector, response variable(s) for implementing the supervised version of (multiblock) PCA or PLS. This is optional. The length of Z is the number of subjects. If there is no input for Z, then unsupervised PLS/PCA is implemented.
#' @param comps numeric vector for the candidates for the numbers of components to be selected.
#' @param lambdaX numeric vector of regularized parameters for X, with a length equal to the number of blocks. If lambdaX is omitted, no regularization is conducted.
#' @param lambdaY numeric vector of regularized parameters for Y, with a length equal to the number of blocks. If lambdaY is omitted, no regularization is conducted.
#' @param eta numeric scalar indicating the parameter indexing the penalty family. This version contains only choice 1.
#' @param type a character, indicating the penalty family. In this version, only one choice is available: "lasso."
#' @param inX a vector or list of numeric vectors specifying the variables in X, always included in the model
#' @param inY a vector or list of numeric vectors specifying the variables in Y, always included in the model
#' @param muX a numeric scalar for the weight of X for the supervised case. 0 <= muX <= 1.
#' @param muY a numeric scalar for the weight of Y for the supervised case. 0 <= muY <= 1. 
#' @param nfold number of folds - default is 5. 
#' @param x an object of class "\code{ncompsearch}", usually, a result of a call to \code{ncompsearch}
#' @param regpara logical, If TRUE, the regularized parameters search is also conducted simultaneously.
#' @param maxrep numeric scalar for the number of iteration.
#' @param minpct minimum candidate parameters defined as a percentile of automatically determined (possible) candidates.
#' @param maxpct maximum candidate parameters defined as a percentile of automatically determined (possible) candidates.
#' @param criterion a character, the evaluation criterion, "CV" for cross-validation, based on a matrix element-wise error, and "BIC" for Bayesian information criteria. The "BIC" is the default.
#' @param intseed seed number for the random number in the parameter estimation algorithm.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return \item{comps}{numbers of components}
#' @return \item{mincriterion}{minimum criterion value}
#' @return \item{criteria}{criterion values}
#' @return \item{optncomp}{optimal number of components with the minimum criteria value}
#' 
#' @examples
#' ##### data #####
#' tmpdata = simdata(n = 50, rho = 0.8, Yps = c(10, 12, 15), Xps = 20, seed=1)
#' X = tmpdata$X; Y = tmpdata$Y 
#' 
#' ##### number of components search #####
#' ncomp1 = ncompsearch(X, Y, comps = c(1, 5, 10*(1:5)), nfold=5)
#' plot(ncomp1)
#'
ncompsearch = function(X, Y=NULL, Z=NULL, comps=1:3, lambdaX=NULL, lambdaY=NULL, eta=1, type="lasso", inX=NULL, inY=NULL, muX = 0, muY = 0, nfold=5, regpara=FALSE, maxrep=3, minpct=0, maxpct=1, criterion=c("BIC", "CV")[1], intseed=1){

if(criterion=="BIC" & !regpara){
cves = msma(X, Y, Z, comp = max(comps), lambdaX=lambdaX, lambdaY=lambdaY, eta=eta, type=type, inX=inX, inY=inY, muX = muX, muY = muY, intseed=intseed)$bic
comps = 1:max(comps)
}else{
cv2 = lapply(comps, function(x) {
if(regpara){
cveout = regparasearch(X, Y, Z, eta, type, inX=inX, inY=inY, muX = muX, muY = muY, comp=x, nfold=nfold, maxrep=maxrep, minpct=minpct, maxpct=maxpct, criterion=criterion)
cve = cveout$mincriterion; lambdaX = cveout$optlambdaX; lambdaY = cveout$optlambdaY
}else
{
if(criterion=="CV"){
cve = cvmsma(X, Y, Z, comp=x, lambdaX=lambdaX, lambdaY=lambdaY, eta=eta, type=type, inX=inX, inY=inY, muX = muX, muY = muY, nfold=nfold, intseed=intseed)$err[1]
}else if(criterion=="BIC"){
cve = msma(X, Y, Z, comp=x, lambdaX=lambdaX, lambdaY=lambdaY, eta=eta, type=type, inX=inX, inY=inY, muX = muX, muY = muY, intseed=intseed)$bic[x]
}
}
list(cve = cve, lambdaX=lambdaX, lambdaY=lambdaY)
})
cves = unlist(lapply(cv2, function(x) x$cve))
}

#RRC = exp(diff(log(exp(diff(log(cves))))))

optidx = which.min(cves)
optncomp = comps[optidx]
mincriterion = cves[optidx]
if(criterion=="BIC" & !regpara){
optlambdaX = lambdaX; optlambdaY = lambdaY
}else{
optlambdaX = cv2[[optidx]]$lambdaX; optlambdaY = cv2[[optidx]]$lambdaY
}

out=list(criterion=criterion, comps=comps, mincriterion=mincriterion, criteria=cves, optncomp=optncomp, optlambdaX = optlambdaX, optlambdaY = optlambdaY)
class(out) = "ncompsearch"
out
}

#' @rdname ncompsearch
#' @method print ncompsearch
print.ncompsearch = function(x, ...)
{
cat("Optimal number of components: ")
cat(paste(x$optncomp, "(min CVE),"))
cat("\n")
}

#' @rdname ncompsearch
#' @method plot ncompsearch
#' @export
plot.ncompsearch = function(x,...){
#ylim = range(pretty(x$criteria))
ylab1 = ifelse(x$criterion == "CV", "CV error", "BIC")
plot(x$comps, x$criteria, type="b", xlab="Number of Components", ylab=ylab1,...)
abline(v=x$optncomp, lty=2, col=2)
abline(v=x$optncomp2, lty=2, col=3)
}


#' Regularized Parameters Search
#'
#' Regularized parameters search method for "msma". 
#'
#' This is a function for identifying the regularized parameters of sparseness lambdaX and lambdaY for \code{msma}. The initial range of candidates is computed based on fit, with regularized parameter values of 0. A binary search is conducted for dividing the parameter range into two regions. The representative value for the region is a median value, and the optimal region is selected using the minimum criteria obtained from the fit with that median value. The CV error or BIC can be used as criteria. The selected region is also divided into two region and the same process is iterated by maxrep times. Thus, the final median value in the selected region is set to be the optimal regularized parameter. The search is conducted with combinations of parameters for X and Y. The range of candidates for regularized parameters can be restricted, with a percentile of the limit (minimum or maximum) for the range.
#'
#' @name regparasearch
#' @aliases regparasearch
#' @rdname regparasearch
#' @docType methods
#' @export
#'
#' @param X a matrix or list of matrices indicating the explanatory variable(s). This parameter is required.
#' @param Y a matrix or list of matrices indicating objective variable(s). This is optional. If there is no input for Y, then PCA is implemented.
#' @param Z a vector, response variable(s) for implementing the supervised version of (multiblock) PCA or PLS. This is optional. The length of Z is the number of subjects. If there is no input for Z, then unsupervised PLS/PCA is implemented.
#' @param comp numeric scalar for the number of components to be considered.
#' @param eta numeric scalar indicating the parameter indexing the penalty family. This version contains only choice 1.
#' @param type a character, indicating the penalty family. In this version, only one choice is available: "lasso."
#' @param inX a vector or list of numeric vectors specifying the variables in X, always included in the model
#' @param inY a vector or list of numeric vectors specifying the variables in Y, always included in the model
#' @param muX a numeric scalar for the weight of X for the supervised case. 0 <= muX <= 1.
#' @param muY a numeric scalar for the weight of Y for the supervised case. 0 <= muY <= 1. 
#' @param nfold number of folds. Default is 5.
#' @param maxrep numeric scalar for the number of iterations.
#' @param minpct minimum candidate parameters defined as a percentile of automatically determined (possible) candidates.
#' @param maxpct maximum candidate parameters defined as a percentile of automatically determined (possible) candidates.
#' @param criterion a character, the evaluation criterion, "CV" for cross-validation, based on a matrix element-wise error, and "BIC" for Bayesian information criteria. The "BIC" is the default.
#' @param x an object of class "\code{regparasearch}", usually, a result of a call to \code{regparasearch}
#' @param intseed seed number for the random number in the parameter estimation algorithm.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return \item{optlambdaX}{Optimal parameters for X}
#' @return \item{optlambdaY}{Optimal parameters for Y}
#' @return \item{mincriterion}{Minimum criterion value}
#' @return \item{criteria}{All resulting criterion values in the process}
#' @return \item{pararange}{Range of candidates parameters}
#'
#' @examples
#' ##### data #####
#' tmpdata = simdata(n = 50, rho = 0.8, Yps = c(10, 12, 15), Xps = 20, seed=1)
#' X = tmpdata$X; Y = tmpdata$Y 
#' 
#' ##### Regularized parameters search #####
#' opt1 = regparasearch(X, Y, comp=1, nfold=5, maxrep=2)
#' opt1
#' fit4 = msma(X, Y, comp=1, lambdaX=opt1$optlambdaX, lambdaY=opt1$optlambdaY)
#' fit4
#' summary(fit4)
#'
#' ##### Restrict search range #####
#' opt2 = regparasearch(X, Y, comp=1, nfold=5, maxrep=2, minpct=0.5)
#' opt2
#'
regparasearch = function(X, Y=NULL, Z=NULL, eta=1, type="lasso", inX=NULL, inY=NULL, muX = 0, muY = 0, comp=1, nfold=5, maxrep=3, minpct=0, maxpct=1, criterion=c("BIC","CV")[1], intseed=1){
##### Requirement #####
if(missing(X)) stop("X should be specified")
if(minpct>=1 | minpct<0) stop("minpct should be 0 <= minpct < 1")
if(maxpct>1 | maxpct<=0) stop("maxpct should be 0 < maxpct <= 1")

##### Candidates for Regularized Parameters #####
intfit = msma(X, Y, Z, comp=1, lambdaX=rep(NULL, length(X)), lambdaY=rep(NULL, length(Y)), eta=eta, type=type, inX=inX, inY=inY, intseed=intseed)
if(!is.null(Y)){intweight = list(X=intfit$wbX, Y=intfit$wbY)}else{intweight = list(X=intfit$wbX)}
lambdaXYcands = lapply(intweight, function(y) do.call(rbind, lapply(y, function(x) cand4lambda(x, 2)[2])))

nblocks = unlist(lapply(lambdaXYcands, nrow))
if(!is.null(Y)){colidx = list(X=1:nblocks[1], Y=nblocks[1]+1:nblocks[2])}else{colidx = list(X=1:nblocks[1])}

totalnblocks = sum(nblocks)
idx = as.matrix(expand.grid(lapply(1:totalnblocks, function(x) c(2, 4))))

pararange = lapply(lambdaXYcands, function(y) apply(y, 1, function(y1) c(y1*minpct, y1*maxpct)))
range1s = do.call(cbind, pararange)

##### Iteration Start #####
criteria = list(); length(criteria) = maxrep
for(repidx in 1:maxrep){
quans = apply(range1s, 2, function(x){quantile(seq(x[1], x[2], length=10000))})
cands = apply(idx, 1, function(x){
unlist(lapply(colidx, function(y){if(length(y)==1){quans[x[y], y]}else{diag(quans[x[y], y])}}))
})

if(is.null(dim(cands))) dim(cands) = c(1, length(cands))

if(!is.null(Y)){
colnames(quans) = rownames(cands) = c(paste("lambdaX", 1:nblocks[1], sep=""), paste("lambdaY", 1:nblocks[2], sep=""))
}else
{
colnames(quans) = rownames(cands) = c(paste("lambdaX", 1:nblocks[1], sep=""))
}

##### CV for number of Regularized Parameters #####
cv3 = apply(cands, 2, function(x){
if(!is.null(Y)){
if(criterion=="CV"){
cvmsma(X, Y, Z, comp = comp, lambdaX=x[colidx$X], lambdaY=x[colidx$Y], eta=eta, type=type, inX=inX, inY=inY, muX = muX, muY = muY, nfold=nfold, intseed=intseed)$err[1]
}else if(criterion=="BIC"){
msma(X, Y, Z, comp = comp, lambdaX=x[colidx$X], lambdaY=x[colidx$Y], eta=eta, type=type, inX=inX, inY=inY, muX = muX, muY = muY, intseed=intseed)$bic[comp]
}
}else
{
if(criterion=="CV"){
cvmsma(X, Y, Z, comp = comp, lambdaX=x[colidx$X], lambdaY=NULL, eta=eta, type=type, inX=inX, inY=inY, muX = muX, muY = muY, nfold=nfold, intseed=intseed)$err[1]
}else if(criterion=="BIC"){
msma(X, Y, Z, comp = comp, lambdaX=x[colidx$X], lambdaY=NULL, eta=eta, type=type, inX=inX, inY=inY, muX = muX, muY = muY, intseed=intseed)$bic[comp]
}
}
})

##### #####
optidx = which.min(cv3)

criteria[[repidx]] = rbind(cands, criterion=cv3)

ranidx = sapply(idx[optidx,], function(x) c(x-1, x+1)); 
range1s = sapply(1:totalnblocks, function(x) quans[ranidx[,x], x])

}

##### #####
optidx = which.min(unlist(lapply(criteria, function(x) min(x["criterion",]))))
cands = criteria[[optidx]]; 
optidx = which.min(cands["criterion",])
mincriterion = cands["criterion",optidx]
cands = cands[!(rownames(cands) %in% "criterion"), ]
if(is.null(dim(cands))) dim(cands) = c(1, length(cands))
if(!is.null(Y)){
optlambda = list(X=cands[colidx$X, optidx], Y=cands[colidx$Y, optidx])
}else
{
optlambda = list(X=cands[colidx$X, optidx])
}
names(criteria) = paste("Step", 1:length(criteria))

##### #####
out=list(criterion=criterion, optlambdaX = optlambda$X, optlambdaY = optlambda$Y, mincriterion=mincriterion, criteria=criteria, pararange=pararange)
class(out) = "regparasearch"
out
}

#' @rdname regparasearch
#' @method print regparasearch
print.regparasearch = function(x, ...)
{
cat("Optimal parameters: \n")
cat("\n")
print(round(x$optlambdaX, 3))
cat("\n")
if(!is.null(x$optlambdaY)){
print(round(x$optlambdaY, 3))
cat("\n")
}
}

#' Parameters Search
#'
#' Combined method for optimizing the number of components and regularized parameters for "msma". 
#'
#' A function for identifying the regularized sparseness parameters lambdaX and lambdaY and the number of components for \code{msma}. Four search methods are available. The "simultaneous" method identifies the number of components by searching the regularized parameters in each component. The "regpara1st" identifies the regularized parameters by fixing the number of components, then searching for the number of components with the selected regularized parameters. The "ncomp1st" method identifies the number of components with a regularized parameter of 0, then searches for the regularized parameters with the selected number of components. The "regparaonly" method searches for the regularized parameters with a fixed number of components.
#'
#' @name optparasearch
#' @aliases optparasearch
#' @rdname optparasearch
#' @docType methods
#' @export
#'
#' @param X a matrix or list of matrices indicating the explanatory variable(s). This parameter is required.
#' @param Y a matrix or list of matrices indicating objective variable(s). This is optional. If there is no input for Y, then PCA is implemented.
#' @param Z a vector, response variable(s) for implementing the supervised version of (multiblock) PCA or PLS. This is optional. The length of Z is the number of subjects. If there is no input for Z, then unsupervised PLS/PCA is implemented.
#' @param search.method a character indicationg search methods, see Details. Default is "simultaneous".
#' @param comp numeric scalar for the number of components to be considered or the maximum canditate number of components.
#' @param eta numeric scalar indicating the parameter indexing the penalty family. This version contains only choice 1.
#' @param type a character, indicating the penalty family. In this version, only one choice is available: "lasso."
#' @param inX a vector or list of numeric vectors specifying the variables in X, always included in the model
#' @param inY a vector or list of numeric vectors specifying the variables in Y, always included in the model
#' @param muX a numeric scalar for the weight of X for the supervised case. 0 <= muX <= 1.
#' @param muY a numeric scalar for the weight of Y for the supervised case. 0 <= muY <= 1. 
#' @param nfold number of folds - default is 5. 
#' @param maxrep numeric scalar for the number of iteration.
#' @param minpct minimum candidate parameters defined as a percentile of automatically determined (possible) candidates.
#' @param maxpct maximum candidate parameters defined as a percentile of automatically determined (possible) candidates.
#' @param criterion a character, the evaluation criterion, "CV" for cross-validation, based on a matrix element-wise error, and "BIC" for Bayesian information criteria. The "BIC" is the default.
#' @param x an object of class "\code{optparasearch}", usually, a result of a call to \code{optparasearch}
#' @param intseed seed number for the random number in the parameter estimation algorithm.
#' @param ... further arguments passed to or from other methods.
#' 
#' @return \item{optncomp}{Optimal number of components}
#' @return \item{optlambdaX}{Optimal parameters for X}
#' @return \item{optlambdaY}{Optimal parameters for Y}
#' @return \item{mincriterion}{Minimum criterion value}
#' @return \item{criteria}{All resulting criterion values in the process}
#' @return \item{pararange}{Range of candidates parameters}
#'
#' @examples
#' ##### data #####
#' tmpdata = simdata(n = 50, rho = 0.8, Yps = c(10, 12, 15), Xps = 20, seed=1)
#' X = tmpdata$X; Y = tmpdata$Y 
#' 
#' ##### Regularized parameters search #####
#' opt1 = optparasearch(X, Y, search.method = "regparaonly", comp=1, nfold=5, maxrep=2)
#' opt1
#' fit4 = msma(X, Y, comp=opt1$optncomp, lambdaX=opt1$optlambdaX, lambdaY=opt1$optlambdaY)
#' fit4
#' summary(fit4)
#'
#' ##### Restrict search range #####
#' opt2 = optparasearch(X, Y, comp=3, nfold=5, maxrep=2, minpct=0.5)
#' opt2
#'
optparasearch = function(X, Y=NULL, Z=NULL, search.method = c("simultaneous", "regpara1st", "ncomp1st", "regparaonly")[1], eta=1, type="lasso", inX=NULL, inY=NULL, muX = 0, muY = 0, comp=1, nfold=5, maxrep=3, minpct=0, maxpct=1, criterion=c("BIC","CV")[1], intseed=1){

comps1 = 1:comp

##### both (regpara first) #####
if(search.method == "regpara1st"){
regpara1 = regparasearch(X=X, Y=Y, Z=Z, comp=comp, muX = muX, muY = muY, nfold=nfold, minpct=minpct, maxpct=maxpct, maxrep=maxrep, criterion=criterion)
params = ncompsearch(X=X, Y=Y, Z=Z, comps = comps1, lambdaX=regpara1$optlambdaX, lambdaY=regpara1$optlambdaY, muX = muX, muY = muY, nfold=nfold, regpara=FALSE, criterion=criterion)
##### both (ncomp first) #####
}else if(search.method == "ncomp1st"){
#params = ncompsearch(X=X, Y=Y, Z=Z, comps = min(c(min(unlist(lapply(X,dim))), min(unlist(lapply(Y,dim))))), muX = muX, muY = muY, nfold=nfold1, regpara=FALSE, criterion=criterion)
params = ncompsearch(X=X, Y=Y, Z=Z, comps = comps1, muX = muX, muY = muY, nfold=nfold, regpara=FALSE, criterion=criterion)
regpara1 = regparasearch(X=X, Y=Y, Z=Z, comp=params$optncomp, muX = muX, muY = muY, nfold=nfold, minpct=minpct, maxpct=maxpct, maxrep=maxrep, criterion=criterion)
params$optlambdaX = regpara1$optlambdaX
params$optlambdaY = regpara1$optlambdaY
##### both (simultaneous, pct=1/ncomp) #####
}else if(search.method == "simultaneous"){
params = ncompsearch(X=X, Y=Y, Z=Z, comps = comps1, muX = muX, muY = muY, nfold=nfold, regpara=TRUE, maxrep=maxrep, minpct=minpct, maxpct=maxpct, criterion=criterion)
##### regpara only (ncomp prespecified) #####
}else if(search.method == "regparaonly"){
params = list(optncomp = comp)
regpara1 = regparasearch(X=X, Y=Y, Z=Z, comp=params$optncomp, muX = muX, muY = muY, nfold=nfold, minpct=minpct, maxpct=maxpct, maxrep=maxrep, criterion=criterion)
params$optlambdaX = regpara1$optlambdaX
params$optlambdaY = regpara1$optlambdaY
params$cverrs = regpara1$cverrs
params$mincverr = regpara1$mincverr
}

params$search.method = search.method
class(params) = "optparasearch"

params
}

#' @rdname optparasearch
#' @method print optparasearch
print.optparasearch = function(x, ...)
{
cat("Optimal number of components: ")
cat(paste(x$optncomp, "(min CVE),"))
cat("\n")
cat("Optimal parameters: \n")
cat("\n")
print(round(x$optlambdaX, 3))
cat("\n")
if(!is.null(x$optlambdaY)){
print(round(x$optlambdaY, 3))
cat("\n")
}
}


#' Internal functions
#'
#' These are internal functions for \code{msma}
#' 
#' These are not intended for use by users.
#'
#' @rdname msma-internal
#' @keywords internal
msma_OneComp = function(X, Y, Z=NULL, lambdaX=NULL, lambdaY=NULL, eta=1, type="lasso", inX=NULL, inY=NULL, muX = 0, muY = 0, dmode = "PLS", verbose=FALSE, intseed=1){

Ynb = length(Y); Xnb = length(X)

if(any(is.null(lambdaX))){ lambdaX = rep(0, Xnb)}
if(any(is.null(lambdaY))){ lambdaY = rep(0, Ynb)}

if(dmode == "PCA"){ lambdaY=lambdaX}

if(is.null(Z)) Z = rep(0, nrow(X[[1]]))
muXY = 1 - muX - muY
if(muXY < -1.0e-16) stop(paste("Your muX and muY were", muX, "and", muY, ". muX or muY should be smaller values to be 0 <= 1 - muX - muY <= 1."))

#svd0 = svd(do.call(cbind, Y)); ssY = normvec(svd0$u[,1])
#ssY = normvec(rep(1, nrow(Y[[1]])))
set.seed(intseed); ssY = normvec(rnorm(nrow(X[[1]])))

if(length(X) == 1){wsX = 1}else{ set.seed(intseed); wsX = normvec(rnorm(length(X)))}
if(length(Y) == 1){wsY = 1}else{ set.seed(intseed); wsY = normvec(rnorm(length(Y)))}

##### Iteration #####
itr = 1
repeat{

### Xside (block) ###
XwsX = lapply(1:length(X), function(x) X[[x]] * wsX[x]) #XwsX = X
argus = lapply(XwsX, function(x){ tmp = c(t(x) %*% (muXY * ssY + muX * Z)); tmp[is.na(tmp)] = 0; tmp })
wbX = lapply(1:Xnb, function(x) normvec(sparse(normvec(argus[[x]]), lambdaX[x], eta, type, inX[[x]])))
sbX = lapply(1:Xnb, function(x) c(X[[x]] %*% wbX[[x]]) )
sbX2 = do.call(cbind, sbX)

### Xside (super) ###
wsX = normvec(t(sbX2) %*% (muXY * ssY + muX * Z))
ssX = sbX2 %*% wsX

### Yside (block) ###
YwsY = lapply(1:length(Y), function(x) Y[[x]] * wsY[x]) #YwsY = Y
argvs = lapply(YwsY, function(x){ tmp = c(t(x) %*% (muXY * ssX + muY * Z)); tmp[is.na(tmp)] = 0; tmp})
wbY = lapply(1:Ynb, function(y) normvec(sparse(normvec(argvs[[y]]), lambdaY[y], eta, type, inY[[y]])))
sbY = lapply(1:Ynb, function(x) c(Y[[x]] %*% wbY[[x]]) )
sbY2 = do.call(cbind, sbY)

### Yside (super) ###
wsY = normvec(t(sbY2) %*% (muXY * ssX + muY * Z))
ssY = sbY2 %*% wsY

### Convergence ###
if(itr > 2){
#dif = drop(crossprod(ssX - oldssX))
difX = drop(crossprod(ssX - oldssX))
difY = drop(crossprod(ssY - oldssY))
dif = max(difX, difY)
if(verbose) print(c(itr, dif))
if(dif < 0.0001 | itr > 20){break}
}
oldssX = ssX
oldssY = ssY
itr = itr + 1
}

### Output ###
list(wbX=wbX, sbX=sbX, wbY=wbY, sbY=sbY, ssX=ssX, wsX=wsX, ssY=ssY, wsY=wsY)
}

#' @rdname msma-internal
normvec = function(a){if(all(a == 0)){a}else{ c(a) / sqrt(drop(crossprod(c(a))))}}

#' @rdname msma-internal
sparse = function(x, lam, eta=1, type="lasso", inidx=NULL){
if(all(x == 0)){ lam = 0 }else if(lam >= max(abs(x))){ 
tmpdiff = diff(sort(unique(abs(x)))); tmpdiff = tmpdiff[tmpdiff > 1e-8]
lam = max(abs(x)) - tmpdiff[length(tmpdiff)]/2
}
lam = rep(lam, length(x))
if(!is.null(inidx)) lam[inidx] = 0
sout = sapply(1:length(x), function(x1) {
if(type=="lasso"){
sign(x[x1]) * max(c(abs(x[x1]) - lam[x1], 0))
}else 
if(type=="hard"){
x[x1] * (abs(x[x1]) > lam[x1])
}
})
sout
}

#' @rdname msma-internal
cand4lambda = function(x, len){
max1 = max(abs(x)) - diff(sort(abs(x)))[length(x)-1]/2
seq(0, max1, length=len)
}

#' @rdname msma-internal
searchseq = function(orgseq, opt){
idx = which(orgseq == opt)
if(idx == 1){x0 = orgseq[idx]; x1 = orgseq[idx+1]; p1 = c(0.25, 0.5)}else
if(idx == length(orgseq)){x0 = orgseq[idx-1]; x1 = orgseq[idx]; p1 = c(0.5, 0.75)}else
{x0 = orgseq[idx-1]; x1 = orgseq[idx+1]; p1 = c(0.25, 0.75)}
quantile(seq(x0, x1, length=10000), probs = p1)
}

#' @rdname msma-internal
project = function(x, y){
if(is.null(ncol(y)))
{c(x %*% cbind(y) / drop(crossprod(y)))}else
{
crossy = t(y) %*% y
i = 0; judge = TRUE
while(judge){
invcrossy = try(solve(crossy + 1e-18*10^i * diag(ncol(crossy))), silent = TRUE)
judge = (class(invcrossy) == "try-error")
i = i + 1
}
invcrossy %*% t(y) %*% x
}
}

#' @rdname msma-internal
matserr = function(X1, X2) mean(sapply(1:length(X1), function(x) mean((X1[[x]]-X2[[x]])^2)))

#' Generate Test Data Sets
#'
#' This is a function for generating multiblock data based on the multivariable normal distribution
#' 
#' The output is a list of matrices.
#'
#' @name simdata
#' @aliases simdata
#' @rdname simdata
#' @docType methods
#' @export
#'
#' @param n a numeric scalar for sample size.
#' @param rho a numeric scalar. Correlation coefficient for all matrices.
#' @param Yps a numeric vector indicating the numbers of columns for Y. The length of the vector corresponds to the number of blocks.
#' @param Xps a numeric vector indicating the numbers of columns for X. The length of the vector corresponds to the number of blocks.
#' @param seed a seed number for generating random numbers for reproducibility. Should be changed in an iterative study.
#'
#' @return \item{X}{Simulated X, which has a list form}
#' @return \item{Y}{Simulated Y, which has a list form}
simdata = function(n = 100, rho = 0.8, Yps = c(100, 120, 150), Xps = 500, seed=1){

if(length(n) > 1) stop("n should be scalar")
if(length(rho) > 1) stop("rho should be scalar")

#require("mvtnorm")

set.seed(seed)
Y = lapply(Yps, function(p){
sigma = matrix(rho, p, p)
diag(sigma) = 1
rmvnorm(n, rep(0, p), sigma)
})

X = lapply(Xps, function(p){
sigma = matrix(rho, p, p)
diag(sigma) = 1
rmvnorm(n, rep(0, p), sigma)
})

list(X=X, Y=Y)
}

#' @rdname msma-internal
strsimdata = function(n = 100, ncomp=5, Xps = 10, Yps = FALSE, rho=0.8, Z = FALSE, seed=1, minpct = 0.25, maxpct = 0.75){
#tmp = list(n = 100, ncomp=3, Yps = c(100, 120, 150), Xps = 10, rho=0.8, Z = FALSE, seed=1, minpct = 0.25, maxpct = 0.75); attach(tmp)

##### Generate X #####
idxZeroX = lapply(Xps, function(p){ 
tmp = lapply(1:ncomp, function(x) sort(sample(p, p*runif(1, minpct, maxpct))))
allzero = as.numeric(names(which(sort(table(unlist(tmp))) == ncomp)))
k = which.max(unlist(lapply(tmp, length)))
tmp[[k]] = tmp[[k]][!(tmp[[k]] %in% allzero)]
tmp
})

nZeroX = lapply(idxZeroX, function(x) unlist(lapply(x, length)))

SUX = matrix(rnorm(n*ncomp), n, ncomp)
UX = lapply(Xps, function(p) SUX/Xps)
VX = lapply(1:length(Xps), function(x){ v2=do.call(rbind, lapply(1:ncomp, function(y){ v = rnorm(Xps[[x]]); v[idxZeroX[[x]][[y]]]=0; v})); rownames(v2) = paste0("comp",1:ncomp); v2})

X = lapply(1:length(Xps), function(x) UX[[x]] %*% ginv(t(VX[[x]])))

##### Generate Y #####
if(!Yps[1]){idxZeroY = NULL; nZeroY = NULL; Y=NULL; VY=NULL}else
{
idxZeroY = lapply(Yps, function(p){ 
tmp = lapply(1:ncomp, function(x) sample(p, p*runif(1, minpct, maxpct)))
allzero = as.numeric(names(which(sort(table(unlist(tmp))) == ncomp)))
k = which.max(unlist(lapply(tmp, length)))
tmp[[k]] = tmp[[k]][!(tmp[[k]] %in% allzero)]
tmp
})

nZeroY = lapply(idxZeroY, function(x) unlist(lapply(x, length)))

SUY = apply(SUX, 2, function(x) rnorm(n, rho*x, sqrt(1-rho^2)))
UY = lapply(Yps, function(p) matrix(rnorm(n*ncomp), n, ncomp))
VY = lapply(1:length(Yps), function(x) do.call(rbind, lapply(1:ncomp, function(y){ v = rnorm(Yps[[x]]); v[idxZeroY[[x]][[y]]]=0; v})))

Y = lapply(1:length(Yps), function(x) UY[[x]] %*% ginv(t(VY[[x]])))
}

##### Generate Z #####
if(Z){
UX2 = do.call(cbind, UX)
ZcoefX = rnorm(ncol(UX2), ncol(UX2):1, 0.01)
UX2b = UX2 %*% ZcoefX
UY2b = 0
if(Yps[1]){
UY2 = do.call(cbind, UY)
ZcoefY = rnorm(ncol(UY2), ncol(UY2):1, 0.01)
UY2b = UY2 %*% ZcoefY
}
UXYb = UX2b + UY2b
UXYb = scale(UXYb, scale=FALSE)
Zprob = exp(UXYb)/(1+exp(UXYb))
Z = sapply(Zprob, function(x) rbinom(1, 1, x))
}

list(X=X, Y=Y, Z=Z, nZeroX=nZeroX, idxZeroX=idxZeroX, nZeroY=nZeroY, idxZeroY=idxZeroY, VX=VX, VY=VY)
}

#' @rdname msma-internal
ginv=function(X, tol = sqrt(.Machine$double.eps))
{
## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}
