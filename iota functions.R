### Extends a vector. Only works correctly if distances between values are fixed.
# x is the vector
# n is the number of values to be added
# Returns extended vector.
extendVec <- function(x, n) {
  d <- x[length(x) - 1] - x[length(x)]
  e <- x[length(x)] - (1:n) * d
  c(x, e)
}

### Uses extendVec() to extend the lambda vector of glmnet.
# x is the lambda vector of glmnet
# n is the number of lambdas to be added
# Returns an extended lambda vector.
GLMnetLambdas <- function(x, n) exp(extendVec(log(x), n))

### Extract coefs from a cv.glmnet() fit
# fit is the cv.glmnet() fit
# Returns coefficient matrix.
GLMnetExtractCoefs <- function(fit) t(do.call(cbind, lapply(coef(fit, s = fit$lambda.min), matrix)))

### Reparametrize a coefficient matrix of glmnet to a Reference Category (Multinomial) model.
# b is the coefficient matrix
# Returns a reparametrized matrix.
GLMnetAsRefCat <- function(b) b - matrix(b[1, ], nrow = nrow(b), ncol = ncol(b), byrow = T)

### Brier score
# p is the predicted probability matirx
# y is the matrix containing the indicators of disease presence.
# Returns Brier score.
Brier <- function(p, y) c("Brier" = mean(colMeans((p - y) ^ 2)))

# Compute performance measures.
# y is the observed outcome vector
# yi is the matrix of observed disease indicators.
# p is the predicted probability matrix
# d is an optional factor for rounding digits
# Returns a data.frame with 2 calibration slopes, the pdi and the Brier score.
Performance <- function(y, yi, p, d = c(2, 3, 4)) {
  data.frame(c(round(PolCalAll(y, p, ref = 3)$slopes[c(1,4)], d[1]), 
               round(PdiDF(y, p, cats = colnames(p))[4], d[2]),
               round(Brier(p, yi), d[3]) ) )
}