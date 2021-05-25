##
# Given a sequence of tuning parameters (lambda), a sequence of 
# mean loss measures obtained via cross-validation, and a 
# sequence of standard errors of loss measures obtained via
# cross-validation, return the tuning parameter (and its index)
# associated with the least loss as well as the tuning parameter
# (and its index) associated with the 1-standard-error-rule
##
getmin.lambda <- function(lambda, cv.mean, cv.se) {
  cv.min      <- min(cv.mean, na.rm = T)
  cv.mean.idx <- cv.mean <= cv.min
  
  lambda.min <- max(lambda[cv.mean.idx], na.rm = T)
  lambda.min.idx <- match(lambda.min, lambda)
  
  cv.1se.min <- (cv.mean + cv.se)[lambda.min.idx]
  cv.1se.idx <- cv.mean <= cv.1se.min
  lambda.1se <- max(lambda[cv.1se.idx], na.rm = T)
  lambda.1se.idx <- match(lambda.1se, lambda)
  
  list(
    "lambda.min"     = lambda.min,
    "lambda.min.idx" = lambda.min.idx,
    "lambda.1se"     = lambda.1se,
    "lambda.1se.idx" = lambda.1se.idx)
}