error.bars <- function(x, upper, lower, width = 0.01, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  #range(upper, lower)
}
plot.cv.biglasso.cox <- function(cv.obj, sign.lambda = 1, ...) {
  xlab <- expression(log(lambda))
  if (sign.lambda < 0) xlab <- expression(-log(lambda))
  
  plot.args <- list(
    x = sign.lambda * log(cv.obj$lambda),
    y = cv.obj$test.loss,
    ylim = range(cv.obj$test.loss.lo, cv.obj$test.loss.hi),
    xlab = xlab,
    ylab = "Partial Likelihood Deviance",
    type = "n"
  )
  new.args <- list(...)
  
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cv.obj$lambda),
             cv.obj$test.loss.hi,
             cv.obj$test.loss.lo,
             col = 'gray60')
  points(sign.lambda * log(cv.obj$lambda), cv.obj$test.loss,
         pch = 20, col = 'red')
  axis(side   = 3,
       at     = sign.lambda * log(cv.obj$lambda),
       labels = apply(cv.obj$model.fit$beta, 2, function(z) sum(z != 0)),
       tick   = F,
       line   = 0)
  abline(v = sign.lambda * log(cv.obj$lambda.min), lty = 3)
  abline(v = sign.lambda * log(cv.obj$lambda.1se), lty = 3)
  invisible()
}
plot.compare.cv <- function(cv.bl, cv.gn, sign.lambda = 1, ...) {
  ### NOTE: This only works for the all-in-R version I wrote in cv-biglasso-cox
  ### See plot.compare.cv2 for a version which works with the cv.biglasso output
  # draw glmnet and biglasso cv plots ontop of one another for easier comparison
  # assumes that the tuning parameter used over the CV procedure is identical
  xlab <- expression(log(lambda))
  if (sign.lambda < 0) xlab <- expression(-log(lambda))
  
  lambda <- cv.bl$lambda
  
  plot.args <- list(
    x    = sign.lambda * log(lambda),
    xlim = range(sign.lambda * log(lambda)),
    ylim = range(cv.bl$test.loss.lo, cv.bl$test.loss.hi,
                 cv.gn$cvlo, cv.gn$cvup),
    xlab = xlab,
    ylab = "Partial Likelihood Deviance",
    type = "n"
  )
  do.call("plot", plot.args)
  grid()
  error.bars(sign.lambda * log(lambda),
             cv.bl$test.loss.hi,
             cv.bl$test.loss.lo,
             col = rgb(1, 0, 0, 0.6))
  error.bars(sign.lambda * log(lambda),
             cv.gn$cvup,
             cv.gn$cvlo,
             col = rgb(0, 0, 1, 0.6))
  points(sign.lambda * log(lambda)[seq(1, length(lambda), 2)], 
         cv.bl$test.loss[seq(1, length(lambda), 2)], 
         pch = 20, col = 'red')
  points(sign.lambda * log(lambda)[seq(1, length(lambda), 2)], 
         cv.gn$cvm[seq(1, length(lambda), 2)], 
         pch = 17, col = 'blue')
  abline(v = log(cv.bl$lambda.min), col = 'red', lty = 'dashed')
  abline(v = log(cv.bl$lambda.1se), col = 'red', lty = 'dashed')
  abline(v = log(cv.gn$lambda.min), col = 'blue', lty = 'dotted')
  abline(v = log(cv.gn$lambda.1se), col = 'blue', lty = 'dotted')
  
  axis(side   = 3,
       at     = sign.lambda * log(cv.bl$lambda),
       labels = apply(cv.bl$model.fit$beta, 2, function(z) sum(z != 0)),
       tick   = F,
       line   = 0.5, col.axis = 'red')
  axis(side   = 3,
       at     = sign.lambda * log(cv.gn$lambda),
       labels = cv.gn$nzero,
       tick   = F,
       line   = -0.5, col.axis = 'blue')
  
  
  legend("topright", legend = c("cv.biglasso.cox", "cv.glmnet"),
         pch = c(20, 17), col = c("red", "blue"), lty = c('dashed', 'dotted'), seg.len = 2)
  
  invisible()
}


plot.cv <- function(cvm, cvsd, lambda, nzero = NA, sign.lambda = 1, ...) {
  # reproduce glmnet figures for cross-validation outputs
  # but instead of supplying the library-specific object 
  # (e.g. "cv.glmnet"), supply the cross-validated mean 
  # error (cvm) and standard deviation (cvsd)
  cvlo <- cvm - cvsd
  cvup <- cvm + cvsd
  
  ## NOTE: this function is written in another file and ought to be loaded
  ## It's probably best if I do this procedure within this function (or
  ## write the function in this file).
  minlams <- getmin.lambda(lambda, cvm, cvsd)
  
  xlab <- expression(log(lambda))
  if (sign.lambda < 0) xlab <- expression(-log(lambda))
  
  plot.args <- list(
    x = sign.lambda * log(lambda),
    y = cvm,
    ylim = range(cvlo, cvup),
    xlab = xlab,
    ylab = "Partial Likelihood Deviance", # assumes Cox model for label
    type = "n"
  )
  new.args <- list(...)
  
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(lambda), cvup, cvlo, col = 'gray60')
  points(sign.lambda * log(lambda), cvm, pch = 20, col = 'red')
  axis(side   = 3,
       at     = sign.lambda * log(lambda),
       labels = nzero, #apply(cv.obj$model.fit$beta, 2, function(z) sum(z != 0)),
       tick   = F,
       line   = 0)
  abline(v = sign.lambda * log(minlams$lambda.min), lty = 3)
  abline(v = sign.lambda * log(minlams$lambda.1se), lty = 3)
  invisible()
}

plot.compare.cv2 <- function(cv.bl, cv.gn, sign.lambda = 1, ...) {
  # same as plot.compare.cv but works with cv.biglasso output
  xlab <- expression(log(lambda))
  if (sign.lambda < 0) xlab <- expression(-log(lambda))
  
  cv.bl.lo <- cv.bl$cve - cv.bl$cvse
  cv.bl.up <- cv.bl$cve + cv.bl$cvse
  
  lambda.bl <- cv.bl$lambda
  lambda.gn <- cv.gn$lambda
  ## NOTE: this function is written in another file and ought to be loaded
  ## It's probably best if I do this procedure within this function (or
  ## write the function in this file).
  minlams.bl <- getmin.lambda(lambda.bl, cv.bl$cve, cv.bl$cvse) # glmnet already provides this in the output
  
  plot.args <- list(
    x    = sign.lambda * log(lambda.gn), # I think using the glmnet output is safer since it will not omit saturated lambdas
    xlim = range(sign.lambda * log(lambda.bl), sign.lambda * log(lambda.gn)),
    ylim = range(cv.bl.lo, cv.bl.up, cv.gn$cvlo, cv.gn$cvup),
    xlab = xlab,
    ylab = "Partial Likelihood Deviance", # assumes Cox model for label
    type = "n"
  )
  do.call("plot", plot.args)
  grid()
  error.bars(sign.lambda * log(lambda.bl), cv.bl.up, cv.bl.lo, col = rgb(1, 0, 0, 0.6))
  error.bars(sign.lambda * log(lambda.bl), cv.gn$cvup, cv.gn$cvlo, col = rgb(0, 0, 1, 0.6))
  points(sign.lambda * log(lambda.bl)[seq(1, length(lambda.bl), 2)], 
         cv.bl$cve[seq(1, length(lambda.bl), 2)], 
         pch = 20, col = 'red')
  points(sign.lambda * log(lambda.gn)[seq(1, length(lambda.gn), 2)], 
         cv.gn$cvm[seq(1, length(lambda.gn), 2)], 
         pch = 17, col = 'blue')
  abline(v = log(minlams.bl$lambda.min), col = 'red', lty = 'dashed')
  abline(v = log(minlams.bl$lambda.1se), col = 'red', lty = 'dashed')
  abline(v = log(cv.gn$lambda.min), col = 'blue', lty = 'dotted')
  abline(v = log(cv.gn$lambda.1se), col = 'blue', lty = 'dotted')
  
  axis(side   = 3,
       at     = sign.lambda * log(cv.bl$fit$lambda),
       labels = apply(cv.bl$fit$beta != 0, 2, sum),
       tick   = F,
       line   = 0.5, col.axis = 'red')
  axis(side   = 3,
       at     = sign.lambda * log(lambda.gn),
       labels = cv.gn$nzero,
       tick   = F,
       line   = -0.5, col.axis = 'blue')
  
  legend("topright", legend = c("cv.biglasso.cox", "cv.glmnet"),
         pch = c(20, 17), col = c("red", "blue"), lty = c('dashed', 'dotted'), seg.len = 2)
  
  invisible()
}

