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

plot.compare.cv <- function(cv.obj, sign.lambda = 1, ...) {
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

