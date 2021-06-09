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

