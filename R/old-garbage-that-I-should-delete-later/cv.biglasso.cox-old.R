cv.biglasso.cox.old <- function(x, y, lambda, nfolds, grouped = F, parallel = F, ...) {
  # NOTE: The 'grouped' argument replicates the functionality of the 'grouped'
  # argument in cv.glmnet for Cox models. Specifically, if grouped = F then it
  # calculates the cross-validation error by the (weighted) mean of the test data
  # partial-likelihoods using the coefficients from the training data. Meanwhile, 
  # if grouped = T then it calculates the error as the (weighted) mean of the
  # difference between the full-data (train+test) and training data 
  # partial-likelihoods, both using the coefficients obtained from the training data.
  # (I believe the ladder is the Verweij & Van Houwelingen method... 
  #  see https://arxiv.org/pdf/1905.10432.pdf for a comparison [Basic vs. V&VH])
  #
  n <- nrow(y)
  
  # calculate the biglasso model over the sequence of (default) tuning parameters lambda
  if (missing(lambda)) {
    mod.init <- biglasso(X = x, y = y, family = "cox", ...)
    lambda.init <- mod.init$lambda
  } else {
    mod.init <- biglasso(X = x, y = y, lambda = lambda, family = "cox", ...)
    lambda.init <- lambda
  }
  
  # assign observations to cross-validation folds from 1, 2, ..., nfolds
  folds <- sample(cut(1:n, breaks = nfolds, labels = F))
  
  # arrange indices into a list of indicators corresponding to whether data belongs
  # to the test or training set in each of the cross-validation folds
  train.idx <- lapply(1:nfolds, function(i) !(folds %in% i))

  ##### meat of the CV algorithm is here #####
  train.fun <- function(trn) {
    tst <- !trn # test indices given current training fold
    
    # non-contiguous subsetting (via sub.big.matrix() is not currently supported) and so
    # we must create new big.matrix objects
    x.trn <- as.big.matrix(x[trn,,drop=F]); y.trn <- y[trn,,drop=F]
    x.tst <- as.big.matrix(x[tst,,drop=F]); y.tst <- y[tst,,drop=F]

    # compute training model and calculate training error
    mod.trn <- biglasso(X = x.trn, y = y.trn, lambda = lambda.init, family = "cox", ...)
    
    # replicate how Cox model cross-validation is done in glmnet
    wt <- sum(y.tst[,"status"]) # this could be divide-by-zero if the test set is too small (with too many censored observations)
    if (grouped) {
      plfull   <- coxnet.deviance(x = x, y = y, beta = mod.trn$beta)
      plminusk <- coxnet.deviance(x = x.trn, y = y.trn, beta = mod.trn$beta)
      loss.raw <- plfull - plminusk
    } else {
      plk      <- coxnet.deviance(x = x.tst, y = y.tst, beta = mod.trn$beta)
      loss.raw <- plk
    }
    return (list("loss.tst"  = loss.raw/wt,
                 "weight"    = wt))
  }
  if (parallel) {
    numCores <- detectCores()
    cv.loss <- mclapply(train.idx, train.fun, mc.cores = numCores)
    # these next two lines is a whole bunch of work just to replicate the format that
    # I had from the non-parallel case. There is almost certainly a more efficient way 
    # of moving forwards
    cv.loss <- unlist(cv.loss, recursive = F, use.names = F)
    cv.loss <- array(cv.loss, dim = c(2, nfolds), dimnames = list(c("loss.tst", "weight"), NULL))    
  } else {
    cv.loss <- sapply(train.idx, train.fun)
  }
  
  # extract cross-validation loss values and organize into a set of matrices (nlambda x nfolds)
  cv.loss.tst <- do.call("cbind", cv.loss["loss.tst",])
  weights     <- unlist(cv.loss["weight",])
  
  # average over all folds where each entry corresponds to a unique tuning value lambda
  test.loss   <- apply(cv.loss.tst, 1, weighted.mean, w = weights)
  # standard errors for test loss
  test.loss.se <- sqrt(apply(scale(t(cv.loss.tst), test.loss, FALSE)^2, 2, weighted.mean,
                             w = weights)/(nfolds - 1))
  test.loss.hi <- test.loss + test.loss.se
  test.loss.lo <- test.loss - test.loss.se
  
  best.lambdas <- getmin.lambda(mod.init$lambda, test.loss, test.loss.se)  
  
  out <- list()
  out$model.fit    <- mod.init
  out$test.loss    <- test.loss
  out$test.loss.se <- test.loss.se
  out$test.loss.hi <- test.loss.hi
  out$test.loss.lo <- test.loss.lo
  out$lambda       <- lambda.init
  out$foldid       <- folds
  out <- c(out, best.lambdas)
  
  return (out)
}