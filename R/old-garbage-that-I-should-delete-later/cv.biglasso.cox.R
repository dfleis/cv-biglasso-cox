cv.biglasso.cox <- function(x, y, lambda, nfolds, grouped = F, parallel = F, ...) {
  # NOTE: The 'grouped' argument replicates the functionality of the 'grouped'
  # argument in cv.glmnet for Cox models. Specifically, if grouped = F then it
  # calculates the cross-validation error by the (weighted) mean of the test data
  # partial-likelihoods using the coefficients from the training data. Meanwhile, 
  # if grouped = T then it calculates the error as the (weighted) mean of the
  # difference between the full-data (train+test) and training data 
  # partial-likelihoods, both using the coefficients obtained from the training data.
  # (I believe the ladder is the Verweij & Van Houwelingen method... 
  #  see https://arxiv.org/pdf/1905.10432.pdf for a comparison ["Basic" vs. "V&VH"])
  n <- nrow(y)
  
  # calculate the biglasso model over the sequence of (default) tuning parameters lambda
  if (missing(lambda)) {
    mod.init <- biglasso::biglasso(X = x, y = y, family = "cox", ...)
    lambda.init <- mod.init$lambda
  } else {
    mod.init <- biglasso::biglasso(X = x, y = y, lambda = lambda, family = "cox", ...)
    lambda.init <- lambda
  }
  
  # assign observations to cross-validation folds from 1, 2, ..., nfolds
  folds <- sample(cut(1:n, breaks = nfolds, labels = F))
  
  # list of integer vectors containing the row indices of x that will be use
  # to fit the model (making use of the 'row.idx' argument in biglasso::biglasso())
  train.row.idx <- lapply(1:nfolds, function(i) which(folds != i))
  
  
  ##### meat of the CV algorithm is here #####
  train.fun <- function(trn.rows) {
    tst.rows <- !((1:n) %in% trn.rows)
    
    x.trn <- bigmemory::deepcopy(x = x, rows = trn.rows)
    x.tst <- bigmemory::deepcopy(x = x, rows = tst.rows)
    y.trn <- y[trn.rows,,drop=F]
    y.tst <- y[tst.rows,,drop=F]
    
    # compute training model and calculate training error
    mod.trn <- biglasso::biglasso(X = x.trn, y = y.trn, lambda = lambda.init, family = "cox", ...)
    
    ### replicate Cox model cross-validation loss is computed in glmnet::glmnet()
    # calculate the weight associated with the k-th CV fold (done in glmnet())
    # NOTE: this could be lead to a division-by-zero if test set is too small (when all response
    # observations are censored in the test set)
    wt <- sum(y.tst[,"status"]) 
    if (grouped) { # "V&VH cross-validation error"
      plfull   <- glmnet::coxnet.deviance(x = x, y = y, beta = mod.trn$beta)
      plminusk <- glmnet::coxnet.deviance(x = x.trn, y = y.trn, beta = mod.trn$beta)
      loss.raw <- plfull - plminusk
    } else { # "basic cross-validation error" 
      plk      <- glmnet::coxnet.deviance(x = x.tst, y = y.tst, beta = mod.trn$beta)
      loss.raw <- plk
    }
    return (list("loss.tst"  = loss.raw/wt,
                 "weight"    = wt))
  }
  if (parallel) {
    numCores <- detectCores()
    cv.loss <- mclapply(train.row.idx, train.fun, mc.cores = numCores)
    # these next two lines is a whole bunch of work just to replicate the format that
    # I had from the non-parallel case. There is almost certainly a more efficient way 
    # of moving forwards
    cv.loss <- unlist(cv.loss, recursive = F, use.names = F)
    cv.loss <- array(cv.loss, dim = c(2, nfolds), dimnames = list(c("loss.tst", "weight"), NULL))    
  } else {
    cv.loss <- sapply(train.row.idx, train.fun)
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