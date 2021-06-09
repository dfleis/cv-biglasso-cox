##
# Comparing glmnet & the (proposed) biglasso CV outputs
# for penalized Cox models
#
###
#=======================================#
#================ SETUP ================#
#=======================================#
library(biglasso)
library(glmnet)
library(data.table) # fread() so I can read fewer columns while testing, otherwise read.csv is fine
options(datatable.fread.datatable=FALSE) # format the data as a data.frame instead of a data.table

NCOLS <- 1000 # number of columns to load (fewer for quicker tests)

# load data
pt <- proc.time()
dat <- fread("./data/sample_SRTR_cleaned.csv", select = 1:NCOLS, # select columns (must contain columns 1-17)
             na.strings = c("", "NA"), stringsAsFactors = T)
proc.time() - pt

# TEMPORARY: ignore NA/blank/missing cases
dat2 <- dat[complete.cases(dat),]
remove(dat) # remove the original data since it's quite large

# construct design matrix (excluding pID as this appears to be simply an identifier)
pt <- proc.time()
X1 <- model.matrix(~.-1, data = dat2[,1:14]) # confounders (note that we do not have an intercept column in a Cox model's design matirx)
X2 <- as.matrix(dat2[,18:ncol(dat2)] * 1.0)  # other variables/predictors
X <- cbind(X1, X2)
proc.time() - pt

# convert to big.matrix for biglasso
Xbig <- as.big.matrix(X)

# outcome variable & censor indicator (it may be a good idea to convert to a survival::is.Surv object)
y <- as.matrix(dat2[,c("surv", "fail_dc")])
colnames(y) <- c("time", "status")

#===========================================#
#================ RUN TESTS ================#
#===========================================#
# load custom R functions
lapply(list.files("./R/", full.names = T), source)

penalty  <- "enet"
alpha    <- 0.5
lambda   <- exp(seq(-2, -4, length.out = 100))
nfolds   <- 5
grouped  <- T
parallel <- F


pt <- proc.time()
cvout.bl <- cv.biglasso.cox(x        = Xbig,
                            y        = y,
                            penalty  = penalty,
                            alpha    = alpha,
                            lambda   = lambda,
                            nfolds   = nfolds,
                            grouped  = grouped,
                            parallel = parallel)
proc.time() - pt


pt <- proc.time()
cvout.gn <- cv.glmnet(x        = X,                               
                      y        = y,
                      family   = "cox",
                      alpha    = alpha,
                      nfolds   = nfolds,
                      lambda   = lambda,
                      grouped  = grouped,
                      parallel = parallel,
                      trace.it = 1,
                      foldid   = cvout.bl$foldid) 
proc.time() - pt





############ figures

plot.cv.biglasso.cox(cvout.bl)
plot(cv.gn)

plot(cvout.bl$test.loss ~ cvout.bl$lambda, type = 'l', log = 'x')
lines(cvout.gn$cvm ~ cvout.gn$lambda, col = 'red')

plot((cvout.bl$test.loss - cvout.gn$cvm) ~ cvout.bl$lambda, type = 'l', log = 'x')






####



















x <- Xbig
y <- Y


n <- nrow(y)

# calculate the biglasso model over the sequence of (default) tuning parameters lambda
#if (missing(lambda)) {
#  mod.init <- biglasso::biglasso(X = x, y = y, family = "cox", ...)
#  lambda.init <- mod.init$lambda
#} else {
mod.init <- biglasso::biglasso(X = x, y = y, lambda = lambda, family = "cox", penalty = penalty, alpha = alpha)
lambda.init <- lambda
#}

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
  mod.trn <- biglasso::biglasso(X = x.trn, y = y.trn, lambda = lambda.init, family = "cox", penalty = penalty, alpha = alpha)# ...)
  
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
pt <- proc.time()
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
proc.time() - pt

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





# works:
system.time(fit1 <- biglasso(X = Xbig, y = Y, screen = "SSR", family = "cox", penalty = "enet", alpha = 0.5))
# not defined:
system.time(fit2 <- biglasso(X = Xbig, y = Y, screen = "sscox", family = "cox", penalty = "lasso"))
#system.time(fit3 <- biglasso(X = Xbig, y = Y, screen = "scox", family = "cox", penalty = "enet", alpha = 0.5))
#system.time(fit4 <- biglasso(X = Xbig, y = Y, screen = "safe", family = "cox", penalty = "enet", alpha = 0.5))
# not available for elastic-net:
system.time(fit5 <- biglasso(X = Xbig, y = Y, screen = "Adaptive", family = "cox", penalty = "enet", alpha = 0.5))
# broken(?):
system.time(fit6 <- biglasso(X = Xbig, y = Y, screen = "Hybrid", family = "cox", penalty = "enet", alpha = 0.5))
# works (but presumably slower than "SSR"):
system.time(fit7 <- biglasso(X = Xbig, y = Y, screen = "None", family = "cox", penalty = "enet", alpha = 0.5))

plot(fit1)
plot(fit7)


alg.logistic = "Newton"

f <- function(alg.logistic = c("Newton", "MM")) {
  alg.logistic = match.arg(alg.logistic)
  (!identical(penalty, "lasso") || alg.logistic == "MM")
  #c(!identical(penalty, "lasso"), alg.logistic == "MM")
}
f()










