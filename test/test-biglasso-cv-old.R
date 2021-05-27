##
# NOTES:
#   1) At no point do I take full advantage of the bigmemory::big.matrix data structure as
#      I load the entire data into memory anyway. To do... Load the data via read.big.matrix
#      and figure out how to construct the design matrix from there.
#      Note... The final 18, ..., ncol(dat) columns are in a indicator/binary format and so 
#      these variables are already in the format required by the design matrix (interactions
#      notwithstanding). The first 1, ..., 14 columns contain categorical/factor variables 
#      and so they must be converted into the appropriate format (sets of indicators). Fortunately
#      this is not a memory intensive task if we consider the first 14 columns independently.
#      We can construct the design matrix for the first 14 variables then figure out how to 
#      merge them (columnwise) with the other 18, ..., ncol(dat) predictors.
#
###
#=======================================#
#================ SETUP ================#
#=======================================#
library(parallel)
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

penalty <- "enet"
alpha   <- 0.5
lambda  <- exp(seq(-2, -4, length.out = 100))
nfolds  <- 5
grouped <- T

set.seed(124)
pt <- proc.time()
cvout.bl <- cv.biglasso.cox.old(x       = Xbig,
                            y       = y,
                            penalty = penalty,
                            alpha   = alpha,
                            nfolds  = nfolds,
                            grouped = grouped,
                            lambda  = lambda,
                            parallel = T)
proc.time() - pt
#lambda <- cvout.bl$lambda

set.seed(124)
pt <- proc.time()
cvout.gn <- cv.glmnet(x       = X,                               
                      y       = y,
                      family  = "cox",
                      alpha   = alpha,
                      grouped = grouped,
                      nfolds  = nfolds, # type.measure = "deviance",
                      lambda  = lambda,
                      foldid  = cvout.bl$foldid)
                      #gamma = seq(0, 1, length.out = 10),
                      #relax = T)
proc.time() - pt

plot.cv.biglasso.cox(cvout.bl)
plot(cvout.gn)
#M <- sapply(cvout.gn$relaxed$statlist, function(x) x$cvm)
#fields::image.plot(M)
#persp(M)


##### plot both CV outputs in a single figure
sign.lambda <- 1
xlab <- expression(log(lambda))
if (sign.lambda < 0) xlab <- expression(-log(lambda))

plot.args <- list(
  x    = sign.lambda * log(lambda),
  xlim = range(sign.lambda * log(lambda)),
  ylim = range(cvout.bl$test.loss.lo, cvout.bl$test.loss.hi,
               cvout.gn$cvlo, cvout.gn$cvup),
  xlab = xlab,
  ylab = "Partial Likelihood Deviance",
  type = "n"
)
do.call("plot", plot.args)
error.bars(sign.lambda * log(lambda),
           cvout.bl$test.loss.hi,
           cvout.bl$test.loss.lo,
           col = rgb(1, 0, 0, 0.6))
error.bars(sign.lambda * log(lambda),
           cvout.gn$cvup,
           cvout.gn$cvlo,
           col = rgb(0, 0, 1, 0.6))
points(sign.lambda * log(lambda)[seq(1, length(lambda), 2)], 
       cvout.bl$test.loss[seq(1, length(lambda), 2)], 
       pch = 20, col = 'red')
points(sign.lambda * log(lambda)[seq(1, length(lambda), 2)], 
       cvout.gn$cvm[seq(1, length(lambda), 2)], 
       pch = 17, col = 'blue')
legend("topright", legend = c("cv.biglasso.cox", "cv.glmnet"),
       pch = c(20, 17), col = c("red", "blue"))

# #### coefficient figures
# beta.bl <- cvout.bl$model.fit$beta
# beta.gn <- cvout.gn$glmnet.fit$beta
# beta.bl.mat <- as.matrix(beta.bl)
# beta.gn.mat <- as.matrix(beta.gn)
# 
# nz.bl <- apply(beta.bl, 2, function(x) sum(x != 0))
# nz.gn <- apply(beta.gn, 2, function(x) sum(x != 0))
# 
# plot(NA,
#      xlim = range(sign.lambda * log(lambda)),
#      ylim = range(nz.bl, nz.gn),
#      xlab = xlab,
#      ylab = "Nonzero Variables")
# grid(); abline(h = 0, lwd = 1.5, col = 'gray60')
# lines(x = sign.lambda * log(lambda), y = nz.bl, type = 'S', col = 'red', lty = "longdash", lwd = 2)
# lines(x = sign.lambda * log(lambda), y = nz.gn, type = 'S', col = 'blue', lty = "dotted", lwd = 2)
# legend("topright", legend = c("cv.biglasso.cox", "cv.glmnet"),
#        lty = c("longdash", "dotted"), col = c("red", "blue"), lwd = 2)
# 
# 
# plot(NA,
#      xlim = range(sign.lambda * log(lambda)),
#      ylim = range(beta.bl.mat, beta.gn.mat),
#      xlab = xlab,
#      ylab = "Coefficient Value")
# grid(); abline(h = 0, lwd = 1.5, col = 'gray60')
# for (i in 1:nrow(beta.bl.mat)) {
#   lines(beta.bl.mat[i,]~log(lambda), col = 'red')
#   lines(beta.gn.mat[i,]~log(lambda), col = 'blue', lty = 'dotted')
# }


