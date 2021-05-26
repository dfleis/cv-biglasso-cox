##
# Testing how to use bigmemory's functions and data structures in order to
#   1) load data without using a data.frame/matrix (via read.csv, read.table, etc.) 
#      as an intermediate,
#   2) take full advantage of biglasso's memory efficiency, and
#   3) subset a big.matrix object in such a way to construct training and
#      testing cross-validation sets. (DONE IN ANOTHER FILE: See test-big-matrix-cv.R)
# 
# For these tests I assume that the data has already been preprocessed so that 
# the predictors (X) and responses (Y, survival times and censor indicator) 
# are in separate files and that the factor predictors have already been mapped
# to sets of indicators (i.e. X is already the correct design matrix).
##
library(bigmemory)
library(biglasso)

set.seed(124)

# ##### TEST 1 #####
# # generate an arbitrary set of data
# n <- 11
# p <- 7
# X <- matrix(rnorm(n * p), nrow = n)
# colnames(X) <- paste0("Var", 1:p)
# # write to disk
# write.csv(X, "data/bigmemory_test_data.csv", row.names = F)
# # remove existing objects from memory
# remove(X, n, p)
# 
# # load data as big.matrix object
# filename <- "data/bigmemory_test_data.csv"
# Xbig <- read.big.matrix(file.path(filename), header = T)

##### TEST 2 #####
##### Try to run biglasso by directly loading the design matrix as a big.matrix object
filenameX <- "data/sample_SRTR_cleaned_X_NArm.csv"
filenameY <- "data/sample_SRTR_cleaned_Y_NArm.csv"
pt <- proc.time()
Xbig <- read.big.matrix(file.path(filenameX), header = T, type = "double")
proc.time() - pt

# I think Y can be loaded without resorting to big.matrix objects (this is convenient since
# I believe biglasso only permits X to be a big.matrix)
Y <- read.csv(filenameY, header = T)

Y <- as.matrix(Y) # biglasso expects the data to be a two-column matrix
colnames(Y) <- c("time", "status") # biglasso requires the data to be named as such

# compute Cox models
alpha  <- 0.1
lambda <- exp(seq(0, -6, length.out = 1e2))

# run biglasso
pt <- proc.time() # takes ~12 seconds to run the full data (NA observations omitted)
fit.bl <- biglasso(X       = Xbig,
                   y       = Y,
                   family  = "cox",
                   penalty = "enet",
                   alpha   = alpha,
                   lambda  = lambda)
proc.time() - pt 

# compare with glmnet to verify everything is fine
X <- as.matrix(Xbig)

pt <- proc.time() # takes ~19 seconds to run the full data (NA omitted)
fit.gn <- glmnet::glmnet(x       = X,
                         y       = Y,
                         family  = "cox",
                         alpha   = alpha,
                         lambda  = lambda)
proc.time() - pt



nonzero.bl <- apply(fit.bl$beta, 2, function(z) sum(z != 0))
nonzero.gn <- apply(fit.gn$beta, 2, function(z) sum(z != 0))
plot(nonzero.gn ~ lambda, log = 'x', type = 's')
lines(nonzero.bl ~ lambda, type = 's', col = 'red', lty = 'dashed')
# biglasso uses a screening algorithm (something it calls 'SSR') to remove speed up
# calculations by discarding predictors that it estimates will not be selected

#### OLD CODE BELOW... recycling bin/graveyard

xdesc <- bigmemory::describe(Xbig.raw)
#XX <- bigmemory::attach.big.matrix(xdesc)


Xdesc <- describe(Xbig.raw)
XX <- attach.big.matrix(Xbig.raw)
Xdesc@description$nrow
Xdesc@description$ncol

col.idx <- 1:Xdesc@description$nrow %in% 1:14
col.idx


Xdesc@address
str(Xdesc)
  
  
  
  
  
