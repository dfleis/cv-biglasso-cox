##
# Testing how to use bigmemory's functions and data structures in order to
#   3) subset a big.matrix object in such a way to construct training and
#      testing cross-validation sets.
# 
# For these tests I assume that the data has already been preprocessed so that 
# the predictors (X) and responses (Y, survival times and censor indicators) 
# are in separate files and that the factor predictors have already been mapped
# to sets of indicators (i.e. X is already the correct design matrix).
#
# SUBSETTING STRATEGY (building training/test subsets)
#   * Look at the 'row.idx' argument in biglasso::biglasso()
##
library(biglasso)

## load data
filenameX <- "data/sample_SRTR_cleaned_X_NArm.csv" # predictors
filenameY <- "data/sample_SRTR_cleaned_Y_NArm.csv" # response
# predictors
pt <- proc.time()
Xbig <- read.big.matrix(file.path(filenameX), header = T, type = "double") 
proc.time() - pt

# responses
Y <- read.csv(filenameY, header = T)
Y <- as.matrix(Y)
colnames(Y) <- c('time', 'status')

### start subsetting test ###
lambda <- exp(seq(-1, -6, length.out = 100))
nfolds <- 5


n <- nrow(Y)

# assign observations to cross-validation folds from 1, 2, ..., nfolds
folds <- sample(cut(1:n, breaks = nfolds, labels = F))

# list of integer vectors containing the row indices of x that will be use
# to fit the model (making use of the 'row.idx' argument in biglasso::biglasso())
train.row.idx <- lapply(1:nfolds, function(i) which(folds != i))

### NOTE:
###  It turns out that I don't think we can directly use the row.idx argument 
###  since we will need a version of the subsetted data to calculate the deviance
###  values in each fold
###
Xbig.train.1 <- deepcopy(Xbig, rows = train.row.idx[[1]])

# compare outputs for the row.idx method vs. the bigmemory::deepcopy method
pt <- proc.time()
set.seed(124)
fit.fold.1.a <- biglasso(X       = Xbig,
                         y       = Y, 
                         row.idx = train.row.idx[[1]],
                         family  = "cox", 
                         penalty = "enet", 
                         alpha   = 0.5, 
                         lambda  = lambda)
proc.time() - pt

pt <- proc.time()
set.seed(124)
fit.fold.1.b <- biglasso(X       = Xbig.train.1,
                         y       = Y[train.row.idx[[1]],], 
                         family  = "cox", 
                         penalty = "enet", 
                         alpha   = 0.5, 
                         lambda  = lambda)
proc.time() - pt

# both coefficient estimates ought to be equal
sum(fit.fold.1.a$beta != fit.fold.1.b$beta)
sum((fit.fold.1.a$beta - fit.fold.1.b$beta)^2)




