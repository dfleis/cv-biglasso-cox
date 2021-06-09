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

NCOLS <- 5000 # number of columns to load (fewer for quicker tests)

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

set.seed(124)
penalty  <- "enet"
alpha    <- 0.5
lambda   <- exp(seq(-1, -4, length.out = 100))
nfolds   <- 10
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
plot(cvout.gn)
plot.compare.cv(cvout.bl, cvout.gn)
#plot(cvout.gn$cvm - cvout.bl$test.loss ~ lambda, type = 'l', log = 'x', lwd = 2)








