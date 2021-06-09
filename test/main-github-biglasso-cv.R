##
# Comparing glmnet & the (proposed) biglasso CV outputs
# for penalized Cox models
#
###
#=======================================#
#================ SETUP ================#
#=======================================#
#remove.packages("biglasso")
#devtools::install_github("dfleis/biglasso")
library(biglasso)
library(glmnet)
library(data.table) # fread() so I can read fewer columns while testing, otherwise read.csv is fine
options(datatable.fread.datatable=FALSE) # format the data as a data.frame instead of a data.table

set.seed(124)

NCOLS <- 500 # number of columns to load (fewer for quicker tests)

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


# no <- 1000
# nx <- 100
# dat <- coxed::sim.survdata(N = no, xvars = nx)
# X <- as.matrix(dat$xdata)
# y <- cbind("time" = as.numeric(dat$data$y), "status" = dat$data$failed)
# Xbig <- as.big.matrix(X)

#===========================================#
#================ RUN TESTS ================#
#===========================================#
# load custom R functions (mainly plot tools for the following tests)
lapply(list.files("./R/", full.names = T), source)


penalty  <- "enet"
alpha    <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
nfolds   <- 5
lambda   <- exp(seq(-2, -6, length.out = 100))
grouped  <- F # grouped = F is 'basic' CV loss, grouped = T is 'V&VH' CV loss. see https://arxiv.org/pdf/1905.10432.pdf
parallel <- F # "use parallel  to fit each fold. Must register parallel before hand, such as doMC or others"
ncores   <- 1
trace.it <- 1

# fit.bl <- biglasso(X       = Xbig,
#                    y       = y,
#                    family  = "cox",
#                    penalty = "enet",
#                    alpha   = alpha,
#                    lambda  = lambda)
# fit.gn <- glmnet(x      = X,
#                  y      = y,
#                  family = "cox",
#                  alpha  = alpha,
#                  lambda = lambda)

set.seed(124) # necessary for testing as the bigmemory package has an (unintended?) effect on random number generation
pt <- proc.time()
cv.bl <- cv.biglasso(X       = Xbig,                               
                     y       = y,
                     family  = "cox",
                     penalty = penalty,
                     alpha   = alpha,
                     nfolds  = nfolds,
                     lambda  = lambda,
                     grouped = grouped,
                     ncores  = ncores,
                     trace   = as.logical(trace.it))
proc.time() - pt


if (parallel) {doMC::registerDoMC(cores = ncores)}

pt <- proc.time()
cv.gn <- cv.glmnet(x        = X,                               
                   y        = y,
                   family   = "cox",
                   alpha    = alpha,
                   nfolds   = nfolds,
                   lambda   = lambda,
                   grouped  = grouped,
                   parallel = parallel,
                   trace.it = trace.it,
                   foldid   = cv.bl$cv.ind) 
proc.time() - pt

#plot(cv.bl)
#plot(cv.gn)
#plot(cv.bl$cve ~ lambda, log = 'x', type = 'l', lwd = 4, ylim = range(cv.bl$cve, cv.gn$cvm))
#lines(cv.gn$cvm ~ lambda, col = 'red', lty = 'dashed', lwd = 4)

#plot.cv(cv.bl$cve, cv.bl$cvse, lambda, nzero = apply(cv.bl$fit$beta != 0, 2, sum))
#plot.cv(cv.gn$cvm, cv.gn$cvsd, lambda, nzero = cv.gn$nzero)

plot.compare.cv2(cv.bl, cv.gn)

#fields::image.plot(1:ncol(X), -log(lambda), as.matrix(fit.bl$beta) - as.matrix(fit.gn$beta), col = pals::coolwarm(65))

















