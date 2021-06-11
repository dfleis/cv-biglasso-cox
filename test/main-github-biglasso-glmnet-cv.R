###############################################################
# REAL DATA 
#
# Comparing glmnet & the (proposed) biglasso CV outputs
# for penalized Cox models
#
###############################################################
#=======================================#
#================ SETUP ================#
#=======================================#
#remove.packages("biglasso")
#devtools::install_github("dfleis/biglasso")
library(biglasso)
library(glmnet)
library(data.table) # fread() so I can read fewer columns while testing, otherwise read.csv is fine
options(datatable.fread.datatable=FALSE) # format the data as a data.frame instead of a data.table

# load custom R functions (mainly plot tools for the following tests)
source("~/projects/cv-biglasso-cox/R/plot.cv.biglasso.cox.R")
source("~/projects/cv-biglasso-cox/R/getmin.lambda.R")

NCOLS <- 1000 # number of columns to load (fewer for quicker tests)
# load data
filepath <- "./data/sample_SRTR_cleaned.csv"

pt <- proc.time()
if (!is.null(NCOLS)) {
  # select columns (NOTE: must contain columns 1-17)
  dat <- fread(filepath, select = 1:NCOLS, na.strings = c("", "NA"), stringsAsFactors = T)
} else {
  dat <- fread(filepath, na.strings = c("", "NA"), stringsAsFactors = T)
}
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

set.seed(124) # necessary for testing as the bigmemory package has an (unintended?) effect on random number generation
penalty  <- "enet"
alpha    <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
nfolds   <- 10
lambda   <- exp(seq(-2, -4, length.out = 100))
grouped  <- T # grouped = F is 'basic' CV loss, grouped = T is 'V&VH' CV loss. see https://arxiv.org/pdf/1905.10432.pdf
parallel <- T # "use parallel  to fit each fold. Must register parallel before hand, such as doMC or others"
ncores   <- 10
trace.it <- 1

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
(tm.bl <- proc.time() - pt)

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
(tm.gn <- proc.time() - pt)

#===========================================#
#================= FIGURES =================#
#===========================================#



beta.bl <- cv.bl$fit$beta
beta.gn <- cv.gn$glmnet.fit$beta
beta.bl.cv <- cv.bl$fit$beta[,which(cv.bl$lambda == cv.bl$lambda.min)]
beta.gn.cv <- cv.gn$glmnet.fit$beta[,which(cv.gn$lambda == cv.gn$lambda.min)]

### PLOT 2
myclrs <- c(rgb(26/255, 133/255, 255/255, 1), 
            rgb(212/255, 17/255, 89/255, 0.6))
plot(NA, xlim = c(1, ncol(X)), ylim = range(beta.bl.cv, beta.gn.cv), 
     ylab = "Coef. Value", xlab = "Covariate Index",
     main = "Coefficient Values at the Minimal Lambda")
grid(); abline(h = 0, lwd = 2, col = 'gray50'); abline(v = ncol(X1) + 0.5)
for (i in 1:ncol(X)) {
  if (beta.bl.cv[i] != 0) segments(x0 = i, x1 = i, y0 = 0, y1 = beta.bl.cv[i], lwd = 4, col = myclrs[1])
  if (beta.gn.cv[i] != 0) segments(x0 = i, x1 = i, y0 = 0, y1 = beta.gn.cv[i], lwd = 4, col = myclrs[2])
}
legend("topleft", legend = c("biglasso", "glmnet"), col = myclrs, seg.len = 2, lwd = 4)

### PLOT 3
myclrs2 <- c(rgb(26/255, 133/255, 255/255, 0.5), 
             rgb(212/255, 17/255, 89/255, 0.5))
plot(NA, xlab = expression(log(lambda)), ylab = "Fitted Coef.",
     xlim = range(log(lambda)), ylim = range(as.matrix(beta.gn), as.matrix(beta.bl)),
     main = "Elastic-Net Coefficient Paths")
grid(); abline(h = 0, lwd = 2, col = 'gray50')
for (i in 1:ncol(X)) {
  lines(beta.bl[i,] ~ log(lambda), col = myclrs2[1])
  lines(beta.gn[i,] ~ log(lambda), col = myclrs2[2])
}
legend("topright", legend = c("biglasso", "glmnet"), col = myclrs2, seg.len = 2, lwd = 4)

### PLOT 1
plot.compare.cv2(cv.bl, cv.gn)


tm.bl
tm.gn
