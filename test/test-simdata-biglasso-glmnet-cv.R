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

set.seed(124)
n <- 100
p <- 1000
p_nz <- 0.1
beta <- rnorm(p, 0, 1) * rbinom(p, 1, p_nz)
dat <- coxed::sim.survdata(N = n, T = 100000, xvars = p, beta = beta)
X <- as.matrix(dat$xdata)
y <- cbind("time" = as.numeric(dat$data$y), "status" = dat$data$failed)
Xbig <- as.big.matrix(X)

#===========================================#
#================ RUN TESTS ================#
#===========================================#
# load custom R functions (mainly plot tools for the following tests)
lapply(list.files("./R/", full.names = T), source)

penalty  <- "enet"
alpha    <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
nfolds   <- 5
lambda   <- exp(seq(1, -1, length.out = 100))
grouped  <- T # grouped = F is 'basic' CV loss, grouped = T is 'V&VH' CV loss. see https://arxiv.org/pdf/1905.10432.pdf
parallel <- F # "use parallel  to fit each fold. Must register parallel before hand, such as doMC or others"
ncores   <- 1
trace.it <- 1

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

