##
# Compare the Cox regression implemented in biglasso::biglasso() 
# with that of glmnet::glmnet()
#
# NOTES:
#   1) biglasso is substantially quicker for large data (see https://github.com/YaohuiZeng/biglasso) 
#   2) I don't (yet) take advantage of the bigmemory::big.matrix structures for which 
#      biglasso was designed to use. At the moment I load the entire data into memory 
#      and run the model as I would if it were a regular matrix object. One of the selling 
#      features of biglasso is its memory efficiency, being able to run these models 
#      without having to load the data into RAM (via bigmemory).
#      TO DO... 
#       a) Look at the read.big.matrix() function and build the big.matrix() structure
#          without obliterating my machine's RAM.
#       b) Figure out how to construct the design matrix without first converting the 
#          data into a data.frame or matrix object, i.e. build the design matrix directly
#          from the big.matrix object.
##

#=======================================#
#================ SETUP ================#
#=======================================#
library(biglasso)
library(glmnet)
library(data.table) # fread() so I can read fewer columns while testing, otherwise read.csv is fine
options(datatable.fread.datatable=FALSE) # format the data as a data.frame instead of a data.table

NCOLS <- 100 # number of columns to load (fewer for quicker tests)

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

penalty <- "enet" # biglasso requires us to explicitly set penalty = "enet" even if alpha is between (0, 1)
alpha   <- 0.5
lambda  <- exp(seq(0, -5, length.out = 100))

pt <- proc.time()
fit.bl <- biglasso(X       = Xbig,
                   y       = y,
                   family  = "cox",
                   penalty = penalty,
                   alpha   = alpha,
                   lambda  = lambda)
proc.time() - pt

pt <- proc.time()
fit.gn <- glmnet(x       = X,                               
                 y       = y,
                 family  = "cox",
                 alpha   = alpha,
                 lambda  = lambda)
proc.time() - pt




fit.bl$beta
fit.gn$beta
nonzero.bl <- apply(fit.bl$beta, 2, function(b) sum(b != 0))
nonzero.gn <- apply(fit.gn$beta, 2, function(b) sum(b != 0))

plot(NA,
     xlim = range(log(lambda)),
     ylim = range(nonzero.bl, nonzero.gn),
     xlab = expression(log(lambda)),
     ylab = "Number of Variables",
     main = "Number of Variables Retained\nvia Penalized Cox Models")
grid(); abline(h = 0, lwd = 1.5, col = 'gray60')
lines(x = log(lambda), y = nonzero.bl, type = 'S', col = 'red', lty = "longdash", lwd = 2)
lines(x = log(lambda), y = nonzero.gn, type = 'S', col = 'blue', lty = "dotted", lwd = 2)
legend("topright", legend = c("biglasso", "glmnet"),
       lty = c("longdash", "dotted"), col = c("red", "blue"), lwd = 2)




