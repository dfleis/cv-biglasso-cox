##
# Run biglasso (without cross-validation). Note that in this script we 
# load the full data into RAM and construct the design matrix (as opposed
# to constructing the design matrix outside this script in a more
# memory-efficient way).
##

#=======================================#
#================ SETUP ================#
#=======================================#
library(biglasso)

# load data
pt <- proc.time() # takes ~60 seconds to load the sample data
dat <- read.csv("./data/sample_SRTR_cleaned.csv", na.strings = c("", "NA"), stringsAsFactors = T)
proc.time() - pt

###
### TEMPORARY(?): ignore NA/blank/missing cases
###
dat <- dat[complete.cases(dat),] # overwrite the original data with the data omitting any NA observations

# construct design matrix (excluding pID as this appears to be simply an identifier)
pt <- proc.time()
X1 <- model.matrix(~.-1, data = dat[,1:14]) # confounders (note: no intercept column in a Cox model's design matrix)
X2 <- as.matrix(dat[,18:ncol(dat)] * 1.0)  # other predictors
# the '* 1.0' above is just a hack to convert 'int' format to 'num' so that X1 and X2 are formatted the same ,
# otherwise cbind will form a data.frame and not a matrix
X <- cbind(X1, X2) # this should be a fully prepared design matrix
proc.time() - pt

Xbig <- as.big.matrix(X) # turn into a big.matrix for use in biglasso()
remove(X1, X2, X)

# outcome variable & censor indicator (it may be a good idea to convert to a survival::is.Surv object)
Y <- as.matrix(dat[,c("surv", "fail_dc")])
colnames(Y) <- c('time', 'status') # glmnet expects the responses to be labeled in this way

remove(dat)

#===========================================#
#================ RUN TESTS ================#
#===========================================#
penalty <- "enet"
alpha   <- 0.5
lambda  <- exp(seq(-1, -6, length.out = 100))
ncores  <- 8

pt <- proc.time()
fit.bl <- biglasso(X       = Xbig,
                   y       = Y,
                   family  = "cox",
                   penalty = penalty,
                   alpha   = alpha,
                   lambda  = lambda,
                   ncores  = ncores,
                   verbose = T)
proc.time() - pt



#plot(fit.bl)

nonzero.bl <- apply(fit.bl$beta, 2, function(b) sum(b != 0))
plot(NA,
     xlim = range(log(lambda)),
     ylim = range(nonzero.bl),
     xlab = expression(log(lambda)),
     ylab = "Number of Variables",
     main = "Number of Variables Retained\nby biglasso::biglasso()")
grid(); abline(h = 0, lwd = 1.5, col = 'gray60')
lines(x = log(lambda), y = nonzero.bl, type = 'S', col = 'red', lwd = 2)

