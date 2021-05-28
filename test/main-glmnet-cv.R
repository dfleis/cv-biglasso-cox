##
# Fully self-contained script to preprocess the data & perform cross-validated 
# penalized Cox regression via glmnet::cv.glmnet().
#
# NOTES:
#   * I'm receiving convergence issues with the test data for some values of the
#     tuning parameter. Maybe this will be resolved using the full data? This
#     is also causing cv.glmnet to be extremely slow as it will exhaust a
#     maximum number of iterations before giving up without reaching the
#     desired threshold.
#   * The 'parallel' argument in glmnet::cv.glmnet() allows us to run the 
#     cross-validation folds in parallel. As such, I'm not sure having more cores 
#     than folds is beneficial? This is something worth testing.
#     WARNING: Increasing the cores will significantly increase the RAM usage!
#       * cores = 2 peaked at ~15GB of RAM (incl. baseline system usage) (600 seconds)
#       * cores = 4 peaked at ~20GB of RAM (incl. baseline) (480 seconds)
#   * The memory management in this script is currently a bit sloppy. We may note that
#     there are a few copies of the data are left lying around: The full data is loaded
#     into 'dat', then the first 14 columns are extracted and stored into X1, the final
#     18, 19, ..., ncol(dat) columns into are stored in X2, then both X1 and X2 are used
#     to create a single design matrix X = cbind(X1, X2). In principle we may want to 
#     avoid having 3 copies of the same data lying around, but given enough RAM it may
#     not be relevant.
#   * Keep an eye out for any errors originating from the usage of 'model.matrix' with
#     our large data. I've gotten some overflow errors while prototyping some of these
#     scripts.
#   * If we wish, glmnet allows us to cross-validate over the mixing parameter
#     concurrently with the tuning parameter. To do so, set 'relax = T' and specify 
#     'gamma' as a vector of mixing params up for consideration (by default, 
#     'gamma = c(0, 0.25, 0.5, 0.75, 1)').
##

#=======================================#
#================ SETUP ================#
#=======================================#
set.seed(124) # set a seed so the CV foldid's remain fixed for reproducibility
library(glmnet)
library(doMC) # library necessary for glmnet to run cv folds in parallel
library(data.table) # fread() so I can read fewer columns while testing, otherwise read.csv is fine
options(datatable.fread.datatable=FALSE) # format the data as a data.frame instead of a data.table

# # load partial data (NCOLS leading columns)
# NCOLS <- 500 # number of columns to load (fewer for quicker tests)
# pt <- proc.time()
# dat <- fread("./data/sample_SRTR_cleaned.csv", select = 1:NCOLS, # select columns (must contain columns 1-17)
#              na.strings = c("", "NA"), stringsAsFactors = T)
# proc.time() - pt
# load full data (>50000 columns)
pt <- proc.time() ## It takes ~60 seconds to load the full sample (2000 rows x 58000+ columns)
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

# outcome variable & censor indicator (it may be a good idea to convert to a survival::is.Surv object)
Y <- as.matrix(dat[,c("surv", "fail_dc")])
colnames(Y) <- c('time', 'status') # glmnet expects the responses to be labeled in this way

#===========================================#
#================ RUN TESTS ================#
#===========================================#
alpha        <- 0.5 # elastic net penalty (alpha = 0 is ridge and alpha = 1 is lasso)
nfolds       <- 10
lambda       <- exp(seq(-2, -4, length.out = 100))
type.measure <- "deviance" # c("deviance", "C") available error metrics for Cox models
grouped      <- T # grouped = F is 'basic' CV loss, grouped = T is 'V&VH' CV loss. see https://arxiv.org/pdf/1905.10432.pdf
parallel     <- F # "use parallel  to fit each fold. Must register parallel before hand, such as doMC or others"
trace.it     <- 1
# "If trace.it=1, then progress bars are displayed; useful for big 
# models that take a long time to fit. Limited tracing if parallel=TRUE"

if (parallel) {registerDoMC(cores = 2)}

pt <- proc.time()
cv.gn <- cv.glmnet(x            = X,                               
                   y            = Y,
                   family       = "cox",
                   alpha        = alpha,
                   nfolds       = nfolds,
                   lambda       = lambda,
                   type.measure = type.measure,
                   grouped      = grouped,
                   parallel     = parallel,
                   trace.it     = trace.it) 
proc.time() - pt

plot(cv.gn)


