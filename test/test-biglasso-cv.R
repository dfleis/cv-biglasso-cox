##
#
# To do...
#
###
#=======================================#
#================ SETUP ================#
#=======================================#
library(biglasso)
library(glmnet)

### load data
filenameX <- "data/sample_SRTR_cleaned_X_NArm.csv" # predictors
filenameY <- "data/sample_SRTR_cleaned_Y_NArm.csv" # response

# predictors (already processed into its design matrix)
pt <- proc.time()
Xbig <- read.big.matrix(file.path(filenameX), header = T, type = "double") 
proc.time() - pt

# responses
Y <- read.csv(filenameY, header = T)
Y <- as.matrix(Y)
colnames(Y) <- c('time', 'status') # biglasso expects the responses to be labeled in this way




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






















