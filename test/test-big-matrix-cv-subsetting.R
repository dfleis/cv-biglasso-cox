##
# Testing how to use bigmemory's functions and data structures in order to
#   3) subset a big.matrix object in such a way to construct training and
#      testing cross-validation sets.
# 
# For these tests I assume that the data has already been preprocessed so that 
# the predictors (X) and responses (Y, survival times and censor indicator) 
# are in separate files and that the factor predictors have already been mapped
# to sets of indicators (i.e. X is already the correct design matrix).
##
library(biglasso)

# strategy: Look at the 'row.idx' argument in biglasso::biglasso()

##### Try to run biglasso by directly loading the design matrix as a big.matrix object
filenameX <- "data/sample_SRTR_cleaned_X_NArm.csv"
filenameY <- "data/sample_SRTR_cleaned_Y_NArm.csv"
pt <- proc.time()
Xbig <- read.big.matrix(file.path(filenameX), header = T, type = "double")
proc.time() - pt
