##
# The biglasso package expects the data to be formatted as a big.matrix object.
# Such objects are memory-efficient matrices and as a matrix object (as opposed 
# to a data.frame) all entries are assumed to be the same kind of data type (i.e.
# all numeric). 
#
# Columns 16 and 17 contain Cox model responses (survival times and 
# censor indicators) and columns 18,...,ncol(data) columns are binary {0, 1} 
# predictors, both of these columns are fine as-is and can be left unprocessed.
# However, columns 1, ..., 14 contain a confounding predictors, many of which are
# multi-level factors/categorical variables (stored as character variables -- 
# the 'matrix' type can't contain factors).
#
# For this reason we cannot use bigmemory::read.big.matrix directly as it will
# convert any non-numeric data to NA. Therefore, the first 14 columns must first
# be processed into a corresponding set of indicators (column 13 is a continuous 
# random variable and would be fine to leave as-is, but for simplicity it's likely
# easier to just let it tag along). 

# My strategy is to twofold:
#   1) Split the predictors and response into separate files, and
#   2) load the first 14 columns separately, do the appropriate processing, 
#      and recombine with the final 18, ..., ncol(data) columns to form
#      the predictor data file.
#
# Note that by splitting off the first 14 columns the preprocessing takes just a few
# seconds even if the number of rows is >500k.
##
library(bigmemory)
#library(biganalytics) # from the same authors as bigmemory, includes many useful functions?
library(data.table) # fread() so I can read fewer columns while testing, otherwise read.csv is fine
options(datatable.fread.datatable=FALSE) # format the data as a data.frame instead of a data.table

### load first 14 covariates
filename <- "data/sample_SRTR_cleaned.csv"
pt <- proc.time()
datX1 <- fread("./data/sample_SRTR_cleaned.csv", select = 1:14,
              na.strings = c("", "NA"), stringsAsFactors = T)
proc.time() - pt

# create the design matrix (indicator sets) for the first 14 columns
options(na.action="na.pass") # force model.matrix to retain rows with NA data
X1 <- model.matrix(~.-1, data = datX1)
options(na.action="na.omit")

# write to disk
write.csv(X1, "./data/sample_SRTR_cleaned_X1.csv", row.names = F)

### load response variables
pt <- proc.time()
datY <- fread("./data/sample_SRTR_cleaned.csv", select = 16:17,
              na.strings = c("", "NA"), stringsAsFactors = T)
proc.time() - pt

# write to disk
write.csv(Y, "./data/sample_SRTR_cleaned_Y.csv", row.names = F)

####
#### TO DO... What is the best way of splitting the final 18, ..., ncol(data)
#### columns without loading the entire thing into memory? For our 2000 sample
#### this is fine, but when dealing with 50000 observations the data occupies
#### something like (50000 rows * 50000 cols) * 8 bytes per cell ~= 20 gigabytes
####















