##
# Testing how to use bigmemory's functions and data structures in order to
#   1) load data without using a data.frame/matrix (via read.csv, read.table, etc.) 
#      as an intermediate,
#   2) take full advantage of biglasso's memory efficiency, and
#   3) subset a big.matrix object in such a way to construct training and
#      testing cross-validation sets.
##
library(bigmemory)
library(biglasso)

set.seed(124)

# #### TEST 1 ####
# # generate an arbitrary set of data
# n <- 11
# p <- 7
# X <- matrix(rnorm(n * p), nrow = n)
# colnames(X) <- paste0("Var", 1:p)
# # write to disk
# write.csv(X, "data/bigmemory_test_data.csv", row.names = F)
# # remove existing objects from memory
# remove(X, n, p)
# 
# # load data as big.matrix object
# filename <- "data/bigmemory_test_data.csv"
# Xbig <- read.big.matrix(file.path(filename), header = T)

##### TEST 2 ####

filename <- "data/sample_SRTR_cleaned.csv"
pt <- proc.time()
Xbig.raw <- read.big.matrix(file.path(filename), header = T)
proc.time() - pt

str(Xbig.raw[,1:14])







