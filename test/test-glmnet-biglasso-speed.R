##
# Test how glmnet::glmnet() and biglasso::biglasso() scale in 
# performance (speed) as the number of observations/rows increase
##
#=======================================#
#================ SETUP ================#
#=======================================#
library(glmnet)
library(biglasso)

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

# convert Xbig to a matrix object for glmnet
X <- as.matrix(Xbig)

Y <- rbind(Y, Y, Y)
X <- rbind(X, X, X)
Xbig <- as.big.matrix(X)

#X <- X[,1:25000]
#Xbig <- as.big.matrix(X)

#===========================================#
#================ RUN TESTS ================#
#===========================================#
penalty  <- "enet"
alpha    <- 0.5
lambda   <- exp(seq(-1, -6, length.out = 100)) # is it best to keep lambda fixed across both packages?

nreps <- 5
ntimesteps <- 10

idx.max <- floor(seq(1, nrow(Y), length.out = ntimesteps + 1))[-1]
idx.list <- lapply(idx.max, function(mx) (1:mx))

pt <- proc.time()
time.sim <- replicate(nreps, {
  TIMES <- matrix(NA, nrow = length(idx.list), ncol = 3)
  colnames(TIMES) <- c("biglasso", "biglasso.p", "glmnet")
  
  for (i in 1:length(idx.list)) {
    #print(paste0("Iteration ", i, "..."))
    time.bl <- system.time(
      fit.bl <- biglasso::biglasso(X       = Xbig,
                                   y       = Y,
                                   row.idx = idx.list[[i]],
                                   family  = "cox",
                                   penalty = penalty,
                                   alpha   = alpha,
                                   lambda  = lambda)
    )
    time.bl.p <- system.time(
      fit.bl.p <- biglasso::biglasso(X       = Xbig,
                                     y       = Y,
                                     row.idx = idx.list[[i]],
                                     family  = "cox",
                                     penalty = penalty,
                                     alpha   = alpha,
                                     lambda  = lambda,
                                     ncores  = parallel::detectCores())
    )
    time.gn <- system.time(
      fit.gn <- glmnet::glmnet(x      = X[idx.list[[i]],],
                               y      = Y[idx.list[[i]],],
                               family = "cox",
                               alpha  = alpha,
                               lambda = lambda)
    )
    TIMES[i,] <- unname(c(time.bl[3], time.bl.p[3], time.gn[3]))
  }
  return (TIMES)
})
tot.time <- proc.time() - pt

#============================================#
#================ MAKE PLOTS ================#
#============================================#
time.m <- apply(time.sim, c(1,2), mean)
time.se <- apply(time.sim, c(1, 2), sd)/sqrt(nreps)
time.m.hi <- time.m + time.se
time.m.lo <- time.m - time.se

mycol <- pals::cols25(3) #c("blue", "darkgreen", "red")
mylty <- c("longdash", "dotdash", "dotted")

plot(NA, 
     xlim = range(idx.max), 
     ylim = range(0, time.sim)*1.2,
     xlab = "Rows", 
     ylab = "Time (s)",
     main = "Penalized Cox Regression Runtimes\n(100 tuning parameters)")
grid(); abline(h = 0, v = 0, lwd = 2, col = 'gray75')
for (j in 1:ncol(time.m)) {
  lines(time.m[,j] ~ idx.max, type = "l", pch = 21, bg = "white", lwd = 3,
        col = mycol[j], lty = mylty[j])
}
for (i in 1:dim(time.sim)[3]) {
  for (j in 1:dim(time.sim)[2]) {
    points(time.sim[,j,i] ~ idx.max, pch = 21, cex = 1, bg = mycol[j], lwd = 1.5)
  }
}
legend("topleft", lwd = 2, seg.len = 2.5, 
       col = mycol, lty = mylty, #pt.bg = mycol, pch = 21,
       legend = c("biglasso", 
                  paste0("biglasso (ncores = ", parallel::detectCores(), ")"), 
                  "glmnet"))

# ratio/relative performance to biglasso
time.m.r <- apply(time.m, 2, function(x) x/time.m[,"biglasso"])
plot(NA, 
     xlim = range(idx.max), 
     ylim = c(min(time.m.r)/1.1, max(time.m.r)*1.1),
     xlab = "Rows", 
     ylab = "Relative Runtime",
     main = "Penalized Cox Regression Runtimes\n(100 tuning paramters)")
grid(); abline(h = 0, v = 0, lwd = 2, col = 'gray75')
for (j in 1:ncol(time.m)) {
  lines(time.m.r[,j] ~ idx.max, type = "o", pch = 21, bg = "white", lwd = 2,
        col = mycol[j], lty = mylty[j])
}
legend("topleft", lwd = 2, seg.len = 2.5, 
       col = mycol, lty = mylty, #pt.bg = mycol, pch = 21,
       legend = c("biglasso", 
                  paste0("biglasso (ncores = ", parallel::detectCores(), ")"), 
                  "glmnet"))

tot.time



