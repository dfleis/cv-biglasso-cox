##
# Testing performance (compute time) of various libraries and their components 
# involving Cox regression and its cross-validation.
##

library(parallel)
library(biglasso)
library(glmnet)
 
library(microbenchmark)

set.seed(124)

pcensor <- 0.2
H0inv <- function(t, scale) t * 1/scale # inverse cumulative hazard fxn for exp times

penalty <- "enet" # only used in biglasso, glmnet will assume lasso if alpha = 1 and ridge if alpha = 0
alpha <- 0.5
lambda <- exp(seq(-1, -6, length.out = 100))

Ns <- c(100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000)
Ps <- c(25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000)

times.bl <- times.gn <- matrix(NA, nrow = length(Ns), ncol = length(Ps)) 

pt0 <- proc.time()
for (i in 1:length(Ns)) {
  n <- Ns[i]
  #print(paste0("i = ", i, ", n = ", n, collapse = ""))
  for (j in 1:length(Ps)) {
    p <- Ps[j]
    
    print(paste0("(n, p) = (", n, ", ", p, ")", collapse = ""))
    
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    Xbig <- as.big.matrix(X)
    beta <- rnorm(p)
    
    U <- runif(n)
    time <- H0inv(-log(U) * exp(-X %*% beta), scale = 1)
    status <- rbinom(n, 1, 1 - pcensor)
    y <- cbind(time, status)
    colnames(y) <- c("time", "status")
    
    pt <- proc.time()
    mod.gn <- glmnet(x = X, y = y, family = "cox", alpha = alpha, lambda = lambda)
    t.gn <- unname(proc.time() - pt)[3]
    
    pt <- proc.time()
    mod.bl <- biglasso(X = Xbig, y = y, family = "cox", penalty = penalty, alpha = alpha, lambda = lambda)
    t.bl <- unname(proc.time() - pt)[3]
    
    times.gn[i,j] <- t.gn
    times.bl[i,j] <- t.bl
  }
}
proc.time() - pt0

R <- 1 - times.gn/times.bl

# mycols <- colorRampPalette(c("red", "blue"))(length(Ps))
# plot(NA, xlim = range(0, Ns), ylim = range(times.bl, times.gn),
#      ylab = "time", xlab = "n")
# grid(); abline(h = 0, v = 0, lwd = 1.5, col = 'gray60')
# for (i in 1:length(Ps)) {
#   lines(times.gn[,i]~Ns, type = 'o', pch = 15, col = mycols[i])
#   lines(times.bl[,i]~Ns, type = 'o', pch = 19, col = mycols[i])
# }
# legend("topleft", legend = c("glmnet", "biglasso"), pch = c(15, 19), lty = 'solid')
# 
# mycols <- colorRampPalette(c("red", "blue"))(length(Ns))
# plot(NA, xlim = range(Ps), ylim = range(times.bl, times.gn),
#      ylab = "time", xlab = "p")
# grid(); abline(h = 0, v = 0, lwd = 1.5, col = 'gray60')
# for (i in 1:length(Ns)) {
#   lines(times.gn[i,]~Ps, type = 'o', pch = 15, col = mycols[i])
#   lines(times.bl[i,]~Ps, type = 'o', pch = 19, col = mycols[i])
# }
# legend("topleft", legend = c("glmnet", "biglasso"), pch = c(15, 19), lty = 'solid')

fields::image.plot(x = Ns, y = Ps, z = R, col = pals::parula(16), log = 'xy')#, zlim = range(0, R))
