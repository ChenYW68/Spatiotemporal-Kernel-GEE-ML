# Load packages, functions, and environment
library(reticulate)
# use_condaenv("r-reticulate")
# use_condaenv("TFv1-15", required = TRUE)
# conda_list()

stationary <- TRUE
# source("./ML/utils.R")

# Dataset 1. From nonstationary, separable covariance
# Simulate warped location
# set.seed(16)
range_s    <- 1  # c(0.1, 0.5, 1)
range_t    <- 1   # c(0.1, 0.5, 1)
if(range_s == 0.1){
  tab.s <- "slow"
}
if(range_s == 0.5){
  tab.s <- "smid"
}
if(range_s > 0.5){
  tab.s <- "supp"
}

if(range_t == 0.1){
  tab.t <- "tlow"
}
if(range_t == 0.5){
  tab.t <- "tmid"
}
if(range_t > 0.5){
  tab.t <- "tupp"
}

m <- 20             # c(10, 20)
n <- 20              # c(20, 30)

r1 <- 50
layers <- c(AWU(r = r1, dim = 1L, grad = 50, lims = c(-1, 1)),
            AWU(r = r1, dim = 2L, grad = 50, lims = c(-1, 1)),
            RBF_block(res = 1L))
nlayers <- length(layers)
eta <- list()
eta[[1]] <- sin(seq(0, pi, length.out = r1))
eta[[2]] <- c(1, rep(0, r1 - 1))
for(j in 3:(nlayers)) eta[[j]] <- runif(n = 1, min = -1, max = exp(3/2)/2)

# spatial_grid <- expand.grid(seq(-1, 1, length.out = sqrt(n)),
#                  seq(-1, 1, length.out = sqrt(n))) %>% as.matrix()


x.coords <- seq(-1, sqrt(n)/10, , n)#x.0 + (1:(n - 1))*delta
y.coords <- seq(-1, sqrt(n)/10, , n)#

spatial_grid <- matrix(NA, ncol = 1, nrow = 1)
while (nrow(spatial_grid)!=n) {
  spatial_grid <- cbind(sample(x.coords, n, replace = T),
                        sample(y.coords, n, replace = T)
  )
  spatial_grid <- unique(spatial_grid)
}



id.s1 <- as.data.frame(spatial_grid)
id.s1$id <- 1:n
colnames(id.s1) <- c("x", "y", "id")


swarped <- spatial_grid
for(j in 1: (nlayers)) {
  swarped <- layers[[j]]$fR(swarped, eta[[j]]) %>% scal_0_5_mat()
}
id.s <- as.data.frame(swarped)
id.s$id <- 1:n
colnames(id.s) <- c("x", "y", "id")


# time 
layers_temp <- c(AWU(r = 50L, dim = 1L, grad = 200, lims = c(-1, 1)))
nlayers_temp <- length(layers_temp)
eta_t <- list()
for(j in 1:nlayers_temp) {
  eta_t[[j]] <- c(rep(0.2, 15), rep(0.5, 10), rep(0.7, 10), rep(0.2, 15))
}

twarped <- temporal_grid <- as.matrix(seq(-1, sqrt(m)/10, length.out = m))
for(j in 1: (nlayers_temp)){
  twarped <- layers_temp[[j]]$fR(twarped, eta_t[[j]]) %>% scal_0_5_mat()
}

warped <- cbind(matrix(c(rep(swarped[,1], length(twarped)), rep(swarped[,2], length(twarped))), ncol = 2),
                rep(twarped, each = nrow(swarped)))


warped <- as.data.frame(warped)
colnames(warped) <- c("x", "y", "time")
warped <- warped %>% left_join(id.s, by = c("x", "y"))
setorderv(warped, c("id", "time"))



locs <- cbind(matrix(c(rep(spatial_grid[,1], length(temporal_grid)), 
                       rep(spatial_grid[,2], length(temporal_grid))), ncol = 2),
              rep(temporal_grid, each = nrow(spatial_grid)))
locs <- as.data.frame(locs)
colnames(locs) <- c("x", "y", "time")
locs <- locs %>% left_join(id.s1, by = c("x", "y"))
setorderv(locs, c("id", "time"))

# 计算空间距离和时间距离
range(rdist(swarped[, 1:2]))
range(rdist(twarped[, 1]))

spatial_dist <- rdist(warped[,1:2])  # 空间距离
time_dist    <- rdist(warped[,3])    # 时间距离


a          <- 0.1
b          <- 0.1
c          <- 0.5
sigma2     <- 1



# if(stationary){
  spatial_dist <- rdist(locs[,1:2])^2  # 空间距离
  time_dist    <- rdist(locs[,3])^2       # 时间距离
  term_1       <- 1 / (1 + ( time_dist^a / range_t ) )^b
  term_2       <- exp( - (spatial_dist * term_1)^c / range_s )
  Cts.1          <- sigma2 * term_1 * term_2 # covariance matrix
  
  diag(Cts.1) <- diag(Cts.1) + 5E-2
  Qts.1       <- solve(Cts.1)
  
  
  
  spatial_dist <- rdist(warped[,1:2])^2  # 空间距离
  time_dist    <- rdist(warped[,3])^2       # 时间距离
  term_1       <- 1 / (1 + ( time_dist^a / range_t ) )^b
  term_2       <- exp( - (spatial_dist * term_1)^c / range_s)
  Cts          <- sigma2 * term_1 * term_2 # covariance matrix


  
  
# }else{
#   spatial_dist <- rdist(warped[,1:2])^2  # 空间距离
#   time_dist    <- rdist(warped[,3])^2       # 时间距离
#   term_1       <- 1 / (1 + ( time_dist^a / range_t ) )^b
#   term_2       <- exp( - (spatial_dist * term_1)^c / range_s)
#   Cts          <- sigma2 * term_1 * term_2 # covariance matrix
# }

# Cts <- time_cov*spatial_cov + diag(nrow(time_cov))*1E-1
# Cts[1:20, 1:10]
diag(Cts) <- diag(Cts) + 5E-2
Qts       <- solve(Cts)
K <- t(chol(Cts))
save(
  K, 
  Qts, 
  Cts,
  Cts.1,
  Qts.1,
  locs, 
  file = paste0("./data/cov_", 
                m, "_",
                n, "_", 
                tab.s, "_",
                tab.t,
                ".RData"))
