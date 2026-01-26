
f1 <- function(t){
  z <- 0.3 * t^11 * (10 * (1 - t))^6 + 10 * (10 * t)^3 * (1 - t)^10 
  return(z/2)
}

siMuIncF <- function(locs,
                     delta = 0.1,
                     para = list(Phis = 0.3, 
                                 nu = 1, 
                                 sigma.sq.s = 1,
                                 sigma.sq.t = 1, 
                                 # Phit = 0.8,
                                 # rho = 0.1,
                                 tau.sq = 0.5, 
                                 beta = c(1, 1)),
                     nRatio = 0.8, iter = 1, W_ts,
                     Uts = NULL,
                     ini.X = NULL){
  
  n       <- nrow(unique(locs[, 1:2]))
  Nt      <- length((unique(locs[, 3])))
  n.train <- n;
  n.test  <- ceiling(0.2*n)
  n       <- n.train #+ n.test
  
  
  if(is.null(Uts)){
    Uts <- matrix(seq(0, 1, len = n*Nt), nrow = Nt, ncol = n)#
  }
  
  Coords <- unique(locs[, 1:2])
  
  D <- fields::rdist(Coords, Coords)
  range(D)
  
  time <- seq(0, 1,, Nt)
  mu <- matrix(NA, nrow = Nt, ncol = n)
  Px <- 2
  X_ts <-  array(0, dim = c(Px, Nt, nrow(Coords)),
                 dimnames = list(paste0("X", 1:Px), as.character(1:Nt),
                                 c(1:n)              
                 ))
  
  X_ts[1,,] <- runif(n*Nt, 0, 1)
  X_ts[2,,] <- rnorm(n*Nt, 0, 1) + cos(10*pi*Uts)
  
  X_ts[1,,] <- X_ts[1,,] - mean(X_ts[1,,] )
  X_ts[2,,] <- X_ts[2,,] - mean(X_ts[2,,] )
  
  Pz <- 1
  Z_ts <- array(1, dim = c(Pz, Nt, nrow(Coords)),
                dimnames = list(paste0("Z", 1:Pz), as.character(1:Nt, c(1:n))
                                
                ))
  thetaF <- array(0, dim = c(Pz, Nt, nrow(Coords)),
                  dimnames = list(paste0("Z", 1:Pz), as.character(1:Nt), 
                                  c(1:n)
                                  
                  ))
  for(t in 1:Nt){
    X_ts[1, t, ] <-  X_ts[1, t, ]
    if(Px > 1){
   
      X_ts[Px, t, ] <- X_ts[Px, t, ] 

      
      mu[t, ] <- as.vector(t(X_ts[1:Px, t, ]) %*% para$beta)
    }else{
      mu[t, ] <-  X_ts[1, t, ]*para$beta[1]
    }
    
    Z_ts[1, t, ] <- rep(1, n)
    thetaF[1, t, ] <- Z_ts[1, t, ]*(f1(Uts[t, ])) - mean(f1(Uts))
    
  }
  
  
  y <- thetaF[1,,] + mu + W_ts  
  
  
  rownames(y) <- 1:Nt
  colnames(y) <- 1:n
  
  index.test <- sample(1:n, n.test, replace = F)
  data <- list(y_ts = y[, ],
               x_ts = X_ts[,,],
               z_ts = Z_ts[1,,],
               w_ts = W_ts[, ],
               u_ts = Uts,
               thetaF = thetaF[1,,],
               Xbeta = mu[,],
               fix.residuals = y[,] - mu[,])
  
  
  gc()
  return(data)
}
