spT.validation <- function (z = NULL, zhat = NULL, zhat.Ens = NULL, names = FALSE) 
{
  VMSE <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sum(u^2)/length(u), 4)
  }
  RMSE <<- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sqrt(sum(u^2)/length(u)), 4)
  }
  MAE <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- abs(c(zhat - z))
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
  }
  MAPE <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- abs(c(zhat - z))/abs(z)
    u <- x[!is.na(x)]
    u <- u[!is.infinite(u)]
    round(sum(u)/length(u) * 100, 4)
  }
  NMGE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% 
      as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    NMGE <- abs(c(x[[2]] - x[[1]]))/sum(x[[1]])
    round(mean(NMGE), 4)
  }
  BIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
  }
  rBIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/(length(u) * mean(z, na.rm = TRUE)), 4)
  }
  nBIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    round(sum(x, na.rm = T) * 100/sum(z, na.rm = T), 4)
  }
  rMSEP <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    y <- c(mean(zhat, na.rm = TRUE) - z)
    v <- y[!is.na(y)]
    round(sum(u^2)/sum(v^2), 4)
  }
  Coef <<- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% 
      as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    Coef <- suppressWarnings(cor(x[[2]], x[[1]]))
    round(Coef, 4)
  }
  COE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% 
      as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    COE <- 1 - sum(abs(x[[2]] - x[[1]]))/sum(abs(x[[1]] - 
                                                   mean(x[[1]])))
    round(COE, 4)
  }
  FAC2 <<- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% 
      as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    ratio <- x[[2]]/x[[1]]
    len <- length(ratio)
    if (len > 0) {
      FAC2 <- length(which(ratio >= 0.5 & ratio <= 2))/len
    }
    else {
      FAC2 <- NA
    }
    round(FAC2, 4)
  }
  MGE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% 
      as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    MGE <- round(mean(abs(x[[2]] - x[[1]])), 4)
  }
  IOA <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% 
      as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    LHS <- sum(abs(x[[2]] - x[[1]]))
    RHS <- 2 * sum(abs(x[[1]] - mean(x[[1]])))
    if (LHS <= RHS) 
      IOA <- 1 - LHS/RHS
    else IOA <- RHS/LHS - 1
    round(IOA, 4)
  }
  CRPS <<- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    Nt <- nrow(z)
    n <- ncol(z)
    CRPS <- vector()
    for (p in 1:n) {
      test.index <- which(is.na(z[, p])) %>% as.vector()
      if (length(test.index) == 0) {
        z0 <- z[, p] %>% as.vector()
        zt <- zhat[, p] %>% as.vector()
        s <- sd(zt)
      }
      else {
        z0 <- z[-test.index, p] %>% as.vector()
        zt <- zhat[-test.index, p] %>% as.vector()
        s <- sd(zt)
      }
      CRPS[p] <- verification::crps(z0, cbind(zt, s))$CRPS
    }
    CRPS <- round(mean(CRPS, na.rm = T), 4)
  }
  CRPS.ES <- function(z, zhat, zhat.Ens) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    Nt <- nrow(z)
    n <- ncol(z)
    if (!is.null(zhat.Ens)) {
      test.index <- which(is.na(z), arr.ind = T)
      if (nrow(test.index) > 0) {
        z[test.index] <- 0
      }
      CRPS <- ES <- matrix(NA, nrow = nrow(z), ncol = n)
      for (t in 1:Nt) {
        CRPS[t, ] <- crps_sample(z[t, ], zhat.Ens[t, 
                                                  , ])
        ES[t, ] <- es_sample(z[t, ], zhat.Ens[t, , ])
      }
      CRPS[test.index] <- ES[test.index] <- NA
      CRPS <- round(mean(CRPS, na.rm = T), 4)
      ES <- round(mean(ES, na.rm = T), 4)
    }
    else {
      library(scoringRules)
      CRPS <- ES <- vector()
      for (t in 1:Nt) {
        test.index <- which(is.na(z[t, ])) %>% as.vector()
        if (length(test.index) == 0) {
          z0 <- z[t, ] %>% as.vector()
          zt <- zhat[t, ] %>% as.matrix()
        }
        else {
          z0 <- z[t, -test.index] %>% as.vector()
          zt <- zhat[t, -test.index] %>% as.matrix()
        }
        CRPS[t] <- mean(crps_sample(z0, zt))
        ES[t] <- es_sample(z0, zt)
      }
      CRPS <- round(mean(CRPS, na.rm = T), 4)
      ES <- round(mean(ES, na.rm = T), 4)
    }
    return(list(CRPS = CRPS, ES = ES))
  }
  if (names == TRUE) {
    cat("##\n Mean Squared Error (MSE) \n Root Mean Squared Error (RMSE) \n Mean Absolute Error (MAE) \n Mean Absolute Percentage Error (MAPE) \n Normalized Mean Error(NME: %) \n Bias (BIAS) \n Relative Bias (rBIAS) \n Normalized Mean Bias(NMB: %) \n Relative Mean Separation (rMSEP) \n Correlation coefficient (Coef.) \n The fraction of predictions within a factor of two of observations (FAC2) \n Continuous ranked probability score (CRPS) based on sample mean and standard deviation \n CRPS based on empirical distribution function \n Energy score (ES) based on empirical distribution function \n##\n")
  }
  if ((!is.null(z)) & (!is.null(zhat))) {
    out <- NULL
    out$MSE <- VMSE((z), (zhat))
    out$RMSE <- RMSE((z), (zhat))
    out$MAE <- MAE((z), (zhat))
    out$MAPE <- MAPE((z), (zhat))
    out$BIAS <- BIAS((z), (zhat))
    out$rBIAS <- rBIAS((z), (zhat))
    out$Coef <- Coef((z), (zhat))
    out$FAC2 <- FAC2(z, zhat)
    if (!is.null(zhat.Ens)) {
      R <- CRPS.ES(z, zhat, zhat.Ens)
    }
    else {
      R <- CRPS.ES(z, zhat, zhat.Ens = NULL)
    }
    out$CRPS <- CRPS(z, zhat)
    out$CRPS.sample <- R$CRPS
    out$ES.sample <- R$ES
    out$NMGE <- NMGE((z), (zhat))
    out$MGE <- MGE(z, zhat)
    out$IOA <- IOA(z, zhat)
    out$NMB <- nBIAS((z), (zhat))
    out$rMSEP <- rMSEP((z), (zhat))
    out$COE <- COE(z, zhat)
    unlist(out)
  }
  else {
    return(0)
  }
}