Profile_GEE <- function(y_ts = NULL, 
                        x_ts = NULL, 
                        z_ts = NULL,
                        u_ts = NULL,
                        locs = NULL,
                        new.locs = NULL,
                        intercept = F,
                        method = c("LLR"),
                        Kernel = 3,
                        Neighbor = 10,
                        h        = NULL,
                        h.range  = c(1e-2, 1e0),
                        deta.h   = 0.005,
                        nu = 2,
                        nuUnifb = 1,
                        nThreads = 5,
                        nIter    = 5,
                        # Iter     = 1,
                        tab = NULL,
                        Cts = NULL,
                        Qts = NULL
){
  
  loc      <- unique(locs)
  
  n   <- ncol(y_ts)
  Nt  <- nrow(y_ts)
  y   <- as.vector(y_ts)
  if("matrix" %in% class(x_ts)){
    temp <- x_ts
    x_ts <-  array(0, dim = c(1, nrow(x_ts), ncol(x_ts)),
                   dimnames = list(c("x_ts"),
                                   c(1:nrow(x_ts)),
                                   as.character(1:ncol(x_ts))
                   ))
    x_ts[1,,] <- temp
  }
  
  if("matrix" %in% class(z_ts)){
    Z    <- z_ts
    z_ts <-  array(0, dim = c(1, nrow(z_ts), ncol(z_ts)),
                   dimnames = list(c("z_ts"),
                                   c(1:nrow(z_ts)),
                                   as.character(1:ncol(z_ts))
                   ))
    z_ts[1,,] <- Z
  }
  Px <- dim(x_ts)[1]
  Pz <- dim(z_ts)[1]
  lm.da <- data.frame(y = y)
  if(!is.null(Px)){
    for (i in 1:Px) {
      lm.da <- cbind(lm.da, as.vector(x_ts[i,,]))
    }
    colnames(lm.da) <- c("y", paste0("X", 1:Px))
    fmla <- as.formula(paste("y ~ -1 + ", paste(paste0("X", 1:Px),
                                                collapse = "+")))
    fit <- lm(fmla, data = lm.da)
    curr.beta <- beta <- fit$coefficients
    # curr.beta <- beta
    print(beta)
    
    
  }else{
    curr.beta <- beta <- 0
  }
  # 
  # beta 
  if(!is.null(Px)){
    Xbeta <- X_ts_beta(Nt, x_ts, beta)
    fix.residuals <- y_ts - Xbeta
  }else{
    fix.residuals <- y_ts
  }
  # g(t)
  if(is.null(h)){
    H.provide <- F
    h <- vector()
    h[1] <- mean(h.range)
  }else{
    H.provide <- T 
  }
  # if((method %in% c("WI"))){
  #   # bw.1 <- bandwidth.u(Y_ts = fix.residuals,
  #   #                      U_ts = U_ts,
  #   #                      h.range = h.range,
  #   #                      bwmethod = "cv.aic")
  #   h[1] <- h.range[1]
  # }else{
  h[1] <- h.range[1]
  # }
  
  
  # cat("**************************************\n")
  # LLR.gu <- LLR.u( y_ts     = fix.residuals, 
  #                  z_ts     = z_ts[1,,], 
  #                  u_ts     = u_ts, 
  #                  Kernel   = Kernel,
  #                  h        = h[1])
  LLR.gu <- LLR_u_C_optimized(fix.residuals, 
                          z_ts[1,,], 
                          u_ts, 
                          rep(1, n*Nt),
                          Kernel,
                          h[1])
  
  
  beta               <- beta_est_ts(y_ts, x_ts, LLR.gu, Q = NULL, intercept = F)
  fix.residuals      <- beta$fix.residuals   #by t and s
  
  
  iter  <- 0; 
  error <- 1
  ord <- order(u_ts)
  cat("**************************************\n")
  cat(" h[1]: ", h[1], "\n")
  cat("**************************************\n")
  df_test <- NULL
  if(!is.null(new.locs)){
    new.locs <- unique(new.locs)
    df_test <- as.data.frame(new.locs)
    colnames(df_test)[1:3] <- c("s1", "s2", "t")
    df_test$z <- NA
  }
  fix.semi.residuals.old <- 1
  
  
 
  
  while((iter < nIter)&(error >= 1e-3) ){ 
    # load("./data/non_stationary_covariance.RData")
    # if((Iter == 1) & (iter < 2)){
    
    if(method %nin% c("WI")){
      df_train <- as.data.frame(locs)
      df_train$z <- as.vector(fix.residuals)
      colnames(df_train)[1:3] <- c("s1", "s2", "t")
  
      # 
      
      
      layers_spat <- c(RBF_block(res = 1L),
                       RBF_block(res = 1L))
      layers_temp <- c(AWU(r = 50L, dim = 1L, grad = 200, lims = c(-1, 1)))
      
      tens.Corr <- deepspat_main(f = z ~ s1 + s2 + t - 1, # formula
                                 data       = df_train, #data.frame with colnames s1, s2, t
                                 newdata    = df_test,
                                 g          = ~ 1,
                                 family     = "gneiting",#c("matern","exp","gneiting")
                                 stationary = FALSE,
                                 separable  = FALSE,
                                 nu     = 0.5,
                                 nu_t   = 0.5,
                                 m      = 20L, # numbers of neighbor
                                 nsteps = 100L,
                                 range_init  = 0.5,
                                 method      = "REML",
                                 layers_spat = layers_spat,
                                 layers_temp = layers_temp,
                                 learn_rates = init_learn_rates(eta_mean = 0.01)
      )
      
      
      # Qt.0 <- Qts[1:Nt, 1:Nt]
      # Ct.0 <- Cts[1:Nt, 1:Nt]
      
      
      Qts  <- tens.Corr$prec_obs 
      Cts  <- tens.Corr$cov_obs 
      
      Qt.0 <- tens.Corr$prec_obs[1:Nt, 1:Nt]
      Ct.0 <- tens.Corr$cov_obs[1:Nt, 1:Nt]
      Fit.Wts  <- tens.Corr$pred$df_pred[-c(1:nrow(df_test)),]
      Pred.Wts <- tens.Corr$pred$df_pred[1:nrow(df_test),]
      
      
      
      # k$clear_session()
      # tf$compat$v1$reset_default_graph()
      # gpus <- tf$config$experimental$list_physical_devices('GPU')
      # for (gpu in gpus) {
      #   tf$config$experimental$set_memory_growth(gpu, TRUE)
      # }
      # gc()
      # py_gc <- reticulate::import("gc")
      # py_gc$collect()
    }else{
      Qts      <- Cts <- diag(n*Nt)
      Pred.Wts <- 0
    }
    
    
    
    
    
    # image.plot(Cts1)
    # image.plot(Cts)
    # save(Qts, Cts, file = paste0("./simCovariance_", tab, "_", n, ".RData"))
    # }else{
    #   load(paste0("./simCovariance_", tab, "_", n, ".RData"))
    # }
    # Qts <- Q
    # if(!H.provide){
    # h.candiate <- vector()
    # o <- 0
    #  for(r in seq(0.1, 0.5, length = 10)){
    #   o <- o + 1
    #    bw.1 <- choose_width(Y_st = (fix.residuals),
    #                         time = time,
    #                         loc = loc,
    #                         num = r,
    #                         h.range = h.range,
    #                         bwmethod = "cv.aic")
    #    h.candiate[o] <- bw.1$bw.all$bw
    #  }
    # h[2] <- mean(h.candiate)
    
    
    # h[2] <- bandwidth.u(y_ts     = beta$fix.residuals -
    #                               matrix(Fit.Wts$pred_mean, nrow = Nt,
    #                               ncol = n),
    #                   u_ts     = u_ts,
    #                   h.range = c(3e-2, h.range[2]),
    #                   length.out = 10,
    #                   ckertype  = "gaussian",
    #                   bwmethod = "cv.aic",
    #                   ckerorder = 2,
    #                   bwtype = "adaptive_nn")
    
    # cat("**************************************\n")
    # cat(
    #   " h[2]: ",
    #   h[2],
    #   "\n"
    # )
    # cat("**************************************\n")
    # 
    # 1. Update covariance matrix
    # if((method %in% c("WI"))){
    #   method0 <- "WI"
    #   t2 <- proc.time()
    #   Qt.0 <- diag(Nt)
    # }
    
    
    # if(method %in% c("LCt")){
    #   method0 <- method#c("LCst")
    #   Q0      <- Q[1:Nt, 1:Nt]
    # }
    
    
    
    # if(method %in% c("LCst", "WLLR")){
    #   method0 <- "WLLR"
    #   Q0 <- Q#solve(Cts + diag(nrow(Q))*1e-1) #Q 
    # }
    candidate.h <- seq(h.range[1], h.range[2], by = deta.h)
    CV.mse <- NULL
    for(cv.h in c(candidate.h)){
      K         <- 5
      fold_id.s <- sample(rep(1:K, length.out = n))
      fold_id.t <- sample(rep(1:K, length.out = Nt))
      temp.mse  <- NULL
      for(fold in 1:K){
        train_idx.s <- which(fold_id.s != fold)
        test_idx.s  <- which(fold_id.s == fold)
        
        train_idx.t <- which(fold_id.t != fold)
        test_idx.t  <- which(fold_id.t == fold)
        
        ind.s <- ind.t <- NULL
        for(j in unique(train_idx.s)){
          ind.s <- c(ind.s, ((j - 1)*Nt + 1):((j - 1)*Nt + Nt))
        }
        
       
        for (tt in unique(train_idx.t)) {
          ind.t <- c(ind.t, seq(tt, by = Nt, length.out = n))
        }
        
        Qts.train.s <- Qts[ind.s, ind.s]
        Qts.train.t <- Qts[ind.t, ind.t]
        
        
        if(method %in% c("WI", "LLR")){
          # cat(paste0("Estimating the covariate u's function using LLR", "\n"))
        
          Pred.semi.effects.s <- Pred_LLR_u_C(y_ts = beta$fix.residuals[, train_idx.s],
                                          z_ts = z_ts[1, , train_idx.s],
                                          u_ts = u_ts[, train_idx.s],
                                          Q_ts = rep(1, length(train_idx.s)*Nt),
                                          new_z_ts = as.vector(z_ts[1, , test_idx.s]), 
                                          new_u_ts = as.vector(u_ts[, test_idx.s]), 
                                          Kernel = Kernel,
                                          h = cv.h)
          
          Pred.semi.effects.t <- Pred_LLR_u_C(y_ts = beta$fix.residuals[train_idx.t, ],
                                          z_ts = z_ts[1, train_idx.t, ],
                                          u_ts = u_ts[train_idx.t, ],
                                          Q_ts = rep(1, length(train_idx.t)*n),
                                          new_z_ts = as.vector(z_ts[1, test_idx.t, ]), 
                                          new_u_ts = as.vector(u_ts[test_idx.t, ]), 
                                          Kernel = Kernel,
                                          h = cv.h)
        }
        
        if(method %in% c("LLRC")){
          # cat(paste0("Estimating the covariate u's function using LLR", "\n"))
          
          Pred.semi.effects.s <- Pred_LLR_u_C(y_ts = beta$fix.residuals[, train_idx.s],
                                            z_ts = z_ts[1, , train_idx.s],
                                            u_ts = u_ts[, train_idx.s],
                                            Q_ts = diag(Qts.train.s),
                                            new_z_ts = as.vector(z_ts[1, , test_idx.s]), 
                                            new_u_ts = as.vector(u_ts[, test_idx.s]), 
                                            Kernel = Kernel,
                                            h = cv.h)
          
          Pred.semi.effects.t <- Pred_LLR_u_C(y_ts = beta$fix.residuals[train_idx.t, ],
                                            z_ts = z_ts[1, train_idx.t, ],
                                            u_ts = u_ts[train_idx.t, ],
                                            Q_ts = diag(Qts.train.t),
                                            new_z_ts = as.vector(z_ts[1, test_idx.t, ]), 
                                            new_u_ts = as.vector(u_ts[test_idx.t, ]), 
                                            Kernel = Kernel,
                                            h = cv.h)
        }
        
        if(method %in% c("WLLRC")){
          ini.LR <- LLR_u_C_optimized( y_ts = beta$fix.residuals[, train_idx.s],
                                       z_ts = z_ts[1, , train_idx.s],
                                       u_ts = u_ts[, train_idx.s], 
                                       Q_ts = rep(1, Nt*length(train_idx.s)),
                                       Kernel   = Kernel,
                                       h        = cv.h)
          
          Qt <- diag(diag(Qts.train.s))
          Pred.semi.effects.s <- Pred.KNN.WLLR.u(y_ts = beta$fix.residuals[, train_idx.s],
                                               z_ts = z_ts[1, , train_idx.s],
                                               u_ts = u_ts[, train_idx.s],
                                               new.z_ts = as.vector(z_ts[1, , test_idx.s]), 
                                               new.u_ts = as.vector(u_ts[, test_idx.s]),  
                                               Q        = Qt,
                                               S        = ini.LR$St,
                                               Neighbor = length(train_idx.s)*Nt,
                                               Kernel   = Kernel,
                                               h        = cv.h,
                                               nThreads = nThreads)
          
          
          ini.LR <- LLR_u_C_optimized( y_ts = beta$fix.residuals[train_idx.t, ],
                                       z_ts = z_ts[1, train_idx.t, ],
                                       u_ts = u_ts[train_idx.t, ],
                                       Q_ts = rep(1, n*length(train_idx.t)),
                                       Kernel   = Kernel,
                                       h        = cv.h)
          
          Qt <- diag(diag(Qts.train.t))
          Pred.semi.effects.t <- Pred.KNN.WLLR.u(y_ts = beta$fix.residuals[train_idx.t, ],
                                                 z_ts = z_ts[1, train_idx.t, ],
                                                 u_ts = u_ts[train_idx.t, ],
                                                 new.z_ts = as.vector(z_ts[1, test_idx.t, ]), 
                                                 new.u_ts = as.vector(u_ts[test_idx.t, ]),  
                                                 Q        = Qt,
                                                 S        = ini.LR$St,
                                                 Neighbor = length(train_idx.t)*n,
                                                 Kernel   = Kernel,
                                                 h        = cv.h,
                                                 nThreads = nThreads)
          
        }
        
        if(method %in% c("WLLR")){
          # cat(paste0("Estimating the covariate u's function using ", ceiling(Neighbor), "-NN-WLLR", "\n"))
          ini.LR <- LLR_u_C_optimized( y_ts = beta$fix.residuals[, train_idx.s],
                                       z_ts = z_ts[1, , train_idx.s],
                                       u_ts = u_ts[, train_idx.s],
                                       Q_ts = rep(1, Nt*length(train_idx.s)),
                                       Kernel   = Kernel,
                                       h        = h[1])
          Pred.semi.effects.s <- Pred.KNN.WLLR.u(y_ts = beta$fix.residuals[, train_idx.s],
                                               z_ts = z_ts[1, , train_idx.s],
                                               u_ts = u_ts[, train_idx.s],
                                               new.z_ts = as.vector(z_ts[1, , test_idx.s]), 
                                               new.u_ts = as.vector(u_ts[, test_idx.s]),  
                                               Q        = Qts.train.s,
                                               S        = ini.LR$St,
                                               Neighbor = Neighbor,
                                               Kernel   = Kernel,
                                               h        = cv.h,
                                               nThreads = nThreads)
          
          ini.LR <- LLR_u_C_optimized( y_ts = beta$fix.residuals[train_idx.t, ],
                                       z_ts = z_ts[1, train_idx.t, ],
                                       u_ts = u_ts[train_idx.t, ],
                                       Q_ts = rep(1, n*length(train_idx.t)),
                                       Kernel   = Kernel,
                                       h        = h[1])
          Pred.semi.effects.t <- Pred.KNN.WLLR.u(y_ts = beta$fix.residuals[train_idx.t, ],
                                                 z_ts = z_ts[1, train_idx.t, ],
                                                 u_ts = u_ts[train_idx.t, ],
                                                 new.z_ts = as.vector(z_ts[1, test_idx.t, ]), 
                                                 new.u_ts = as.vector(u_ts[test_idx.t, ]), 
                                                 Q        = Qts.train.t,
                                                 S        = ini.LR$St,
                                                 Neighbor = Neighbor,
                                                 Kernel   = Kernel,
                                                 h        = cv.h,
                                                 nThreads = nThreads)
          
        }
        resi <- c(as.vector(beta$fix.residuals[, test_idx.s]) - as.vector(Pred.semi.effects.s$pred.y),
                  as.vector(beta$fix.residuals[test_idx.t, ]) - as.vector(Pred.semi.effects.t$pred.y))
        
        temp.mse <-  c(temp.mse, mean((resi)^2))
      }
      
      CV.mse <-  c(CV.mse, mean(temp.mse))
      
    }
    # plot(candidate.h, CV.mse, xlab = "h", ylab = "RMSE")
    h[2] <- candidate.h[which.min(CV.mse)]
    cat("**************************************\n")
    cat(
      " h[2]: ",
      h[2],
      "\n"
    )
    cat("**************************************\n")
    
    
    if(method %in% c("WI", "LLR")){
      cat(paste0("Estimating the covariate u's function using LLR", "\n"))
      gu <- LLR_u_C_optimized( y_ts = beta$fix.residuals, 
                               z_ts = z_ts[1,,], 
                               u_ts = u_ts, 
                               Q_ts = rep(1, n*Nt),
                               Kernel   = Kernel,
                               h        = h[2])
      
      
      # plot(U_ts[ord], gu$alpha[ord, 1], type = "l", col = "gray30")
      
    }
     #LLRC
    if(method %in% c("LLRC")){
      
      gu <- LLR_u_C_optimized(y_ts = beta$fix.residuals, 
                              z_ts = z_ts[1,,], 
                              u_ts = u_ts, 
                              Q_ts = diag(Qts),
                              Kernel   = Kernel,
                              h        = h[2])
      
    }
     
    #3. Update theta 
    ini.LR <- NULL
    if(method %in% c("WLLRC")){
      ini.LR <- LLR_u_C_optimized(y_ts = beta$fix.residuals, 
                                  z_ts = z_ts[1,,], 
                                  u_ts = u_ts, 
                                  Q_ts = rep(1, n*Nt),
                                  Kernel   = Kernel,
                                  h        = h[1]- 0.01)
      # Qt <- Qt.0
      # for(s in 2:n){
      #   Qt <- bdiag(Qt, Qt.0);
      # }
      Qt <- diag(diag(Qts))
      gu <- KNN.WLLR.u(y_ts = beta$fix.residuals,
                       z_ts = z_ts[1,,],
                       u_ts = u_ts,
                       Q = as.matrix(Qt),
                       S = ini.LR$St,
                       Neighbor = n*Nt,
                       Kernel   = Kernel,
                       h        = h[2],#h.candidate[which.min(Cv)],
                       nThreads = nThreads)
      
    }
    if(method %in% c("WLLR")){
      cat(paste0("Estimating the covariate u's function using ", ceiling(Neighbor), "-NN-WLLR", "\n"))
      ini.LR <- LLR_u_C_optimized( y_ts = beta$fix.residuals, 
                       z_ts = z_ts[1,,], 
                       u_ts = u_ts, 
                       Q_ts = rep(1, n*Nt),
                       Kernel   = Kernel,
                       h        = h[1]- 0.01)
      
      
      
      #  Cv <- vector()
      #  h.candidate <- seq(h.range[1], h.range[2], length = 20)
      #   for(r in 1:length(h.candidate)){
      #   gu <- KNN.WLLR.u(y_ts = beta$fix.residuals,
      #                                z_ts = z_ts[1,,],
      #                                u_ts = u_ts,
      #                                Q = Qts,
      #                                S = ini.LR$St,
      #                                Neighbor = ceiling(Neighbor),
      #                                Kernel   = Kernel,
      #                                h        = h.candidate[r],
      #                                nThreads = nThreads)
      #   Sii <- diag(gu$St)
      #   y <- as.vector(beta$fix.residuals)
      #   Y.hat <- gu$St %*% y
      #   fy <- 0
      #   for(i in 1:(n*Nt)){
      #     fi <- (Y.hat[i] - Sii[i]*y[i])/(1 -  Sii[i])
      #     fy <- fy + (y[i] - fi)^2
      #   }
      #   Cv[r] <- fy/(n*Nt)*(1 + sum(Sii)/(n*Nt))
      # }
      #  plot(h.candidate, Cv)
      #  h[2] <- h.candidate[which.min(Cv)]
      gu <- KNN.WLLR.u(y_ts = beta$fix.residuals,
                       z_ts = z_ts[1,,],
                       u_ts = u_ts,
                       Q = Qts,
                       S = ini.LR$St,
                       Neighbor = ceiling(Neighbor),
                       Kernel   = Kernel,
                       h        = h[2],#h.candidate[which.min(Cv)],
                       nThreads = nThreads)
      # plot(u_ts[ord], ini.LR$alpha[ord, 1], type = "l", col = "gray30")
      # lines(u_ts[ord], gu$alpha[ord, 1], col = "red")
      cat(".....\n")
      # rm(ini.LR)
      # gc()
    }
  
    # # under the assumption on independent space
    # if(method %in% c("WCXt")){
    #   method0 <- method#
    #   Q0 <- Qt.0; 
    #   for(s in 2:n){
    #     Q0 <- bdiag(Q0, Qt.0);
    #   }
    #   theta <- semPar.space.time( y_st = t(fix.residuals), # by s and t
    #                               z_st = aperm(z_ts, c(1, 3, 2)),
    #                               Time = time,
    #                               Q = as.matrix(Q0),
    #                               S = theta$St,
    #                               Neighbor = Nt,#ceil(Neighbor/2),
    #                               Kernel = Kernel[1],
    #                               h = h[2],
    #                               nuUnifb = nuUnifb,
    #                               nu = 1,
    #                               nThreads = nThreads,
    #                               method = method0)
    #   # theta$alpha[, 1]
    #   # theta <- theta_Wang_Space_Inde(y_st = fix.residuals,
    #   #                                z_st = z_st,
    #   #                                Time = time,
    #   #                                Q = Q0,
    #   #                                S0 = theta0$St,
    #   #                                Kernel = Kernel[1],
    #   #                                h = h[2])
    #   # theta1$alpha[, 1]
    # }
    # save(gu, file = "./Simulation_Ux/gu.RData")
    # load("./Simulation_Ux/gu.RData")
    #3. Update beta
    if(!is.null(Px)){
      if(method %in% c("WI")){
        # Q <- diag(n*Nt)
        beta    <- beta_est_ts(y_ts, x_ts, gu, Q =Qts, intercept = intercept)
        beta.sd <- beta$W %*% beta$R %*% Qts %*% t(beta$R) %*% t(beta$W)
      }
      
      # if(method %in% c("WCXt")){
      #   Ct <- Ct.0
      #   for(s in 2:n){
      #     Ct <- bdiag(Ct, Ct.0)
      #   }
      #   
      #   beta    <- beta_est_ts(y_ts, x_ts, gu, Q = Qt, intercept = intercept)
      #   beta.sd <- beta$W %*% beta$R %*% Ct %*% t(beta$R) %*% t(beta$W)
      # }
      
      if(method %in% c("LLR", "LLRC", "WLLRC", "WLLR")){
        beta    <- beta_est_ts(y_ts, x_ts, gu, Q = Qts, intercept = intercept)
        beta.sd <- beta$W %*% beta$R %*% Cts %*% t(beta$R) %*% t(beta$W)
      }
      
      
      fix.residuals  <- beta$fix.residuals  
      # plot(as.vector(fix.residuals))
      Coef.           <- as.matrix(beta$beta)
      
    }else{
      fix.effect.fit <- 0
      fix.residuals  <- y_ts 
    }
    # intercept.est      <- mean(beta$fix.residuals   - gu$y.fit)
    fix.semi.residuals <- mean((fix.residuals - gu$fit.y)^2) #- intercept.est
    
    error     <-  abs(fix.semi.residuals.old - fix.semi.residuals)
    error2     <-  abs(fix.semi.residuals.old - fix.semi.residuals)/fix.semi.residuals.old
    fix.semi.residuals.old <- fix.semi.residuals
    #mean(abs(curr.beta - Coef.)/abs(Coef.))
    curr.beta <- Coef.
    iter      <- iter + 1
    cat("--------------------------------------\n")
    cat("**************************************\n")
    cat(
      " iter: ",
      iter,
      "\n",
      "beta: ",
      Round(Coef., 4),
      "\n",
      "beta.sd: ",
      Round(as.vector(diag(beta.sd)), 4),
      "\n",
      # "var: ",
      # Round(Cs$Var.s, 4),
      # "\n",
      "error: ",
      Round(error, 4),
      "\n",
      "error2: ",
      Round(error2, 4),
      "\n"
    )
    cat("**************************************\n")
    cat("--------------------------------------\n")
    if(iter == 1){error <- 1}
  }
  return(list(
    beta               = Coef., 
    g                  = gu, 
    St                 = gu$St,
    Qts                = Qts,
    beta_sd            = sqrt(diag(beta.sd)),
    ini.LR             = ini.LR,
    h                  = h,
    fix.residuals      = fix.residuals,
    randome.effects    = fix.residuals - gu$fit.y,
    fit.value          = fix.residuals + gu$fit.y,
    w_ts                = Pred.Wts,
    nIter              = iter))
} 

