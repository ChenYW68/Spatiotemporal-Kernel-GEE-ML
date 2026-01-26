# .rs.restartR()   # or restart R manually
# gc()             # garbage collect
# rm(list = ls())  # clear workspace
# Sys.setenv(OMP_THREAD_LIMIT=1)
source("./LoadPackages/RDependPackages.R")
source("./R/spT.validation.R")
source("./Simulation_Ux/Ux_sim_Data.R")
source("./Simulation_Ux/Ux_Profile_GEE.R")
# library(tensorflow)
source("./ML/utils.R")
n  <- 10
Nt <- 10

tab    <- c("slow_tlow", "supp_tlow", "slow_tupp", "supp_tupp", "smid_tmid")[5]
load(paste0("./data/cov_", 
            Nt, "_", n, "_", 
            tab, ".RData"))
tab <- paste0("Est_Table3_", tab)
# n      <- nrow(unique(locs[, 1:2]))
# Nt     <- length((unique(locs[, 3])))
set.seed(1234)

ini.X <- 0#matrix(runif(n*Nt, 0, 1), nrow = Nt, ncol = n)
ini.uts <- ini.X + matrix(runif(n*Nt, 0, 1), nrow = Nt, ncol = n)
test.index <- sort(sample(1:n, ceiling(0.2*n), replace = FALSE))
nThreads <- 50
test.index
# options(device = NULL)
######################################################a########
h.range <- c(2E-2, 6e-2)            # minimum and maximum of bandwidth
M      <- c("WI", "LLR", "LLRC", "WLLRC", "WLLR") # five methods for estimating function g(t)
L      <- c(5)               # index of methods
Kernel <- 3                     # 0. exponential kernel; 1. Bisperia Kernel; 3. Gaussian kernel
nIter  <- 3                 # the number of iteration
Neig0  <- ceiling(10.91 + 0.05*n*Nt)    #seq(2, Nt, by = 2)    # K value of K-nearest-weighted local linear regression (K-WLLR) only when K-WLLR is set
n.sim  <- 300           # the number of simulations
# c(4, 6, Nt)#c(6, Nt)#, n*Nt

Neig0

##############################################################
seed <- ceil(seq(3, 1500,, 1000))
##############################################################
yc     <- c(-6, 6) # for making Figures

##############################################################
par(mfrow = c(1, 2))
##############################################################
##############################################################

for(l in c(L)){ #1, 2, 3, 4, 5; 1, 4, 5
  if(l %in% c(1, 2, 3)){
    Neig <- n*Nt#length((unique(locs[, 3])))
  }else{
    Neig <- Neig0
  }
  method <- M[l]
 
  for(Neighbor in Neig)
  {
    # tab0   <- paste0("ST_", Neighbor)
    # ##############################################################
    # tab       <- paste0(method, "_", n, "_", Nt, "_", tab0)
    # alpha.est <- paste0(method, "_A_", n, "_", Nt, "_", tab0)
    ##############################################################
    iter <- n.k <- 17
    Result <- alpha.e <- NULL
    ##############################################################
    while(iter <= n.sim){
      # if(iter%%5 == 0){
      #   .rs.restartR() 
      # }
      set.seed(seed[n.k])
      para = list(
        beta = c(1, 1),
        # sill = c(sigma.sq.s),
        # scale_s = as.numeric(Phis),
        # scale_t = as.numeric(Phit), # the smoothness in time: (0, 1]
        # sep = 0.5, # space–time interaction: [0, 1]   0.5
        # power_s = 3, #(0, 2]   0.1
        # power_t = 0.5,         # >= 0   0.5
        nugget = 5e-2  #0.1
      )
      
      # if(iter == 1){
      # W_ts <- matrix(K %*% rnorm(nrow(K)), nrow = Nt, ncol = n)
      
      W_ts <- matrix(mvnfast::rmvt(1, mu = rep(0, n*Nt),
                                   sigma = Cts, df  = 5,
                                   ncores = 10),
                     nrow = Nt, ncol = n)
      # 
      # W_ts <- t(matrix(K %*% rnorm(nrow(K)), nrow = n, ncol = Nt))
      #W_ts <- (matrix(K %*% rnorm(nrow(K)), nrow = Nt, ncol = n))
      # }
      simDa <- siMuIncF(locs, para = para, iter = iter, W_ts = W_ts - mean(W_ts), Uts = ini.uts, ini.X = ini.X)
      
      ###################################################################
      inx1 <- which(locs$id %nin%test.index)
      data <- list(y_ts = simDa$y_ts[, -test.index],
                   x_ts     = simDa$x_ts[,, -test.index],
                   z_ts     = simDa$z_ts[, -test.index],
                   u_ts     = simDa$u_ts[, -test.index],
                   w_ts     = simDa$w_ts[, -test.index],
                   locs     = locs[locs$id %nin%test.index, ],
                   Cts      =  Cts[inx1, inx1],
                   Qts      =  Qts[inx1, inx1])
      inx2 <- which(locs$id %in%test.index)
      test <- list(y_ts = simDa$y_ts[, test.index],
                   x_ts     = simDa$x_ts[,, test.index],
                   z_ts     = simDa$z_ts[, test.index],
                   u_ts     = simDa$u_ts[, test.index],
                   w_ts     = simDa$w_ts[, test.index],
                   locs     = locs[locs$id %in%test.index, ],
                   Cts      =  Cts[inx2, inx1]) 
      
      
      
      ord  <- order(data$u_ts)
      u_ts <- data$u_ts[ord]
      
      g.true <- f1(u_ts) - mean(f1(u_ts))
      
      # library(profvis)
      start_time <- proc.time()
      
      tensorflow::tf$compat$v1$reset_default_graph()
      sess <- tf$compat$v1$Session()
      on.exit(sess$close(), add = TRUE)
      
      library(tensorflow)
      k <- keras::backend()
      Fit <- Profile_GEE(y_ts     = data$y_ts,
                         x_ts     = data$x_ts,
                         z_ts     = data$z_ts,
                         u_ts     = data$u_ts,
                         locs     = data$locs,
                         new.locs = test$locs,
                         Kernel   = Kernel,
                         method   = method,
                         Neighbor = Neighbor,
                         h.range  = h.range, #c(quantile(diff(u_ts), 1.5e-1),
                         #quantile(diff(u_ts), 5e-1)),
                         nu       = 2,
                         nuUnifb  = 1,
                         nThreads = nThreads,
                         nIter    = nIter,
                         Cts      = data$Cts,
                         Qts      = data$Qts,
                         tab      = tab#,
                         # Iter     = iter
                         ) 
      k$clear_session()
      tf$compat$v1$reset_default_graph()
      gpus <- tf$config$experimental$list_physical_devices('GPU')
      for (gpu in gpus) {
        tf$config$experimental$set_memory_growth(gpu, TRUE)
      }
      gc()
      py_gc <- reticulate::import("gc")
      py_gc$collect()
      
      if(method %in% c("WI", "LLR")){
        Pred.semi.effects <- Pred_LLR_u_C(y_ts = Fit$fix.residuals,
                                          z_ts = data$z_ts,
                                          u_ts = data$u_ts,
                                          Q_ts = rep(1, nrow(Fit$Qts)),
                                          new_z_ts = as.vector(test$z_ts),
                                          new_u_ts = as.vector(test$u_ts),
                                          Kernel = Kernel,
                                          h = Fit$h[2])
      }
      if(method %in% c("LLRC")){
       
        Pred.semi.effects <- Pred_LLR_u_C(y_ts = Fit$fix.residuals,
                                          z_ts = data$z_ts,
                                          u_ts = data$u_ts,
                                          Q_ts = diag(Fit$Qts),
                                          new_z_ts = as.vector(test$z_ts),
                                          new_u_ts = as.vector(test$u_ts),
                                          Kernel = Kernel,
                                          h = Fit$h[2])
      }
      
      if(method == "WLLRC"){
        Qt <- diag(diag(Fit$Qts))

        Pred.semi.effects <- Pred.KNN.WLLR.u(y_ts = Fit$fix.residuals,
                                             z_ts = data$z_ts,
                                             u_ts = data$u_ts,
                                             new.z_ts = as.vector(test$z_ts),
                                             new.u_ts = as.vector(test$u_ts),
                                             Q = as.matrix(Qt),
                                             S =  Fit$ini.LR$St,
                                             Neighbor = n*Nt,
                                             Kernel   = Kernel,
                                             h        = Fit$h[2],#h.candidate[which.min(Cv)],
                                             nThreads = nThreads)
      }
      if(method == "WLLR"){
        Pred.semi.effects <- Pred.KNN.WLLR.u(y_ts = Fit$fix.residuals,
                                             z_ts = data$z_ts,
                                             u_ts = data$u_ts,
                                             new.z_ts = as.vector(test$z_ts),
                                             new.u_ts = as.vector(test$u_ts),
                                             Q        = Fit$Qts,
                                             S        = Fit$ini.LR$St,
                                             Neighbor = Neighbor,
                                             Kernel   = Kernel,
                                             h        = Fit$h[2],
                                             nThreads = nThreads)
      }
     
     
      
      if(method %in% c("WI")){
        Pred.rand.effects <- 0
      }else{
        # Pred.rand.effects <- as.vector(test$w_ts)
        Pred.rand.effects <- as.vector(Fit$w_ts$pred_mean)
      }
      
      Pred.fixed.effects <- as.vector(X_ts_beta(Nt, test$x_ts, Fit$beta))
      Pred.semi.effects <- Pred.semi.effects$pred.y
   
      #Fit$Pred.Wts$pred_mean#test$Cts %*% Fit$Qts %*% as.vector(Fit$randome.effects)
      Pred.y <- Pred.fixed.effects + as.vector(Pred.semi.effects) + Pred.rand.effects#
      true.y <- as.vector(test$y_ts)
     
      y.Pred.error <- as.vector(spT.validation(Pred.y, true.y)[c(2, 3, 9)])
      y.Pred.error
      
      plot(true.y, Pred.y, xlim = c(range(c(true.y, Pred.y))), ylim = c(range(c(true.y, Pred.y))))
      # test$Cts %*% data$Qts %*% as.vector(W_ts[, -test.index])
      
      # loess_fit <- loess(as.vector(y_ts) ~ rep(simDa$time, times = n), 
      #                    span = 0.4, degree = 1)
      # 
      # v_grid <- rep(simDa$time, times = n)
      # g_est <- predict(loess_fit, newdata = data.frame(v = v_grid))
      # g_true_vals <- f1(v_grid) - mean(f1(v_grid))
      # 
      # # ---------------------------
      # # 5. Plot the true function versus the estimated function
      # # ---------------------------
      # plot(v_grid, g_true_vals, type = "l", col = "blue", lwd = 2,
      #      xlab = "v = z * beta", ylab = "g(v)",
      #      main = "True g(v) (blue) vs. Estimated g(v) (red)")
      # lines(v_grid, g_est, col = "red", lwd = 2)
      # legend("topright", legend = c("True g(v)", "Estimated g(v)"),
      #        col = c("blue", "red"), lwd = 2)
      
      
      # plot(Fit$theta$alpha[, 1], type = "l")
      objects_over_10MB <- sapply(ls(), function(x) object.size(get(x))/1024^2)
      objects_over_10MB <- objects_over_10MB[objects_over_10MB > 10]
      print(sort(objects_over_10MB, decreasing = TRUE))
      # 1. nonparametric component
      {
        
        est.g <- Fit$g$alpha[ord, 1] - mean(Fit$g$alpha[ord, 1])
        h       <- Fit$h
        error      <- as.vector(spT.validation(as.vector(g.true), est.g)[c(2, 3, 9)]) 
        error
        # if(!is.null(simDa$theta)){
        alpha <- data.frame(t         = 1:Nt, 
                            u_ts      = u_ts,
                            true.gt   = g.true,
                            esti.gt   = est.g,
                            H_0       = h[1],
                            H_1       = h[2],
                            iter      = iter)
        
        # }
        
        
        
        
        {
          # plot(model.np.1, plot.errors.method = "bootstrap",
          #      plot.errors.boot.num = 25)
          plot(u_ts,
               g.true, 
               col  = "gray50",
               cex  = 1, 
               ylim = c(yc[1], yc[2]), 
               xlab = "u", 
               ylab = "g(u)",
               pch = 19)
          lines(u_ts, est.g, lwd = 3, col = "red")
          
          points(u_ts, Fit$fix.residuals[ord], cex = 0.5, pch = 21)
          
          
          # for(t1 in 1:Nt){
          #   points(rep(simDa$time[t1], nrow(bw.1$s)), 
          #          simDa$fix.residuals[t1, bw.1$s[, t1]], 
          #          cex = 0.3, pch = 19, col = "blue")
          # }
          text(0.2, yc[1] + 0.2, paste0("n = ", n))
          text(0.6, yc[1] + 0.2, paste0("RMSE = ", Round(error[1], 4)))
          text(0.2, yc[2] - 0.2, paste0("h = ", Round(h[2], 4)))
          text(0.6, yc[2] - 0.2, paste0("iter = ", iter))
          text(0.8, yc[2] - 1.5, paste0(Neighbor,"-NN-", method))
          text(0.2, yc[1] + 1, tab)
        }
        # dev.off()
        # if(l %in% L[1]){
        #   g_t <- model.np.1
        #   save(g_t, file = "./data/g_t.RData")
        # }
        # load(file = "./data/g_t.RData")
        # 
        # pdf(file = paste0("./figure/g_t_LLR_WLLR.pdf"),
        #     width = 8, height = 6)
        # par(mar = c(3.5, 3.5, 2, 0.5) + 0, mgp = c(2, 0.8, 0))
        # par(mfrow = c(1, 1),
        #     cex = 1.2,
        #     cex.axis = 1.3,
        #     cex.lab = 1.2,
        #     cex.main = 1,
        #     lwd = 1)
        # {
        #   # plot(model.np.1, plot.errors.method = "bootstrap",
        #   #      plot.errors.boot.num = 25)
        #   plot(simDa$time, f1.true, col = "black", cex = 2, pch = 20,
        #        ylim = c(yc[1], yc[2]), xlab = "t", ylab = "g(t)")
        #   lines(simDa$time, model.np.1, lwd = 3, col = "blue")
        #   lines(simDa$time, g_t, lwd = 3, col = "red")
        #   for(s in 1:n){
        #     points(simDa$time, simDa$fix.residuals[,s], cex = 0.1, pch = 21)
        #   }
        #   for(t1 in 1:Nt){
        #     points(rep(simDa$time[t1], nrow(bw.1$s)),
        #            simDa$fix.residuals[t1, bw.1$s[, t1]],
        #            cex = 0.3, pch = 19, col = "blue")
        #   }
        #   legend("topleft", legend = paste(c("Sampling points", "LLR", "K-NN-WLLR")),
        #          col = c("black", "blue", "red"),
        #          pch = c(20, NA, NA),  # 1 for points, NA for lines
        #          lty = c(NA, 1, 1),    # NA for points, 1 for lines,
        #          lwd = 5)
        # }
        # dev.off()
      }
      
      # fix.semi.residuals <- t(Fit$fix.semi.residuals)
      {
        end_time <- proc.time()
        run_time <- paste0(if_else(month(Sys.Date()) > 9
                                   , as.character(month(Sys.Date()))
                                   , paste0("0",  as.character(month(Sys.Date()))))
                           , "_", if_else(day(Sys.Date()) > 9
                                          , as.character(day(Sys.Date()))
                                          , paste0("0",  as.character(day(Sys.Date()))))
                           , "_", if_else(hour(Sys.time()) > 9 
                                          , as.character(hour(Sys.time()))
                                          , paste0("0",  as.character(hour(Sys.time()))))
                           
        )
        # 2. parametric component
        # plot(as.vector(simDa$thetaF + simDa$Xbeta), as.vector(Fit$fit.value))
        # points(as.vector(Fit$y_ts), as.vector(Fit$fit.value), col = "red")
        # y.fitted.error <- as.vector(spT.validation(as.vector(Fit$fit.value), as.vector(data$y_ts))[c(2, 3, 9)])
      
        temp0 <- data.frame(Iter         = iter,
                            nIter        = Fit$nIter,
                            H0           = Fit$h[1],
                            H1           = Fit$h[2],
                            Beta.1       = Fit$beta[1],
                            Beta.2       = Fit$beta[2],
                            beta.sd1     = Fit$beta_sd[1],
                            beta.sd2     = Fit$beta_sd[2],
                            g.RMSE       = error[1],
                            g.MAE        = error[2],
                            g.CRPS       = error[3],
                            # y.RMSE       = y.fitted.error[1],
                            # y.MAE        = y.fitted.error[2],
                            # y.CRPS       = y.fitted.error[3],
                            y.RMSE       = y.Pred.error[1],
                            y.MAE        = y.Pred.error[2],
                            y.CRPS       = y.Pred.error[3],
                            Neighbor     = Neighbor,
                            elapsed      = round((end_time - start_time)[[3]], 3))
        
        alpha$run_time <- run_time
        temp0$run_time <- run_time
        alpha$Method   <- paste0(M[l], "_", Neighbor)
        temp0$Method   <- paste0(M[l], "_", Neighbor)
        
        
        if(iter > 1){
          load(paste0("./Result/",tab, "_Nt-", 
                      Nt, "_n-", n, "_",
                      method, "-", 
                      Neighbor, ".RData"))
        }
        alpha.e          <- rbind(alpha.e, alpha)
        Result           <- rbind(Result, temp0)
        rownames(Result) <- NULL
        # cat("..................", tab, "..................\n")
        col <- c("H0", 
                 "H1", 
                 "Beta.1", 
                 "Beta.2", 
                 "g.RMSE",
                 "y.RMSE",
                 # "test.y.RMSE", 
                 # "test.y.RMSE",
                 # "test.fix.RMSE", 
                 # "test.g.RMSE",
                 # "test.random.RMSE",
                 "Neighbor",
                 "nIter")
        da <- as.data.frame(Result)[, colnames(Result) %in% col]
        setDT(da)
        print(round((da), 4))
        cat("..................", paste0(M[l], "_", Neighbor), "..................\n")
        print(round(colMeans(da), 4))
        cat(".................. iter = ", iter, "..................\n\n")
        iter <- iter + 1
      }
      n.k <- n.k + 1
      unique(Result$Method)
      save(Result, alpha.e, 
           file = paste0("./Result/",tab, "_Nt-", 
                         Nt, "_n-", n, "_",
                         method, "-", 
                         Neighbor, ".RData"))
    }
  }
  
}





