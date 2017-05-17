# Importation of ourside functions
#' @importFrom stats runif rnorm lm arima.sim var acf rexp spec.taper
#' @importFrom utils tail read.table
NULL
if(F){
  ###############################################################################
  # 3*3 Gaussian mixture simulation
  data1 <- array(dim = c(2,30,3))
  for(ep in 1:3){
    data1[,,ep] <- t(rbind(MASS::mvrnorm(10, mu = c(0,10), Sigma = 5*diag(2)),
                           MASS::mvrnorm(10, mu = c(10,0), Sigma = 5*diag(2)),
                           MASS::mvrnorm(10, mu = c(5,5), Sigma = 5*diag(2))))
  }
  data1[,,2] <- data1[,c(11:30,1:10),2]
  data1[,,3] <- data1[,c(3:30,1,2),3]

  data2 <- array(dim = c(2,30,3))
  for(ep in 1:3){
    data2[,,ep] <- t(rbind(MASS::mvrnorm(10, mu = c(0,5), Sigma = 4*diag(2)),
                           MASS::mvrnorm(10, mu = c(5,0), Sigma = 4*diag(2)),
                           MASS::mvrnorm(10, mu = c(5,5), Sigma = 4*diag(2))))
  }
  data2[,,2] <- data2[,c(11:30,1:10),2]
  data2[,,3] <- data2[,c(3:30,1,2),3]

  data3 <- array(dim = c(2,30,3))

  for(ep in 1:3){
    data3[,,ep] <- t(rbind(MASS::mvrnorm(10, mu = c(0,.20), Sigma = .10*diag(2)),
                           MASS::mvrnorm(10, mu = c(.20,0), Sigma = .10*diag(2)),
                           MASS::mvrnorm(10, mu = c(.20,.20), Sigma = .10*diag(2))))
  }
  data3[,,2] <- data3[,c(11:30,1:10),2]
  data3[,,3] <- data3[,c(6:30,1:5),3]
  list_dat = list(data1,data2,data3)

  start <- Sys.time()
  MIC(list_dat, 3, 10000) -> test
  Sys.time() - start
  #################################################################################
  ## TS simulation

  # id <- 1
  # newdir <- paste0('dataset',id)
  # cat(newdir,'\n')
  # dir.create(newdir)
  # setwd(newdir)
  ################## dk search ##################
  set.seed(123)
  sim <- MIC_sim(alpha = 0.9,  nsub = 10, fs = 200, segs = 10)
  start <- Sys.time()
  out <- dk_search(sim$Data, max_d = 10, n.iter = 10000, par.spectrum = c(100,50,256), par.win=c(4,2))
  out
  Sys.time() - start
  ################### simple fit ################
  for(i in 1:300){
    set.seed(555)
    ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 200)
    list_data <- lapply(ts_sim$Data, function(x) MIC_prep(X = x, d = 4,
                            par.spec = c(100,50,100), par.win = c(2, 0), unit_len = F))
    # ## plots
    # ep <- list_data[[2]][,,3]; par(mfrow = (c(2,2)));
    # plot((ep[1,]), ep[2,], col = ts_sim$Ci[,2]);plot((ep[1,]), ep[3,], col = ts_sim$Ci[,2])
    # plot((ep[1,]), ep[4,], col = ts_sim$Ci[,2]);plot(ep[2,], ep[3,], col = ts_sim$Ci[,2])

    start <- Sys.time()
    output <- MIC(data = list_data, K = 4, nit = 10000, drop = F); Sys.time() - start

    ## truth and estimated:
    ts_sim$C;align(refL = ts_sim$C, L = output$clustS, type = 'vec')
    t(ts_sim$Ci);for(i in 1:dim(ts_sim$Ci)[2]){
      output$clustC[i,] = align(ts_sim$Ci[,i],output$clustC[i,],type='vec')
    }
    output$clustC
  }
  ################### TS and specs ##############
  ts <- arima.sim(n = 100, model = list(ar = pars(20, fs = 100)))
  ts <- ts - mean(ts)
  spec1 <- Mod(fft(ts))^2/length(ts) ### raw spec
  spec2 <- spec.pgram(ts, plot = F, demean = F, detrend = F)$spec  ### spec.pgram
  spec3 <- spec.parzen(ts, a = 99, nn=50) ### spec.parzen
  matplot(cbind(spec1[2:51], spec2, spec3), type = 'l')
  ################### SpecSim, EigLap and EigLapSph #####
  sim <- MIC_sim(alpha = 1, nsub = 1, Ct = c(1,1,1,2,2,2,3,3,3), fs = 200, segs = 10)$Data[[1]]
  sim <- apply(sim, c(1,3), function(x) lm(x~c(1:200))$residuals)
  test <- SpecSim(sim, 60, 50, 10,0,100)
  for(i in 1:2){
    EigLap(test, 2, (i==2)) -> tplot
    plot(t(tplot[,,1]), col = rep(c(1,2,3),c(3,3,3))) #Consider Cartesan-spherical transform
  }

  specs <- apply(sim, c(2,3), function(x) specParzen(x,100,50,100)/sum(specParzen(x,100,50,100)))
  specs <- apply(specs, c(1,2), mean)
  matplot(specs, type = 'l', col = c(1,1,1,2,2,2,3,3,3))

}
