# Importation of ourside functions
#' @importFrom stats runif rnorm lm arima.sim var acf rexp
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

  id <- 1
  newdir <- paste0('dataset',id)
  cat(newdir,'\n')
  dir.create(newdir)
  setwd(newdir)
  set.seed(1234)
  sim <- MIC_sim(alpha = 0.9,  nsub = 10, Ct = rep(c(1,2,3),rep(30,3)), fs = 200, segs = 10)
  start <- Sys.time()
  search <- dk_search(sim$Data, max_d = 10, n.iter = 10000, par.spectrum = c(50,50,256))
  search
  Sys.time() - start
  ## simpler
  for(i in 1:300){
    set.seed(5)
    ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 200)
    list_data <- lapply(ts_sim$Data, function(x) MIC_prep(X = x, d = 4,
                                                          par.spectrum = c(50, 50, 100), par.win = c(3, 1)))
    start <- Sys.time()
    output <- MIC(data = list_data, K = 4, nit = 10000); Sys.time() - start
    ts_sim$C;as.vector(c(1:4) %*% output$clustS);t(ts_sim$Ci)
    t(apply(output$clustC, c(3), function(x) as.vector(c(1:4) %*% x)))

  }
}
