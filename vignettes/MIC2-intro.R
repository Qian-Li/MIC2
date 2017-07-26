## ---- echo = TRUE, eval=FALSE--------------------------------------------
#  devtools::install_github("Qian-Li/MIC2")     ## to install
#  library(MIC2)                                 ## loading

## ---- echo = F, eval = T-------------------------------------------------
# 3*3 Gaussian mixture simulation
  set.seed(12345)
  library(MIC2)
  data1 <- array(dim = c(2,30,3))
  for(ep in 1:3){
    data1[,,ep] <- t(rbind(MASS::mvrnorm(10, mu = c(0,10), Sigma = 3*diag(2)),
                           MASS::mvrnorm(10, mu = c(10,0), Sigma = 5*diag(2)),
                           MASS::mvrnorm(10, mu = c(5,5), Sigma =  3*diag(2))))
  }
  data1[,,2] <- data1[,c(11:30,1:10),2]
  data1[,,3] <- data1[,c(3:30,1,2),3]

  data2 <- array(dim = c(2,30,3))
  for(ep in 1:3){
    data2[,,ep] <- t(rbind(MASS::mvrnorm(10, mu = c(0,5), Sigma = 3*diag(2)),
                           MASS::mvrnorm(10, mu = c(5,0), Sigma = 2*diag(2)),
                           MASS::mvrnorm(10, mu = c(5,5), Sigma = 3*diag(2))))
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

## ---- echo = F, eval = T, fig.height= 6, fig.width= 6--------------------
par(mfrow = c(3,3), oma = c(0,0,0,0))
for(sub in 1:3){
  for(ep in 1:3){
    par(mar = c(1,1,0,0))
    plot(t(list_dat[[sub]][,,ep]), pch = 19, col = rep(1:3,rep(10,3)), xlab = "", ylab = "")
  }
}

## ---- echo = T, eval = T-------------------------------------------------
## -- Model fitting for 3 clusters, 10000 MCMC samples
MIC(list_dat, K = 3, nit = 10000) -> out

## ------------------------------------------------------------------------
## -- Optimal group labels:
out$clustS
## -- individual coherence to the group
out$alpha

## ---- echo = F, eval = T, fig.height= 6, fig.width= 6--------------------
par(mfrow = c(3,3), oma = c(0,0,0,0))
for(sub in 1:3){
  for(ep in 1:3){
    par(mar = c(1,1,0,0))
    plot(t(list_dat[[sub]][,,ep]), pch = out$clustL[[sub]][ep,]+14, 
        col = out$clustC[sub,], xlab = "", ylab = "", cex = 1.3)
  }
}

## ---- echo=T, eval=T-----------------------------------------------------
set.seed(123)
## -- 20 segments of signals sampled at 200Hz frequency, for 5 highly coherent individuals
sim <- MIC_sim(alpha = 0.9,  nsub = 5, fs = 200, segs = 20)

## ---- echo = T, eval = F-------------------------------------------------
#  ## -- optimal model search, with a smoothing setting of window = 4, overlap of 2.(takes approx 3mins)
#  out <- dk_search(sim$Data, max_d = 10, n.iter = 5000, par.win = c(4, 2))
#  ## ## -- search log -- ##
#  ## D = 2: best_K = 3 COH = 0.8568124
#  ## -----------------------------------------
#  ## D = 3: best_K = 4 COH = 0.8476597
#  ## -----------------------------------------
#  ## D = 4: best_K = 4 COH = 0.8485574
#  ## -----------------------------------------

## ---- echo = T, eval = T-------------------------------------------------
## -- Spectral Estimation and Eigen-Laplacian transformation
list_data <- lapply(sim$Data, function(x) MIC_prep(X = x, d = 4,
                                  par.spec = c(80, 50, 100), par.win = c(3, 1)))
## -- Model fitting (takes approx 17secs)
output <- MIC(data = list_data, K = 4, nit = 5000)

## ---- echo = T, eval = T-------------------------------------------------
## -- Realign estimated labels against truth
est_S = align(refL = sim$C, L = output$clustS)
n     = length(sim$Data)
est_C = lapply(1:n, function(i) output$clustC[i,] = align(refL = est_S, L = output$clustC[i,]))
tru_C = lapply(1:n, function(i) sim$Ci[,i])

## ---- echo = T, eval = T-------------------------------------------------
## -- Group label accuracy:
sum(est_S == sim$C) / length(est_S)
## -- Individual labels accuracy:
mapply(function(x,y) sum(align(refL = x, L = y) == x)/length(est_S), tru_C, est_C)

## ---- echo=F, eval = T---------------------------------------------------
## -- a fake label estimate
label <- rep(c(1,2,3), c(32, 48, 44)); colors = c("#4477AA", "#DDCC77", "#CC6677")

## ---- echo=T, eval=T, fig.width=5, fig.height=5--------------------------
par(oma = rep(0,4), mar = rep(0.5,4)); EEGplot(clust = label, color = colors)

