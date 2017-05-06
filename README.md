# MIC2
Multilevel Integrative Clustering (Rcpp v2), packaged by Qian Li

# Installation:
```r
devtools::install_github("Qian-Li/MIC2")
```

See our manuscript for more details!

# Quick Example

Time series simulation:
```r
library(MIC2)
ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 200)
```

Data preparation:
```r
# Spectral estimation on simulated time series.
list_data <- lapply(ts_sim$Data, function(x) MIC_prep(X = x, d = 4,
                                  par.spectrum = c(50, 50, 100), par.win = c(3, 1)))
```


MIC modle fitting:
```r
# Running time approx.?
start <- Sys.time()
output <- MIC(data = list_data, K = 4, nit = 10000); 
Sys.time() - start
```

Clustering result:
```r
# (DONT RUN!) Accuracy relative to the true group labels
sum(clust_align(ts_sim$C, output$Cbest, type = 'vec') == ts_sim$C) / 40
```
