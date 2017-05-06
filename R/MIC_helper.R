#' MIC_sim Helper: pars
#'
#' \code{pars} determines the AR(2) coefficient with desired oscillation properties: peak location, peak width and
#'  sampling frequency.
#'
#' @param eta numeric, peak location
#' @param M   numeric greater than 1, narrower as M approximates 1
#' @param fs  sampling frequency
#' @return vector of AR(2) coefficients
#' @seealso \code{\link{MIC_sim}} for its usage
#'
pars<- function(eta, M = 1.1, fs){
  phi1 <- - 1 / M ^ 2
  phi2 <- 2 * cos(2 * pi * eta / fs) / M
  return(c(phi2, phi1))
}

#' Total Variation Distance
#'
#' \code{TVD} returns the Total Variation Distance between two densities.
#'
#' @param w  A vector of sampled points.
#' @param f1 Density 1 evaluated at values w.
#' @param f2 Density 2 evaluated at values w.
#' @return A scalar, in range [0,1] if both f1 and f2 are normalized densities.
#'
#'
TVD <- function(w, f1, f2){
  #Compute TV distance using the trapezoidal rule
  if (length(w) != length(f1) | length(w) != length(f2)) stop("w, f1 and f2 must have the same length")
  n <- length(w)
  int <- (w[2:n]-w[1:(n-1)])%*%(pmin(f1,f2)[2:n]+pmin(f1,f2)[1:(n-1)])/2
  return(as.double(1-int))
}

#' Estimate Spectral Density using Parzen Lag Window
#'
#' \code{spec.parzen} calculate the spectral estimate based on a Fourier transform of a
#'   truncated and Parzen window smoothed auto-covariance function (ACF).
#'
#' The raw periodogram \code{spec.pgram} is not a consistent estimator of the spectral density,
#'   therefore a class of lag window estimators are considered as surrogates in practice achieving
#'   smootheness and consistency.
#'
#' Parzen window estimator works the best when the true spectrum is continuous, and specifically
#'   have peak concentrations at certain frequencies. Such estimators operates on times series
#'   that are presumably zero-mean and stationary, so that demean and detrending are highly
#'   recommended on \code{x} before implementation.
#'
#' @param x A univariate time series.
#' @param a Max lag to truncate the sample ACF estimates, default \code{a=100}. It cannot exceed
#'   the number of observations in \code{x}.
#' @param dt 1/freq, unit sampling interval of \code{x}.
#' @param w0 The minimal frequency of interest.
#' @param wn The maximal frequency of interest, default as freq/2.
#' @param nn Resolution of the spectral estimates (number of points on frequency domain).
#'
#' @return A \code{nn}-by-2 matrix, with the first column being the estimated frequencies,
#'   and the second column being the esimated spectrum.
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' spec <- spec.parzen(x, a = 50, nn = 50)
#'
#'
#' plot(spec, type = 'l')
#' }
#' @export

spec.parzen <- function( x,
                         a = 100,
                         dt = 1/length(x),
                         w0 = 10^(-5),
                         wn = 1/(2*dt),
                         nn = 512){
  # Compute the smoothed periodogram using a Parzen window
  if (a >= (length(x) - 1)) stop("The bandwidth needs to be smaller than TS length")
  ## parzen window smoother
  kp <- numeric(2 * a - 1)
  tt <- seq(-1, 1, length = 2*a - 1)
  kp[abs(tt) < .5] <- 1 - 6*abs(tt[abs(tt) < .5])^2 + 6*(abs(tt[abs(tt) < .5])^3)
  kp[abs(tt) >= .5]<-2*(1 - abs(tt[abs(tt) >= .5]))^3
  ## acf
  Cov <- as.vector(acf(x, lag.max = a-1, type = "covariance", plot=F)$acf) * kp[a:(2*a-1)]
  time <- 2*pi*seq(0, length(x) * dt, by = dt)[1: a]
  w <- seq(w0, wn, length.out = nn)
  SS <- sapply(w, function(x) (4*dt) * sum(cos(x * time)[-1] * Cov[-1]) + (2*dt) * Cov[1])
  ## FFT can double the speed.
  return(SS)
}

#' Epoch-wise pairwise Distance Matrices
#'
#' \code{ep_dismat} calculates strucally smoothed pairwise distances across all objects.
#'   It is an essential prerequisite step of MIC implementation.
#'
#' This procedure consists of distance characterization across all objects to be clustered on
#'   and further smoothing into epochs.
#'
#' The distance metric is based on the \code{\link{TVD}} between a pair of spectral
#'   densities, calculated by the Parzen lag window estimator \code{\link{spec.parzen}}.
#'
#' Epochs are constructed by combining adjacent units (segments), which are further allowed
#'   to overlap with each other at certain amount.
#'
#' @param X 3-d array of input, organized as:
#'   No.objects * No.observations * No.segments. Here objects refer to the clustering units;
#'   observations refer to the recorded time series at each segment; and segments are pieces of recordings
#'   that are repeatedly observed/measured with same length.
#' @param sf integer, sampling frequency
#' @param par.spectrum vector, parameter for \code{\link{spec.parzen}} in the order of \code{c(a, wn, nn)}.
#' @param window vector, epoch-smoothing setting in the order of window size and overlap size.
#'   For example, non-smoothing setting is equivalent to \code{window = c(1, 0)}, and \code{window = c(3, 1)}
#'   combines adjacent \strong{3 segnments as} epochs allowing for \strong{1 segment} overlapping in between.
#' @return a list of objects with the following components:
#'   \item{\code{diss_array}}{3-d array of pairwise distance matrices, No.objects*No.objects*No.epochs}
#'   \item{\code{ave_spec}}{3-d array of epochwise average spectrum see \code{\link{spec.parzen}}}
#'   \item{\code{fw}}{vector, frequencies used for spectral estimate}
#' @seealso \code{\link{MIC_prep}} for usage
#'

ep_dismat <- function(X,
                      sf,
                      par.spectrum = c(100, sf/2, 512),
                      window = c(1,0)){
  if (length(dim(X)) != 3) stop("X needs to be a 3-D array!") else {
    dim(X)[1] -> nt
    dim(X)[2] -> nc
    dim(X)[3] -> ns
  }
  dt<-1/sf
  length.w<-512
  np<-length(par.spectrum)

  if (np==1) {a <- par.spectrum[1]; wn <- sf/2; length.w <- 512}
  if (np==2) {a <- par.spectrum[1]; wn <- par.spectrum[2]; length.w <- 512}
  if (np==3) {a <- par.spectrum[1]; wn <- par.spectrum[2]; length.w <- par.spectrum[3]}

  ## seg-wise normalized spectrum
  Spec <- apply(X, c(2,3), function(channel)
    spec.parzen(channel, a = a, dt = dt, wn = wn, nn = length.w) / var(channel))
  fw <- seq(10e-5, wn, length.out = length.w)

  ## Sliding window averaging.
  if(window[1]==1){
    a.Spec <- Spec
    ne <- ns
    incre = 1
  } else {
    incre <- window[1] - window[2]
    ne <- floor((ns-window[1])/incre)+1
    a.Spec <- array(NA, dim = c(length.w, nc, ne))
    for(k in 1:ne){
      start <- (k - 1) * incre + 1
      end <- start + window[1] - 1
      a.Spec[, , k]<-apply(Spec[, , start:end], c(1, 2), mean)
    }
  }
  rm(Spec)

  ## Epoch distance mat
  MatDiss<-array(NA, c(nc, nc, ne))
  for (k in 1:ne){
    for (i in 1:(nc-1)) {
      for (j in (i+1):nc) {
        MatDiss[i, j, k] <- TVD(w = fw, a.Spec[, i, k], a.Spec[, j, k])
        MatDiss[j, i, k] <- MatDiss[i, j, k]
      }
    }
    diag(MatDiss[, , k]) <- 0
    # Mat <- MatDiss[,,k]
    # Mat [lower.tri(Mat)] <- t(Mat) [lower.tri(Mat)] #Symmetry
    # MatDiss[, , k] <- Mat
    # rm(Mat)
  }
  rm(a.Spec)
  return(MatDiss)
}

#' Eigen-Laplacian Dimension Reduction
#'
#' \code{eigen_lap} carries out eigen-decomposition on the graph Laplacian matrices, which
#'   operates on distance matrices and yields low dimensional representations of original objects.
#'
#' This procedure is a direct implementation of the spectral clustering technique explicated by [Ng 2001].
#' It is also a preprocessing step for \code{\link{MIC}}, following distance characterization by \code{\link{ep_dismat}}.
#'
#' @references Ng, Andrew Y., Michael I. Jordan, and Yair Weiss. \emph{On spectral clustering: Analysis and an algorithm.}
#'   Advances in neural information processing systems 2 (2002): 849-856.
#'
#' @param X 3-d array, suggest \code{diss_array} from \code{\link{ep_dismat}}, symmetric distance matrices stacked into a 3-D array.
#' @param D integer, dimension of the output, corresponding to the number of eigenvectors extracted from graph Laplacian.
#' @return A list of objects with the following components:
#'   \item{\code{eig_data}}{A 3-d array of the eigen-Laplacian representation, No.objects * D * No.epochs}
#'   \item{\code{eig_value}}{matrix, No.epochs * No.objects, complete eigen-values of graph Laplacian matrices}
#'

eigen_lap <- function(X, D){
  nc <- dim(X) [1]
  ne <- dim(X) [3]
  eig_data <- array(NA, c(D, nc, ne))
  eig_val <- matrix(NA, nrow = ne, ncol = nc)
  for (j in 1:ne){
    simMat <- 1 - X[, , j]
    diag(simMat) <- 0
    dia_mat <- diag(1 / sqrt(rowSums(simMat)))
    L <- dia_mat %*% simMat %*% dia_mat
    L_eig <- eigen(L, symmetric = TRUE)
    #vec <- L_eig$vectors
    eig_val [j, ] <- L_eig$values
    eigscore <- L_eig$vectors [, 1:D]
    Ln <- diag(eigscore %*% t(eigscore))
    eig_data[, , j] <- t(eigscore) %*% diag(1 / sqrt(Ln))
    rm(simMat, dia_mat, L, L_eig, eigscore, Ln)
  }
  return(list(eig_data = eig_data, eig_value = eig_val))
}

