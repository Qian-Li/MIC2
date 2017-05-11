#' Data preparation for MIC
#'
#' \code{MIC_prep} prepares stuctured data for \code{MIC}. It sequentially carries out signal detrending, epoch smoothing,
#'   distance characterization and dimensional reduction sequentiall.
#'
#' @param X 3-d array, organized as No.objects * No.observations * No.segments, see simulated example by \code{\link{MIC_sim}}
#' @param d integer, dimensionality of the eigen-Laplacian representation
#' @param exclude vector of integers, indicating objects to be excluded
#' @param spec.lag integer, maximal lag in ACF estimation, see \code{\link{spec.parzen}}
#' @param spec.mfreq integer, maximal frequency to investigate, default at \code{sf/2}
#' @param spec.len integer, length of spectral estimates, see \code{\link{spec.parzen}}
#' @param par.win vector, epoch smoothing parameters see \code{\link{ep_dismat}}
#' @param polar boolean, to use polar coordinates (default \code{F})
#' @return List of data matrices, each with No.objects rows and \code{d} columns.
#' @examples
#' \dontrun{
#' # Simulated data:
#' ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 100)
#'
#' # Data preparation on subject 1
#' sub1 <- MIC_prep(ts_sim$Data[[1]], d = 3, par.win = c(3, 1), polar = T)
#'
#'
#' # Data structure: D / No.channels / No.Epochs
#' dim(sub1)
#'
#'
#' # To visualize preprocessed data on Epoch 1
#' ep1 <- t(sub1[,,1])
#' plot(ep1, col = ts_sim$Ci[,1], xlab = 'dim1', ylab = 'dim2')
#' }
#'@seealso \code{\link{MIC}} for its usage and \code{\link{MIC_sim}} for time series simulation.
#'
#' @export
MIC_prep <- function(X, d,
                     exclude = c(),
                     spec.lag = NULL, spec.mfreq = NULL, spec.len = 512,
                     par.win = c(1, 0), polar = F){
  if (length(exclude) != 0) X <- X[- exclude, , ]         #chanel exclusion
  segs <- dim(X) [3]; fs <- dim(X) [2]; nc <- dim(X) [1]  #extract consts
  # Detrending, might consider Cpp alternative?
  X       <- apply(X, c(1,3), function(x) lm(x ~ c(1:fs))$residuals)
  # Spectral estimation parameters:
  if(is.null(spec.lag)) spec.lag = fs-1
  if(is.null(spec.mfreq)) spec.mfreq = floor(fs/2)
  # Cpp based spectral similarity
  ep_sim  <- SpecSim(ts = X, lag = spec.lag, wn = spec.mfreq, win = par.win[1], overlap = par.win[2], specN = spec.len)
  if(polar){
    ep_eig <- EigLapSph(data = ep_sim, D = d)
  } else{
    ep_eig  <- EigLap(data = ep_sim, D = d, normal = F)
  }
  return(ep_eig)
}
