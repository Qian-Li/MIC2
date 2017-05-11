#' MIC_sim Helper: pars
#'
#' \code{pars} determines the AR(2) coefficient with desired oscillation properties: peak location, peak width and
#'  sampling frequency.
#'
#' @param eta numeric, peak location
#' @param M   numeric, dispersion of power with narrower peak as M approximates 1
#' @param fs  sampling frequency
#' @return vector of AR(2) coefficients
#' @seealso \code{\link{MIC_sim}} for its usage
#'

pars<- function(eta, M = 1.1, fs){
  phi1 <- - 1 / M ^ 2
  phi2 <- 2 * cos(2 * pi * eta / fs) / M
  return(c(phi2, phi1))
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
#' @param x vector, a univariate time series.
#' @param a integer, max lag to truncate the sample ACF estimates, default \code{length(x)-1} (no truncation).
#'   It has to be smaller than the length fo time seires \code{x}.
#' @param nn integer, resolution of the spectral estimates (number of points on frequency domain).
#'
#' @return List of the following items:
#'   \item{\code{freq}}{frequencies where PSD is evaluated at}
#'   \item{\code{spec}}{estimated correlogram with truncation \code{a} and Parzen window smoothed ACF}
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' spec <- spec.parzen(x, a = 50, nn = 50)
#'
#' plot(spec.pgram(x, plot = "F")$spec, ylab = 'spec')
#' lines(spec$spec)
#' }
#' @export

spec.parzen <- function( x,
                         a  = length(x)-1,
                         nn = 512){
  # Detrrend and Demean:
  x <- lm(x~c(1:length(x)))$residuals
  # Compute the smoothed periodogram using a Parzen window
  if (a >= length(x)){
    warning("The lag a is too big, reset to 'length(x)-1'.")
    a = length(x) - 1
  }
  freq = seq(0,0.5,length.out = nn+1)[-1]
  spec = specParzen(ts = x, lag = a, maxf = floor(length(x)/2), outn = nn)
  return(list(freq = freq, spec = spec))
}
