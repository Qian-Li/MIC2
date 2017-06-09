#' 128-EEGNet channel plots
#'
#' \code{EEGplot} plots the clustering results of 124 channels on a 2D projected scalp plot.
#'
#' @param clust vector, cluster labels for 124 channels.
#' @param color vector, colors for K clusters in the plot.
#' @examples
#' \dontrun{
#' # The channel plot:
#'
#' EEGplot(rep(1,124))
#'
#' }
#'
#' @importFrom graphics box plot points
#' @export
EEGplot <- function(clust, color = NULL)
{ # A few checks
  K <- max(clust); if(length(clust)!=124) stop("Length of clust needs to be 124!")
  if(is.null(color)) color <- rep('black', K)
  if(length(color)!=K) stop("Number of colors need to match clusters(max(clust))!")
  # Initiate
  plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(-1,1),ylim=c(-1,1))
  # Box
  box()
  points(x = EEGcoords$x, y = EEGcoords$y, pch = 19, col = color[clust])
}
