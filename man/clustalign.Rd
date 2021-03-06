% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{clustalign}
\alias{clustalign}
\title{Cpp function for aligning indicator matrices}
\usage{
clustalign(now, ref)
}
\arguments{
\item{now, }{matrix with its rows to be realigned}

\item{ref, }{matrix to be used as reference}
}
\value{
no value returned, and matrix \code{now} is aligned on the spot
}
\description{
\code{clustalign} performs indicator matrices alignment to achieve maximal
  concordance of \code{now} to the reference \code{ref}, by swapping the rows
  of matrix \code{now}. Please refer to \code{\link{align}} for general vector
  and matrix alignment.
}
