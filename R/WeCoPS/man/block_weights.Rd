\name{block_weights}
\alias{block_weights}
\title{Usual block weightings, relative to the default weighting}
\description{
Weights of the usual weightings, expressed relative to the default weighting (blocks of unit
total inertia). Only ratios between weights matter.
}
\usage{
block_weights(blocks, method = "default")
}
\arguments{
\item{blocks}{a list of blocks preprocessed with \code{\link{prep_blocks}}.}
\item{method}{\code{"default"} (all weights 1); \code{"mfa"} (inverse of the first eigenvalue of
each block, as in multiple factor analysis); \code{"rv"} (inverse of
\eqn{\|X_kX_k'\|_F}, the RV normalisation of STATIS); \code{"statis"} (RV normalisation times the
weights of the STATIS compromise); \code{"nvar"} (inverse of the number of variables); or a
vector of \eqn{K} positive weights, returned unchanged.}
}
\value{A numeric vector of length \eqn{K}.}
\examples{
b <- load_wines()$blocks
sapply(c("default", "rv", "mfa", "statis", "nvar"), function(m) block_weights(b, m))
}
