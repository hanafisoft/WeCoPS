\name{wecops}
\alias{wecops}
\alias{print.wecops}
\title{Weighted CPCA read in two bases}
\description{
Weighted consensus PCA with its principal components and its common and proper components:
saliences, labels and motifs in both bases, mixing of the principal components, four-part
decomposition of each block, and canonical factorisations \eqn{X_k \approx B D_k V_k'} on both
bases.
}
\usage{
wecops(blocks, R, weights = "default", alpha = 0.05, nstart = 10, seed = 0,
       block_names = NULL)
\method{print}{wecops}(x, ...)
}
\arguments{
\item{blocks}{a list of blocks preprocessed with \code{\link{prep_blocks}}.}
\item{R}{number of retained components.}
\item{weights}{a weighting name or a vector of weights; see \code{\link{block_weights}}.}
\item{alpha}{salience threshold of the labels.}
\item{nstart, seed}{starts of the concentration; see \code{\link{ccp}}.}
\item{block_names}{optional names of the blocks.}
\item{x}{a \code{wecops} object.}
\item{...}{unused.}
}
\value{
An object of class \code{wecops}: a list with \code{U}, \code{saliences_pc}, \code{labels_pc},
\code{motif_pc}, \code{V}, \code{saliences}, \code{T}, \code{labels}, \code{motif},
\code{mixing}, \code{parts}, \code{factorization} (\code{ccp} and \code{pc}), \code{eigenvalues}
and \code{weights}.
}
\examples{
a <- wecops(load_simulated(1)$blocks, R = 5)
a
a$motif_pc
a$motif
}
