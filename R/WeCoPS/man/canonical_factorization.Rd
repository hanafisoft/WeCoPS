\name{canonical_factorization}
\alias{canonical_factorization}
\title{Canonical factorisation of the blocks}
\description{
Canonical factorisation \eqn{X_k \approx B D_k V_k'} of each block on a basis \eqn{B}:
\eqn{d_{kh} = \|X_k'b_h\|} (so that \eqn{d_{kh}^2} is the salience \eqn{s_{kh}}) and
\eqn{V_k = X_k'BD_k^{-1}} (unit-norm loadings).
}
\usage{
canonical_factorization(blocks, B)
}
\arguments{
\item{blocks}{a list of blocks preprocessed with \code{\link{prep_blocks}}.}
\item{B}{an \eqn{n \times R} orthonormal basis (principal components or common and proper
components).}
}
\value{A list with, for each block, \code{D} (\eqn{R \times R} diagonal matrix) and \code{V}
(\eqn{J_k \times R} loadings).}
\examples{
a <- wecops(load_wines()$blocks, R = 3)
a$factorization$ccp[[1]]$D
}
