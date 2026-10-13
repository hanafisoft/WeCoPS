\name{cpca}
\alias{cpca}
\alias{saliences}
\title{Weighted consensus PCA and relative saliences}
\description{
\code{cpca} computes the \eqn{R} dominant eigenvectors of
\eqn{M(w) = \sum_k w_k X_k X_k'}, that is the principal components of the concatenated table
\eqn{[\sqrt{w_1}X_1 | \cdots | \sqrt{w_K}X_K]}. \code{saliences} returns the relative saliences
\eqn{s_{kh} = u_h' X_k X_k' u_h}, the share of the inertia of block \eqn{k} carried by component
\eqn{h}.
}
\usage{
cpca(blocks, weights = NULL, R = 2)
saliences(U, blocks)
}
\arguments{
\item{blocks}{a list of blocks preprocessed with \code{\link{prep_blocks}}.}
\item{weights}{\code{NULL} (default weighting), a weighting name or a vector of weights; see
\code{\link{block_weights}}.}
\item{R}{number of retained components.}
\item{U}{an \eqn{n \times R} matrix of orthonormal components.}
}
\value{
\code{cpca}: a list with \code{U} (principal components), \code{eigenvalues} (all eigenvalues of
\eqn{M(w)}, weights normalised to sum 1), \code{explained} (eigenvalues as a share of their sum)
and \code{weights}. \code{saliences}: a \eqn{K \times R} matrix.
}
\examples{
b <- load_simulated(1)$blocks
a <- cpca(b, R = 5)
round(saliences(a$U, b), 2)
}
